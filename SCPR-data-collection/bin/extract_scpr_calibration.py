#!/usr/bin/env python3
"""Extract CSCPR calibration values from downloaded physical-unit .dat files.

Offline, Python standard library only. Independently processes both gauges:
minimum after the downward threshold crossing -> 12 consecutive 100-s bins
of pressure minus the nearest internal barometer -> one-sigma trimmed means.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys
import tempfile


DEFAULT_INPUT_ROOT = Path(__file__).resolve().parent.parent / "sample_data"   # SCPR-data-collection/sample_data
KPA_PER_PSI = 6.894757293168361
BIN_SECONDS = 100
N_BINS = 12
DURATION = BIN_SECONDS * N_BINS
MIN_BIN_SAMPLES = 90  # Nominal 1 Hz; tolerate timestamp jitter and small losses.
MAX_GAUGE_GAP_SECONDS = 5.0
SAMPLE_PADDING_SECONDS = 60
RETURN_HOLD_SECONDS = 30
SAMPLE_BAROMETER_MAX_LAG_SECONDS = 1.0  # Strictly less than 1 s; otherwise leave blank.


@dataclass(frozen=True)
class Sample:
    time: datetime
    pressure: float
    temperature: float


def utc_time(value: str) -> datetime:
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError(f"Timestamp must specify UTC or its offset: {value}")
    return result.astimezone(timezone.utc)


def iso(value: datetime) -> str:
    return value.isoformat(timespec="microseconds").replace("+00:00", "Z")


def round_second(value: datetime) -> datetime:
    """Round UTC timestamps to the nearest second, with half-second ties up."""
    return (value + timedelta(microseconds=500_000)).replace(microsecond=0)


def load_samples(path: Path, *, auxiliary: bool = False) -> list[Sample]:
    """Read current three-column gauge or nine-column auxiliary .dat schema."""
    samples = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            try:
                if len(fields) != (9 if auxiliary else 3):
                    raise ValueError("unexpected column count; use the physical-unit .dat files")
                pressure = float(fields[2 if auxiliary else 1])
                temperature = float(fields[3 if auxiliary else 2])
                if not (math.isfinite(pressure) and math.isfinite(temperature)):
                    raise ValueError("non-finite pressure or temperature")
                samples.append(Sample(utc_time(fields[0]), pressure, temperature))
            except ValueError as exc:
                raise ValueError(f"{path.name}:{number}: {exc}") from exc
    if not samples:
        raise ValueError(f"No measurements in {path}")
    # Stable sort retains all records, including coincident logger timestamps.
    return sorted(samples, key=lambda sample: sample.time)


def filter_barometer_startup(samples: list[Sample], raw_paths: list[Path],
                             settling_seconds: float) -> tuple[list[Sample], dict]:
    """Exclude the PTB110 settling interval after a logged SensePowerOn call.

    Both times are logger arrival times, so this is a conservative operational
    quality rule, not proof of the exact hardware power-on instant. Source files
    remain unchanged; the exclusion applies to barometer selection only.
    """
    audit = {
        "settling_seconds": settling_seconds,
        "rule": "Exclude auxiliary observations with 0 <= time - logged SensePowerOn < settling_seconds",
        "specification": "https://docs.vaisala.com/api/khub/documents/gjV4l~_i5iS2AELD6NvlVQ/content",
        "raw_logs": [], "power_on_times_utc": [], "excluded_samples": [],
    }
    if settling_seconds == 0:
        return samples, audit
    if not raw_paths:
        raise ValueError("Missing raw SCPR log for barometer startup checks. Keep the downloaded raw/ folder, "
                         "or explicitly use --barometer-settling-seconds 0 to disable this check.")
    starts = set()
    for path in raw_paths:
        raw = path.read_bytes()
        audit["raw_logs"].append({"path": str(path), "sha256": hashlib.sha256(raw).hexdigest()})
        for line in raw.decode("utf-8", errors="replace").splitlines():
            match = re.match(r"^<log-ts:([^>]+)>.*\bSensePowerOn\b", line)
            if match:
                value = match.group(1)
                timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
                if timestamp.tzinfo is None:
                    timestamp = timestamp.replace(tzinfo=timezone.utc)
                starts.add(timestamp.astimezone(timezone.utc))
    starts = sorted(starts)
    audit["power_on_times_utc"] = [iso(t) for t in starts]
    kept = []
    for sample in samples:
        index = bisect_right(starts, sample.time) - 1
        elapsed = (sample.time - starts[index]).total_seconds() if index >= 0 else math.inf
        if 0 <= elapsed < settling_seconds:
            audit["excluded_samples"].append({"time_utc": iso(sample.time), "barometer_psi": sample.pressure,
                                              "power_on_utc": iso(starts[index]), "elapsed_seconds": elapsed})
        else:
            kept.append(sample)
    if not kept:
        raise ValueError("No barometer observations remain after startup checks")
    return kept, audit


def detect_events(samples: list[Sample], threshold: float) -> tuple[list[int], list[dict]]:
    """Match the paper's threshold/minimum event definition, with coverage checks."""
    times = [sample.time for sample in samples]
    accepted: list[int] = []
    rejected: list[dict] = []
    for i in range(1, len(samples)):
        if not (samples[i - 1].pressure >= threshold > samples[i].pressure):
            continue
        crossing = times[i]
        if accepted and crossing < times[accepted[-1]] + timedelta(seconds=DURATION):
            continue
        search_end = bisect_right(times, crossing + timedelta(seconds=100))
        origin = min(range(i, search_end), key=lambda j: samples[j].pressure)
        t0 = times[origin]
        pre_time = t0 - timedelta(minutes=5)
        pre = bisect_left(times, pre_time)
        low_end = bisect_right(times, t0 + timedelta(minutes=10))
        end = bisect_left(times, t0 + timedelta(seconds=DURATION))
        reason = ""
        if times[0] > pre_time or (times[pre] - pre_time).total_seconds() > 5:
            reason = "missing pre-calibration coverage at t0 - 300 s"
        elif samples[pre].pressure <= threshold:
            reason = "pressure was already below the threshold at t0 - 300 s"
        elif times[-1] < t0 + timedelta(seconds=DURATION - 1):
            reason = "file ends before the complete 20-minute calibration interval"
        elif max(s.pressure for s in samples[origin:low_end]) > threshold:
            reason = "pressure does not remain low for the first 10 minutes"
        elif any(s.pressure >= threshold for s in samples[low_end:end]):
            reason = "pressure returns to ambient before the 20-minute interval ends"
        if reason:
            rejected.append({"crossing_utc": iso(crossing), "t0_utc": iso(t0), "reason": reason})
        else:
            accepted.append(origin)
    return accepted, rejected


def nearest_sample(time: datetime, samples: list[Sample], times: list[datetime]) -> Sample:
    j = bisect_left(times, time)
    if j == 0:
        return samples[0]
    if j == len(times) or time - times[j - 1] <= times[j] - time:
        return samples[j - 1]  # Earlier observation wins an equal-distance tie.
    return samples[j]


def summarize_event(
    channel: str, event_in_day: int, samples: list[Sample], origin: int,
    auxiliary: list[Sample], max_baro_gap: float | None,
) -> tuple[list[list], list[list], dict]:
    t0 = samples[origin].time
    end = t0 + timedelta(seconds=DURATION)
    aux_times = [sample.time for sample in auxiliary]
    rows: list[list] = []
    first = bisect_left([sample.time for sample in samples], t0)
    event = []
    for sample in samples[first:]:
        if sample.time >= end:
            break
        event.append(sample)
    gaps = [(b.time - a.time).total_seconds() for a, b in zip(event, event[1:])]
    gaps.append((end - event[-1].time).total_seconds())
    if max(gaps) > MAX_GAUGE_GAP_SECONDS:
        raise ValueError(f"{channel} at {iso(t0)}: gauge data gap exceeds {MAX_GAUGE_GAP_SECONDS:g} s")
    for sample in event:
        baro = nearest_sample(sample.time, auxiliary, aux_times)
        dt = (baro.time - sample.time).total_seconds()
        if max_baro_gap is not None and abs(dt) > max_baro_gap:
            raise ValueError(f"{channel} at {iso(sample.time)}: nearest barometer is {abs(dt):.3f} s away")
        elapsed = (sample.time - t0).total_seconds()
        rows.append([iso(sample.time), elapsed, sample.pressure, sample.temperature,
                     iso(baro.time), baro.pressure, dt, sample.pressure - baro.pressure,
                     int(elapsed // BIN_SECONDS) + 1, 0])
    bins = []
    for number in range(1, N_BINS + 1):
        positions = [i for i, row in enumerate(rows) if row[8] == number]
        if len(positions) < MIN_BIN_SAMPLES:
            raise ValueError(f"{channel} at {iso(t0)}: bin {number} has fewer than {MIN_BIN_SAMPLES} native samples")
        values = [rows[i][7] for i in positions]
        mean = statistics.fmean(values)
        sigma = statistics.pstdev(values)  # Population SD, ddof=0, one pass.
        kept = [i for i in positions if sigma == 0 or abs(rows[i][7] - mean) <= sigma]
        if not kept:
            raise ValueError(f"{channel} at {iso(t0)}: no retained samples in bin {number}")
        for i in kept:
            rows[i][9] = 1
        bins.append([number, (number - 1) * BIN_SECONDS, number * BIN_SECONDS,
                     len(values), len(kept), statistics.fmean(rows[i][7] for i in kept), mean, sigma])
    means = [row[5] for row in bins]
    sd = statistics.stdev(means)  # Sample SD across the 12 means, ddof=1.
    summary = {
        "channel": channel, "event_in_day": event_in_day, "t0_utc": iso(t0),
        "native_samples": len(rows), "retained_samples": sum(row[4] for row in bins),
        "bins_available": len(bins), "calibration_pressure_psi": statistics.fmean(means),
        "bin_means_sd_psi": sd, "sasagawa_se_psi": sd / math.sqrt(N_BINS - 1),
        "max_abs_barometer_dt_s": max(abs(row[6]) for row in rows),
        "barometer_pairs_over_2s": sum(abs(row[6]) > 2 for row in rows),
    }
    return rows, bins, summary


def write_dat(path: Path, comments: list[str], columns: list[str], rows: list[list]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for comment in comments:
            handle.write("# " + comment + "\n")
        handle.write("# " + "\t".join(columns) + "\n")
        for row in rows:
            fields = ["" if v is None or (isinstance(v, float) and not math.isfinite(v))
                      else format(v, ".17g") if isinstance(v, float) else str(v) for v in row]
            handle.write("\t".join(fields) + "\n")


def paired_events(summaries: list[dict]) -> list[tuple[dict, dict]]:
    p1 = sorted((s for s in summaries if s["channel"] == "PARO1"), key=lambda s: s["t0_utc"])
    p2 = sorted((s for s in summaries if s["channel"] == "PARO2"), key=lambda s: s["t0_utc"])
    if len(p1) != len(p2):
        raise ValueError("Different numbers of PARO1/PARO2 events; cannot write paired summary")
    for first, second in zip(p1, p2):
        if abs((utc_time(second["t0_utc"]) - utc_time(first["t0_utc"])).total_seconds()) > 100:
            raise ValueError("PARO1/PARO2 calibration minima are more than 100 s apart")
    return list(zip(p1, p2))


def compact_summary(summaries: list[dict], auxiliary: list[Sample]) -> list[dict]:
    """One row per paired event; PARO1 defines the shared time and barometer window."""
    rows = []
    for first, second in paired_events(summaries):
        start = utc_time(first["t0_utc"])
        end = start + timedelta(seconds=DURATION)
        barometer = [sample.pressure for sample in auxiliary if start <= sample.time < end]
        if not barometer:
            raise ValueError(f"No internal barometer samples in the interval starting {iso(start)}")
        center = start + timedelta(seconds=DURATION / 2)
        rounded = round_second(center)
        rows.append({
            "time": rounded.isoformat(timespec="seconds").replace("+00:00", "Z"),
            "cal_Paro1": first["calibration_pressure_psi"],
            "cal_Paro2": second["calibration_pressure_psi"],
            "barom": statistics.fmean(barometer),
        })
    return rows


def external_pressure_return(samples: list[Sample], t0: datetime, threshold: float) -> datetime:
    """First upward crossing after the 1200-s interval, sustained for 30 s."""
    times = [sample.time for sample in samples]
    start = max(1, bisect_left(times, t0 + timedelta(seconds=DURATION)))
    for i in range(start, len(samples)):
        if not (samples[i - 1].pressure < threshold <= samples[i].pressure):
            continue
        stop = bisect_left(times, times[i] + timedelta(seconds=RETURN_HOLD_SECONDS))
        if stop == len(samples):
            break
        if any(s.pressure < threshold for s in samples[i:stop + 1]):
            continue
        if any((times[j] - times[j - 1]).total_seconds() > MAX_GAUGE_GAP_SECONDS
               for j in range(i, stop + 1)):
            continue
        return times[i]
    raise ValueError(f"No sustained return to external pressure after {iso(t0)}")


def rounded_stream(samples: list[Sample], start: datetime, end: datetime) -> tuple[dict, dict]:
    """Use native readings on a rounded 1-s grid; neither average nor fill gaps."""
    selected: dict[datetime, Sample] = {}
    native_count = 0
    for sample in samples:
        tick = round_second(sample.time)
        if tick < start:
            continue
        if tick > end:
            break
        native_count += 1
        previous = selected.get(tick)
        # Stable input order handles coincident timestamps; earlier time wins ties.
        if previous is None or (abs(sample.time - tick), sample.time) < (abs(previous.time - tick), previous.time):
            selected[tick] = sample
    n_seconds = int((end - start).total_seconds()) + 1
    return selected, {
        "native_samples": native_count, "occupied_seconds": len(selected),
        "rounding_collisions": native_count - len(selected),
        "missing_seconds": n_seconds - len(selected),
    }


def combined_samples(summaries: list[dict], gauges: dict[str, list[Sample]],
                     auxiliary: list[Sample], threshold: float) -> tuple[list[list], list[dict]]:
    """Both gauges and barometer from 60 s before minima to 60 s after return."""
    output = {}
    windows = []
    aux_times = [sample.time for sample in auxiliary]
    for number, pair in enumerate(paired_events(summaries), 1):
        minima = {s["channel"]: utc_time(s["t0_utc"]) for s in pair}
        returns = {channel: external_pressure_return(gauges[channel], t0, threshold)
                   for channel, t0 in minima.items()}
        exact_start = min(minima.values()) - timedelta(seconds=SAMPLE_PADDING_SECONDS)
        exact_end = max(returns.values()) + timedelta(seconds=SAMPLE_PADDING_SECONDS)
        start, end = round_second(exact_start), round_second(exact_end)
        streams, counts = {}, {}
        for channel, samples in gauges.items():
            if round_second(samples[0].time) > start or round_second(samples[-1].time) < end:
                raise ValueError(f"{channel}: input does not cover the full padded interval {iso(start)} to {iso(end)}")
            streams[channel], counts[channel] = rounded_stream(samples, start, end)
        barometer_lags, matched_lags, barometer_times = [], [], set()
        for offset in range(int((end - start).total_seconds()) + 1):
            tick = start + timedelta(seconds=offset)
            p1, p2 = (streams[ch].get(tick) for ch in ("PARO1", "PARO2"))
            baro = nearest_sample(tick, auxiliary, aux_times)
            lag = (baro.time - tick).total_seconds()
            barometer_lags.append(lag)
            matched = abs(lag) < SAMPLE_BAROMETER_MAX_LAG_SECONDS
            if matched:
                matched_lags.append(lag)
                barometer_times.add(baro.time)
            output[tick] = [tick.isoformat(timespec="seconds").replace("+00:00", "Z"),
                            p1.pressure if p1 else math.nan, p2.pressure if p2 else math.nan,
                            p1.temperature if p1 else math.nan, p2.temperature if p2 else math.nan,
                            baro.pressure if matched else None]
        windows.append({
            "event_in_day": number, "t0_utc": {ch: iso(t) for ch, t in minima.items()},
            "external_pressure_return_utc": {ch: iso(t) for ch, t in returns.items()},
            "unrounded_start_utc": iso(exact_start), "unrounded_end_utc": iso(exact_end),
            "start_utc": start.isoformat(timespec="seconds").replace("+00:00", "Z"),
            "end_utc": end.isoformat(timespec="seconds").replace("+00:00", "Z"),
            "row_count": int((end - start).total_seconds()) + 1, "streams": counts,
            "barometer_matching": {
                "method": "nearest eligible native auxiliary timestamp; ties earlier; blank if abs(lag) >= 1 s",
                "max_lag_seconds_exclusive": SAMPLE_BAROMETER_MAX_LAG_SECONDS,
                "matched_rows": len(matched_lags), "blank_rows": len(barometer_lags) - len(matched_lags),
                "unique_source_timestamps": len(barometer_times),
                "max_abs_lag_seconds": max((abs(lag) for lag in matched_lags), default=None),
                "max_candidate_abs_lag_seconds": max(abs(lag) for lag in barometer_lags),
                "first_source_time_utc": iso(min(barometer_times)) if barometer_times else None,
                "last_source_time_utc": iso(max(barometer_times)) if barometer_times else None,
            },
        })
    return [output[tick] for tick in sorted(output)], windows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", type=date.fromisoformat, default=date(2020, 5, 4))
    parser.add_argument("--input-dir", type=Path, help="Day folder containing the physical .dat files")
    parser.add_argument("--out-dir", type=Path, help="Default: INPUT_DIR/calibration")
    parser.add_argument("--threshold-psi", type=float, default=15400.0 / KPA_PER_PSI,
                        help="Downward-crossing threshold (default: 15400 kPa expressed in psi)")
    parser.add_argument("--max-baro-gap-seconds", type=float,
                        help="Optional calibration-statistics lag limit; fail if exceeded. Default: paper's uncapped nearest match")
    parser.add_argument("--barometer-settling-seconds", type=float, default=1.0,
                        help="Exclude barometer readings this many seconds after logged SensePowerOn (default: 1; 0 disables)")
    args = parser.parse_args()
    if not math.isfinite(args.threshold_psi):
        parser.error("--threshold-psi must be finite")
    if not math.isfinite(args.barometer_settling_seconds) or args.barometer_settling_seconds < 0:
        parser.error("--barometer-settling-seconds must be finite and nonnegative")
    if args.max_baro_gap_seconds is not None and (
        not math.isfinite(args.max_baro_gap_seconds) or args.max_baro_gap_seconds < 0
    ):
        parser.error("--max-baro-gap-seconds must be finite and nonnegative")
    input_dir = (args.input_dir or DEFAULT_INPUT_ROOT / str(args.date)).expanduser().resolve()
    out_dir = (args.out_dir or input_dir / "calibration").expanduser().resolve()
    source_paths = {channel: input_dir / f"SCPR_{channel}_{args.date:%Y%m%d}_physical.dat"
                    for channel in ("PARO1", "PARO2", "SCPR")}
    auxiliary = load_samples(source_paths["SCPR"], auxiliary=True)
    raw_logs = sorted((input_dir / "raw").glob(f"SCPRAA301_SCPR_*_{args.date:%Y%m%d}T0000_UTC.dat"))
    auxiliary, barometer_quality = filter_barometer_startup(auxiliary, raw_logs, args.barometer_settling_seconds)
    artifacts = []
    rejected = {}
    gauges = {}
    for channel in ("PARO1", "PARO2"):
        samples = load_samples(source_paths[channel])
        gauges[channel] = samples
        if any(sample.time.date() != args.date for sample in samples):
            raise ValueError(f"{channel}: input has samples outside the requested UTC day")
        origins, rejected[channel] = detect_events(samples, args.threshold_psi)
        if not origins:
            detail = "; ".join(record["reason"] for record in rejected[channel])
            raise ValueError(f"No complete calibration found for {channel} on {args.date}. {detail}")
        for number, origin in enumerate(origins, 1):
            artifacts.append(summarize_event(channel, number, samples, origin, auxiliary, args.max_baro_gap_seconds))
    summaries = [artifact[2] for artifact in artifacts]
    summary_rows = compact_summary(summaries, auxiliary)
    sample_rows, sample_windows = combined_samples(summaries, gauges, auxiliary, args.threshold_psi)
    sample_comments = [
        "UTC; absolute pressures and internal barometer in psi; temperatures in degC.",
        "Window: earliest gauge minimum - 60 s through latest gauge external-pressure return + 60 s (inclusive).",
        f"Return: upward crossing of {args.threshold_psi:.17g} psi, sustained for {RETURN_HOLD_SECONDS} s; maximum gauge gap 5 s.",
        "Gauge times and window bounds rounded to nearest second (half up) before selection on a common 1-s grid.",
        "Gauge same-second collisions: native reading closest to the whole second; ties use earlier time, then input order.",
        "Barometer: nearest native reading from the full input day to each output second; ties use earlier time.",
        "Barometer is blank if the nearest eligible observation is 1 s or more away from the output time.",
        f"Barometer quality: exclude the first {args.barometer_settling_seconds:g} s after logged SensePowerOn; see calibration_manifest.json.",
        "No interpolation, averaging or barometer subtraction. Missing numeric values are empty tab-separated fields.",
    ]
    summary_comments = [
        "UTC; all pressures in psi. cal_Paro1/cal_Paro2 are pressure minus internal barometer.",
        "time: PARO1 interval midpoint (t0+600 s), rounded to nearest second (half up).",
        "barom: arithmetic mean of native internal barometer samples in PARO1 [t0,t0+1200 s).",
    ]
    comments = [
        "CSCPR calibration: physical-unit input; no downloads; pressures in psi, temperatures in degC.",
        "t0: minimum within 100 s after downward crossing; gauge-specific UTC arrival times.",
        "Interval [t0, t0+1200 s); 12 non-overlapping 100-s bins, native unsmoothed samples.",
        "Gauge coverage: at least 90 native samples per 100-s bin; no gap longer than 5 s.",
        "Subtract nearest internal barometer; ties use earlier sample. No interpolation or default lag cutoff.",
        "barometer_dt_s = barometer_time_utc - time_utc (signed seconds).",
        "One-pass trimming per bin: retain |P-baro - bin_mean| <= population SD (ddof=0).",
        "Calibration pressure = equal-weight mean of all 12 trimmed bin means; no DWT-reference subtraction.",
        "bin_means_sd: sample SD (ddof=1); sasagawa_se = bin_means_sd/sqrt(11).",
        "event_in_day is local detection order, not the FAQ or paper event number.",
    ]
    provenance = {
        "date_utc": str(args.date), "method": comments, "threshold_psi": args.threshold_psi,
        "max_barometer_gap_seconds": args.max_baro_gap_seconds,
        "sources": {channel: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for channel, path in source_paths.items()},
        "events": summaries, "rejected_candidates": rejected,
        "summary_definition": summary_comments, "summary_rows": summary_rows,
        "sample_definition": sample_comments, "sample_windows": sample_windows,
        "barometer_quality": barometer_quality,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    # Compute and validate both channels before replacing any existing result.
    with tempfile.TemporaryDirectory(prefix=".calibration-", dir=out_dir) as temporary:
        staging = Path(temporary)
        columns = ["time", "cal_Paro1", "cal_Paro2", "barom"]
        write_dat(staging / "SCPR_calibration_summary.dat", summary_comments,
                  columns, [[summary[key] for key in columns] for summary in summary_rows])
        write_dat(staging / "SCPR_calibration_samples.dat", sample_comments,
                  ["time", "pressure1_psi", "pressure2_psi", "temperature1_degC", "temperature2_degC", "barom_psi"],
                  sample_rows)
        (staging / "calibration_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        for path in staging.iterdir():
            path.replace(out_dir / path.name)
    for summary in summaries:
        print(f"{summary['channel']} t0={summary['t0_utc']}: "
              f"{summary['calibration_pressure_psi']:.12f} psi; "
              f"{summary['retained_samples']}/{summary['native_samples']} samples retained; "
              f"max |barometer lag|={summary['max_abs_barometer_dt_s']:.3f} s")
    print(f"Saved: {out_dir}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
