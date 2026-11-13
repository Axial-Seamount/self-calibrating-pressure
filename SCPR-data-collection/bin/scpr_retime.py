"""Reconstruct sample times for SCPR records that piweb logged without timestamps.

Some periods were logged with no <log-ts:...> arrival stamps (2024-04-11 19:22 to
2024-04-17 16:19 and from 2026-10-17 20:56 UTC). The instrument kept recording: the
gauges write one line per second and every auxiliary (SCPR) record carries the
instrument clock (UNIX_SEC), which ran without resets through both gaps.

Method (chosen and validated on 179 timestamped days by stripping the stamps and
comparing with the converter's normal times; see README):
  * instrument-clock offset (clock - UTC arrival): a straight line between clock resets,
    so it is predicted from the daily offsets in the surrounding days' manifests (both
    sides of a gap when available). Prediction error: median 0.01 s, 90% within 0.06 s;
  * SCPR (auxiliary) records: time = instrument clock - predicted offset;
  * gauge records: one per second from the start of a midnight file (+0.5 s); outages
    taken from the instrument clock (the seconds missing there); whole seconds corrected
    from the clock only when it disagrees by > 1.5 s (file switched off midnight, e.g.
    2026-10-19). Runs in restart files (e.g. _T2056) are anchored on the clock instead;
  * acceptance: gauge and SCPR runs must agree (line counts within MAX_COUNT_DIFF,
    clock gaps <= MAX_CLOCK_GAP_S); otherwise the run is reported, not reconstructed.
  Validation: gauges 97.6% of gauge-days with every sample within 0.75 s (98.8% with
  >= 99% of samples within 0.75 s; worst 22 s); SCPR 95.3% of days fully within 0.75 s.
Reconstructed days are labelled in the manifest ("timing_reconstructed").
"""

from __future__ import annotations

from datetime import date
import json
from pathlib import Path
import re

import numpy as np

MIDNIGHT_OFFSET_S = 0.5          # first record of a midnight file arrives within its first second
MAX_COUNT_DIFF = 1               # gauge lines vs SCPR records in the same run
MAX_CLOCK_GAP_S = 5              # total seconds missing from the instrument clock in the run
CLOCK_OVERRIDE_S = 1.5           # clock must disagree with the midnight anchor by more than this
RESET_JUMP_S = 3.0               # a day-to-day offset change larger than this is a clock reset
VALIDATION = ("179 timestamped days, stamps removed: gauges 97.6% of gauge-days with every sample within "
              "0.75 s of the normal timing (98.8% with >= 99% of samples; worst 22 s); SCPR channels 95.3% of "
              "days fully within 0.75 s; clock-offset prediction median error 0.01 s")

STAMPED = re.compile(r"^<log-ts:")
GAUGE_RECORD = re.compile(r"^\*\d{4},")
SCPR_RECORD = re.compile(r"^1 \d{9,} ")
FILE_TIME = re.compile(r"_\d{8}T(\d{4})_UTC\.dat$")


def offset_history(manifest_dir: Path) -> list[tuple[int, float]]:
    """(date ordinal, daily median clock offset) from manifests of days with real timestamps."""
    out = []
    for path in sorted(manifest_dir.glob("scpr_mseed_*.json")):
        try:
            m = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        scpr = m.get("streams", {}).get("SCPR", {})
        if scpr.get("clock_offset_s") and not scpr.get("reconstructed_records"):
            out.append((date.fromisoformat(m["date_utc"]).toordinal(), float(scpr["clock_offset_s"]["median"])))
    return out


def _side_fit(days: np.ndarray, vals: np.ndarray, before: bool) -> tuple[np.ndarray, np.ndarray]:
    """Keep only the part of one side that is on the same clock as the gap (no reset)."""
    cut = np.flatnonzero(np.abs(np.diff(vals)) > RESET_JUMP_S)
    if len(cut):
        days, vals = (days[cut[-1] + 1:], vals[cut[-1] + 1:]) if before else (days[:cut[0] + 1], vals[:cut[0] + 1])
    return days, vals


def predict_offset(history: list[tuple[int, float]], day: date, n_days: int = 30) -> tuple[float | None, str]:
    """Clock offset (s) for `day` from straight-line fits to the nearest timestamped days
    before and after it (up to n_days each side, however far away; the clock drifts in a
    straight line between resets, e.g. -0.555 s/day in 2026)."""
    if not history:
        return None, "no offset history"
    d = np.array([h[0] for h in history], float)
    v = np.array([h[1] for h in history])
    o = day.toordinal()
    preds = []
    for before in (True, False):
        idx = np.flatnonzero(d < o)[-n_days:] if before else np.flatnonzero(d > o)[:n_days]
        dd, vv = _side_fit(d[idx], v[idx], before)
        if len(dd) >= 3:
            b, a = np.polyfit(dd - o, vv, 1)
            preds.append((float(a), float(np.min(np.abs(dd - o))), "before" if before else "after"))
    if not preds:
        return None, "fewer than 3 timestamped days on the same clock"
    if len(preds) == 2 and abs(preds[0][0] - preds[1][0]) <= 2.0:
        w = np.array([1.0 / p[1] for p in preds])
        return (float(np.dot(w, [p[0] for p in preds]) / w.sum()),
                f"interpolated between timestamped days {preds[0][1]:.0f} d before and {preds[1][1]:.0f} d after")
    return preds[0][0], f"extrapolated from timestamped days {preds[0][2]} (nearest {preds[0][1]:.0f} d away)"


def unstamped_runs(paths: list[Path], record: re.Pattern, markers: tuple[str, ...] = ()) -> list[dict]:
    """Maximal runs of untimestamped record lines, per file, in order.

    Each run: file name, whether it starts at the first record of a midnight file,
    its record payloads, and (index, text) of marker lines found inside it.
    """
    runs = []
    for path in paths:
        midnight = (m := FILE_TIME.search(path.name)) is not None and m[1] == "0000"
        current, seen_record = None, False
        with path.open("rb") as handle:
            for raw in handle:
                text = raw.decode("ascii", "replace").rstrip("\r\n").replace("\0", "")
                if STAMPED.match(text):
                    seen_record = True
                    current = None
                    continue
                if record.match(text):
                    if current is None:
                        current = {"file": path.name, "at_midnight_start": midnight and not seen_record,
                                   "payloads": [], "markers": []}
                        runs.append(current)
                    current["payloads"].append(text)
                    seen_record = True
                elif current is not None and any(k in text for k in markers):
                    current["markers"].append((len(current["payloads"]) - 1, text))
    return runs


def reconstruct(gauge_runs: dict[str, list[dict]], scpr_runs: list[dict], day0: float,
                stamped_scpr: tuple[np.ndarray, np.ndarray],
                predicted_offset: float | None = None) -> tuple[dict, list[dict]]:
    """Times for untimestamped runs.

    gauge_runs: {"PARO1": runs, "PARO2": runs}; scpr_runs: SCPR runs (payloads "1 <unix> ...");
    stamped_scpr: (arrival times, instrument clock) of the day's timestamped SCPR records,
    used to anchor runs that do not start a midnight file.
    Returns ({channel: [(times, payloads)]}, report).
    """
    out: dict[str, list] = {"PARO1": [], "PARO2": [], "SCPR": []}
    report = []
    for k, srun in enumerate(scpr_runs):
        unix = np.array([int(p.split()[1]) for p in srun["payloads"]], dtype=float)
        rel = unix - unix[0]                                   # seconds since the run's first record
        clock_gap = int(rel[-1] - (len(rel) - 1)) if len(rel) else 0
        missing = rel - np.arange(len(rel))                    # seconds missing before each record
        entry = {"run": k, "file": srun["file"], "scpr_records": len(unix), "clock_gap_s": clock_gap}
        t_st, u_st = stamped_scpr
        if len(t_st) and not srun["at_midnight_start"]:
            offset = float(np.median((u_st - t_st)[-3600:]))   # clock - arrival, last hour of stamps
            source = f"{offset:+.2f} s measured from that day's stamps"
        elif predicted_offset is not None:
            offset, source = predicted_offset, f"{predicted_offset:+.2f} s predicted from neighbouring days"
        else:
            offset, source = None, ""
        clock_anchor = None if offset is None else unix[0] - offset
        if srun["at_midnight_start"]:
            gauge_anchor = day0 + MIDNIGHT_OFFSET_S
            how = "gauges: midnight file start + 0.5 s"
            if clock_anchor is not None and abs(clock_anchor - gauge_anchor) > CLOCK_OVERRIDE_S:
                shift = round(clock_anchor - gauge_anchor)
                gauge_anchor += shift
                how += f", shifted {shift:+d} s by the instrument clock"
        elif clock_anchor is not None:
            gauge_anchor, how = clock_anchor, "gauges: instrument clock"
        else:
            entry["result"] = "not reconstructed: no clock offset (no stamps that day, no offset history)"
            report.append(entry)
            continue
        if clock_anchor is None:
            entry["result"] = "not reconstructed: no clock offset for the auxiliary records"
            report.append(entry)
            continue
        anchor = clock_anchor                                  # SCPR records: their own clock
        how += f"; SCPR: instrument clock, offset {source}"
        entry["anchor"] = how
        gauges_ok = {}
        for ch in ("PARO1", "PARO2"):
            runs = gauge_runs.get(ch, [])
            if k >= len(runs):
                gauges_ok[ch] = "no matching gauge run"
                continue
            n = len(runs[k]["payloads"])
            diff = n - len(unix)
            entry[f"{ch}_records"] = n
            if abs(diff) > MAX_COUNT_DIFF or clock_gap > MAX_CLOCK_GAP_S:
                gauges_ok[ch] = f"not reconstructed: {ch} {diff:+d} lines vs SCPR, clock gap {clock_gap} s"
                continue
            idx = np.minimum(np.arange(n), len(missing) - 1)
            out[ch].append((gauge_anchor + np.arange(n) + missing[idx], runs[k]["payloads"]))
            gauges_ok[ch] = "reconstructed"
        if clock_gap <= MAX_CLOCK_GAP_S:
            out["SCPR"].append((anchor + rel, srun["payloads"], srun["markers"]))
            entry["SCPR"] = "reconstructed"
        else:
            entry["SCPR"] = f"not reconstructed: clock gap {clock_gap} s"
        entry.update(gauges_ok)
        report.append(entry)
    return out, report
