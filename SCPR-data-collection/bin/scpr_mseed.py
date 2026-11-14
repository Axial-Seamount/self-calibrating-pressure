#!/usr/bin/env python3
"""Convert one UTC day of raw piweb SCPR logger files to MiniSEED for OO.AXCC1.

Sources (http://piweb.ooirsn.uw.edu/scpr/data/): SCPRAA301_PARO1 and _PARO2
(the two Paroscientific gauges) and SCPRAA301_SCPR (the auxiliary log). The
channel list follows the AXCC1 rows of `COSZO Instrument Data Channels.xlsx`:
location 40 holds gauge 1, LVDT1 and the single auxiliary sensors (plus LCE and
LEF); location 41 holds gauge 2 and LVDT2.

Timing. Logger timestamps are data-center arrival times. For each source stream:
  1. exact duplicate records are dropped;
  2. buffered bursts after a network outage are put back in their grid slots,
     but only when the burst exactly fills the outage (see repair_bursts);
  3. the day's sample interval is fitted (gap_algorithms 'anomaly') and traces
     are split with PREST's 'timing' segmenting at 0.75 sample (see SEGMENT_TOL),
     so every written sample lies within 0.75 s of its (repaired) timestamp and
     MiniSEED readers never re-join the traces.
Values are written unaveraged and uninterpolated: pressures in Pa, periods in
s, voltages in V, temperatures in degC, tilts in degrees, as in the spreadsheet.

Output: an SDS tree, <out>/YEAR/OO/AXCC1/CHA.D/OO.AXCC1.LOC.CHA.D.YEAR.DOY (one
file per channel-day), plus a JSON manifest per day.

Requires numpy and ObsPy (conda env ooi_env).
    python scpr_mseed.py --date 2020-05-04            # raw files in sample_data/2020-05-04/raw
    python scpr_mseed.py --date 2020-05-04 --download # fetch missing raw files from piweb first
"""

from __future__ import annotations

import argparse
from datetime import date, datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
from obspy import Stream, Trace, UTCDateTime, read

import gap_algorithms as ga
import scpr_retime as retime
import scpr_paths as paths
from download_scpr_day import BASE_URL, DEFAULT_OUT as SAMPLE_ROOT, GAUGES, discover_files, download, source_pattern


NETWORK, STATION = "OO", "AXCC1"
DATA_QUALITY = "D"
REC_LEN = 512
SP_NOMINAL = 1.0
# Readers join traces contiguous within 1/2 sample; allow a little for the stored rate precision.
READBACK_TOL_SAMPLES = 0.51
# Trace split threshold for PREST's 'timing' segmenting (fraction of a sample). PREST uses
# 0.5, but then every split leaves a gap of almost exactly -1/2 sample, right at the
# MiniSEED readers' join tolerance: ObsPy/libmseed join some of them and re-time the
# following samples, by up to 1 s after two joins (2022-07-22, clock rate off by 55 ppm).
# At 0.75 every split gap is well outside the join tolerance, so readers get back exactly
# the times written; samples stay within 0.75 s of their (arrival-time) stamps.
SEGMENT_TOL = 0.75
# Bump when the conversion changes: scpr.py then re-converts days made by an older version.
CONVERTER_VERSION = 4   # 4: LE9 piston height (LVDT coefficients)
DEFAULT_MSEED = paths.MSEED
DEFAULT_MANIFEST = paths.MANIFEST
LOG_LINE = re.compile(r"^<log-ts:([^>]+)>\s?(.*)$")
EPOCH = np.datetime64("1970-01-01T00:00:00", "us")

# Readme - Cabled SCPR Data files.docx (forms confirmed against the document's equation images).
BARO_KPA_PER_V, BARO_OFFSET_KPA = 23.9836, 49.7607      # P = a*V + b
PRT_DEGC_PER_V, PRT_OFFSET_DEGC = 1000.0, 245.8         # T = a*V - b
TILT = {"X": (10.377, 2.479), "Y": (10.312, 2.551)}     # deg = a*(V - b)
# LVDT piston height: Z[mm] = a*(V - b), a = sensitivity (mm/V), b = null (V output at
# 0 mm). Values from Glenn Sasagawa (Nov 2026; his processing uses lvdt*_sens = 1.222,
# lvdt*_null = 0.25); to be re-checked after the 2027 recovery. None skips LE9.
LVDT: dict[int, tuple[float, float] | None] = {1: (1.222, 0.25), 2: (1.222, 0.25)}

# LEF calibration flag: 1 from the selector valve switching to the dead-weight
# tester until it switches back to seawater. "Starting Calibration Run" follows
# the switch by ~2 s and is the fallback start if the switch line is missing.
CAL_ON = "Switching SV to DWT"
CAL_ON_FALLBACK = "Starting Calibration Run"
CAL_OFF = "Exiting: SV to Sea"
# Manual valve commands typed by an operator (e.g. deployment, 2018-07-06 07:50-08:17)
# print these lines without a <log-ts:...> prefix; they take the time of the
# preceding <OOI-CMD time> block (the response follows within ~1 s).
MANUAL_ON, MANUAL_OFF = "Switching SV to DWT", "Switching SV to SEA"
OOI_CMD = re.compile(r"^<OOI-CMD (\d{4}-\d\d-\d\dT[\d:.]+)Z?\s")

AUX_COLUMNS = ("unix", "lvdt1", "lvdt2", "baro", "prt", "rpm", "xtilt", "ytilt")

# The unit each channel is written in, from the conversions below. This, not the
# spreadsheet, is the reference that check_metadata.py tests the StationXML against.
CHANNEL_UNITS = {
    "LDO": "Pa", "LDI": "Pa",                      # gauge / internal barometer pressure
    "LE1": "s", "LE2": "s",                        # gauge pressure / temperature periods
    "LE3": "V", "LE4": "V", "LE5": "V", "LE7": "V", "LE8": "V",
    "LE6": "rpm",                                  # mass-bell rotation
    "LE9": "m",                                    # piston height (lvdt_m: mm * 1e-3)
    "LK1": "C", "LKI": "C",                        # gauge / PRT temperature, degC
    "LAX": "degree", "LAY": "degree",              # tilt
    "LCE": "s",                                    # clock offset
    "LEF": "count",                                # calibration flag 0/1
}


# ── parsing ─────────────────────────────────────────────────────────────────
UNSTAMPED_RECORD = re.compile(r"^(\*\d{4},|1 \d{9,} )")   # a gauge or type-1 record line without <log-ts>


def read_log(paths: list[Path], tally: list[dict] | None = None):
    """(timestamp text, payload) for every <log-ts:...> line, in file order.

    Lines without an arrival timestamp cannot be timed and are not returned. Some piweb
    periods were logged that way (2024-04-11 19:22 .. 2024-04-17 partway, and from
    2026-10-17 20:56 UTC), sometimes for a whole file, sometimes for part of one. They
    are never dropped silently: if `tally` is given, one entry per file is appended with
    its stamped and untimestamped record counts, and the caller records it.
    """
    for path in paths:
        stamped = unstamped = 0
        with path.open("rb") as handle:
            for raw in handle:
                text = raw.decode("ascii", "replace").rstrip("\r\n")
                # Some logger lines were partly overwritten with NUL bytes and run into the
                # next record ("2024-12-\0\0...\0<log-ts:2024-12-12T04:06:08.339> *0002,...",
                # 2024-12-12): keep the last complete record on the line.
                if "\0" in text or text.count("<log-ts:") > 1:
                    text = text.replace("\0", "")
                    text = text[text.rfind("<log-ts:"):] if "<log-ts:" in text else text
                match = LOG_LINE.match(text)
                if match:
                    stamped += 1
                    yield match[1], match[2]
                elif UNSTAMPED_RECORD.match(text):
                    unstamped += 1
        if tally is not None:
            tally.append({"file": path.name, "stamped_lines": stamped, "untimestamped_records": unstamped})


def require_stamped(channel: str, n_records: int, tally: list[dict], why: str = "") -> None:
    """A channel-day with records only in untimestamped form (and not re-timed) cannot be converted."""
    lost = sum(f["untimestamped_records"] for f in tally)
    if not n_records and lost:
        raise ValueError(f"{channel}: all {lost} records lack <log-ts:...> arrival timestamps "
                         f"(piweb logged them untimestamped) and could not be re-timed"
                         + (f": {why}" if why else ""))


def operator_marks(paths: list[Path]) -> list[tuple[str, int]]:
    """(OOI-CMD time, 1 = to DWT / 0 = to sea) for valve messages without <log-ts>."""
    marks, command_time = [], None
    for path in paths:
        with path.open("rb") as handle:
            for raw in handle:
                line = raw.decode("ascii", "replace").rstrip("\r\n")
                cmd = OOI_CMD.match(line)
                if cmd:
                    command_time = cmd[1]
                elif command_time and not line.startswith("<log-ts:"):
                    if MANUAL_ON in line:
                        marks.append((command_time, 1))
                    elif MANUAL_OFF in line:
                        marks.append((command_time, 0))
    return marks


def gauge_row(payload: str, channel: str) -> list[float] | None:
    parts = payload[1:].split(",")
    try:
        row = [float(x) for x in parts[1:]] if len(parts) == 4 and parts[0] == GAUGES[channel]["id"] else None
    except ValueError:
        return None
    return row if row is not None and np.isfinite(row).all() else None


def scpr_row(payload: str) -> list[float] | None:
    fields = payload.split()
    try:
        row = [float(x) for x in fields[1:]] if len(fields) == 9 and fields[0] == "1" else None
    except ValueError:
        return None
    return row if row is not None and np.isfinite(row).all() else None


def add_reconstructed(rec: dict, pieces: list, parse, ncol: int, day0: float) -> int:
    """Merge reconstructed (times, payloads) into a parsed record set; returns records added."""
    t_new, v_new = [], []
    for times, payloads, *_ in pieces:
        for t, p in zip(times, payloads):
            row = parse(p)
            if row is not None:
                t_new.append(t)
                v_new.append(row)
    if not t_new:
        return 0
    t = np.concatenate([rec["t"], np.array(t_new)])
    v = np.vstack([rec["values"].reshape(-1, ncol), np.array(v_new, dtype=float)])
    before = len(rec["t"])
    rec["t"], rec["values"], n_dup = day_records(t, v, day0)
    rec["duplicates_dropped"] += n_dup
    return len(rec["t"]) - before


def to_epoch(stamps: list[str]) -> np.ndarray:
    """Logger timestamps (UTC, no zone suffix) to float seconds since 1970."""
    if any("+" in s or s.count("-") > 2 for s in stamps):
        raise ValueError("unexpected time-zone offset in logger timestamp")
    values = np.array([s.rstrip("Z") for s in stamps], dtype="datetime64[us]")
    return (values - EPOCH).astype(np.int64) / 1e6


def day_records(t: np.ndarray, values: np.ndarray, day0: float) -> tuple[np.ndarray, np.ndarray, int]:
    """Records of the UTC day, stably time-sorted, exact duplicate records removed."""
    keep = (t >= day0) & (t < day0 + 86400)
    t, values = t[keep], values[keep]
    order = np.argsort(t, kind="stable")
    t, values = t[order], values[order]
    same = np.zeros(len(t), dtype=bool)
    same[1:] = (np.diff(t) == 0) & np.all(values[1:] == values[:-1], axis=1)
    return t[~same], values[~same], int(same.sum())


def parse_gauge(paths: list[Path], channel: str, day0: float, allow_unstamped: bool = False) -> dict:
    """Gauge records: pressure (kPa), pressure period (us), temperature period (us)."""
    ident = GAUGES[channel]["id"]
    stamps, rows, malformed, tally = [], [], [], []
    for stamp, payload in read_log(paths, tally):
        if not payload.startswith("*"):
            continue
        # The log has occasional truncated records (e.g. "*0001,5934,5.8547578"):
        # skip anything that is not "*ID,pressure,period,period" with finite numbers.
        parts = payload[1:].split(",")
        try:
            row = [float(x) for x in parts[1:]] if len(parts) == 4 and parts[0] == ident else None
        except ValueError:
            row = None
        if row is None or not np.isfinite(row).all():
            malformed.append(f"{stamp} {payload[:60]}")
            continue
        stamps.append(stamp)
        rows.append(row)
    if not allow_unstamped:
        require_stamped(channel, len(rows), tally)
    values = np.array(rows, dtype=float).reshape(-1, 3)
    t, values, n_dup = day_records(to_epoch(stamps), values, day0)
    return {"t": t, "values": values, "duplicates_dropped": n_dup, "malformed_records": malformed,
            "log_files": tally}


def parse_scpr(paths: list[Path], day0: float, allow_unstamped: bool = False) -> dict:
    """Type-1 auxiliary records plus the calibration on/off log lines.

    Operator menus are echoed into the same log (e.g. "1    Auto Calibration ..."),
    so a line counts as a record only if it has 9 fields that all parse as numbers.
    """
    stamps, rows, event_stamps, events, tally = [], [], [], [], []
    for stamp, payload in read_log(paths, tally):
        fields = payload.split()
        if len(fields) == 9 and fields[0] == "1":
            try:
                row = [float(x) for x in fields[1:]]
            except ValueError:
                continue
            stamps.append(stamp)
            rows.append(row)
        elif CAL_ON in payload or CAL_ON_FALLBACK in payload or CAL_OFF in payload:
            event_stamps.append(stamp)
            events.append(0 if CAL_OFF in payload else 1)
    if not allow_unstamped:
        require_stamped("SCPR", len(rows), tally)
    values = np.array(rows, dtype=float).reshape(-1, len(AUX_COLUMNS))
    if not np.isfinite(values).all():
        raise ValueError("SCPR: non-finite auxiliary value")
    t, values, n_dup = day_records(to_epoch(stamps), values, day0)
    for stamp, mark in operator_marks(paths):
        event_stamps.append(stamp)
        events.append(mark)
    et = to_epoch(event_stamps) if event_stamps else np.array([])
    inday = (et >= day0) & (et < day0 + 86400)
    marks = sorted(zip(et[inday].tolist(), np.array(events)[inday].tolist())) if len(et) else []
    return {"t": t, "values": values, "duplicates_dropped": n_dup, "cal_marks": marks, "log_files": tally}


# ── timing ──────────────────────────────────────────────────────────────────
def repair_bursts(t: np.ndarray, sp: float = SP_NOMINAL) -> tuple[np.ndarray, list[dict]]:
    """Put buffered samples back into the outage they belong to.

    After a network outage the logger delivers the held samples in one burst
    (spacing << sp), stamped on arrival. If the burst exactly fills the outage
    (burst count == round(span/sp) - 1, where span runs from the last sample
    before the outage to the first regular sample after the burst), the burst
    samples are spaced evenly across that span. Anything else is left alone;
    the timing segmenter then keeps those samples near their stamps.
    """
    t = t.copy()
    repaired = []
    for i in np.flatnonzero(np.diff(t) > 1.5 * sp):
        j = k = i + 1
        while k + 1 < len(t) and t[k + 1] - t[k] < 0.5 * sp:
            k += 1
        count = k - j + 1
        if count < 2 or k + 1 >= len(t):
            continue
        span = t[k + 1] - t[i]
        if count == round(span / sp) - 1:
            moved = t[i] + np.arange(1, count + 1) * span / (count + 1)
            repaired.append({"outage_start_utc": iso(t[i]), "samples": int(count),
                             "max_shift_s": round(float(np.max(np.abs(moved - t[j:k + 1]))), 3)})
            t[j:k + 1] = moved
    return t, repaired


def segment(t: np.ndarray) -> dict:
    """Fitted sample rate and trace splits (PREST: gap_algo anomaly, mseed_segmenting timing)."""
    rel = t - t[0]
    result = ga.detect_gaps("anomaly", rel, SP_NOMINAL)
    splits, starts = ga.mseed_segments("timing", rel, result, tol=SEGMENT_TOL)
    starts = [t[0] + s for s in starts]
    grid = np.concatenate([s + np.arange(len(piece)) * result.sp
                           for s, piece in zip(starts, np.split(t, splits))])
    diag = result.diagnostics
    return {"sr": result.sr, "sp": result.sp, "splits": splits, "starts": starts,
            "summary": {"samples": int(len(t)), "dt_true_s": result.sp,
                        "true_missing": int(diag.get("true_missing", 0)),
                        "jitter_sigma_ms": round(float(diag.get("sigma_ms", np.nan)), 3),
                        "traces": len(splits) + 1,
                        "max_abs_grid_minus_timestamp_s": round(float(np.max(np.abs(grid - t))), 4)}}


# ── conversions ─────────────────────────────────────────────────────────────
def paro_temperature(period_us: np.ndarray, channel: str) -> np.ndarray:
    """Paroscientific temperature (degC) from the temperature period (us)."""
    c = GAUGES[channel]
    u = period_us - c["U0"]
    return c["Y1"] * u + c["Y2"] * u**2 + c["Y3"] * u**3


def lvdt_m(volts: np.ndarray, sensor: int) -> np.ndarray | None:
    coeff = LVDT[sensor]
    return None if coeff is None else coeff[0] * (volts - coeff[1]) * 1e-3


def gauge_channels(channel: str, values: np.ndarray) -> dict[tuple[str, str], np.ndarray]:
    loc = "40" if channel == "PARO1" else "41"
    p_kpa, p_period, t_period = values.T
    return {(loc, "LDO"): p_kpa * 1000.0, (loc, "LE1"): p_period * 1e-6,
            (loc, "LE2"): t_period * 1e-6, (loc, "LK1"): paro_temperature(t_period, channel)}


def aux_channels(values: np.ndarray) -> tuple[dict[tuple[str, str], np.ndarray], list[str]]:
    col = dict(zip(AUX_COLUMNS, values.T))
    out = {("40", "LE3"): col["lvdt1"], ("41", "LE3"): col["lvdt2"],
           ("40", "LE4"): col["baro"], ("40", "LE5"): col["prt"], ("40", "LE6"): col["rpm"],
           ("40", "LE7"): col["xtilt"], ("40", "LE8"): col["ytilt"],
           ("40", "LDI"): (BARO_KPA_PER_V * col["baro"] + BARO_OFFSET_KPA) * 1000.0,
           ("40", "LKI"): PRT_DEGC_PER_V * col["prt"] - PRT_OFFSET_DEGC,
           ("40", "LAX"): TILT["X"][0] * (col["xtilt"] - TILT["X"][1]),
           ("40", "LAY"): TILT["Y"][0] * (col["ytilt"] - TILT["Y"][1])}
    skipped = []
    for sensor, loc in ((1, "40"), (2, "41")):
        z = lvdt_m(col[f"lvdt{sensor}"], sensor)
        if z is None:
            skipped.append(f"{loc}.LE9 (LVDT{sensor} coefficients unknown)")
        else:
            out[(loc, "LE9")] = z
    return out, skipped


def clock_offset(t: np.ndarray, unix: np.ndarray, day0: float) -> tuple[np.ndarray, np.ndarray]:
    """LCE: per-minute median of SCPR clock minus timestamp, stamped at minute centres."""
    minute = np.floor((t - day0) / 60).astype(int)
    offset = unix - t
    bins = np.unique(minute)
    med = np.array([np.median(offset[minute == b]) for b in bins])
    return day0 + bins * 60.0 + 30.0, med


def calibration_flag(marks: list[tuple[float, int]], day0: float) -> tuple[np.ndarray, list[dict]]:
    """LEF at whole seconds: 1 between CAL_ON and CAL_OFF log lines, else 0.

    A day whose first marker is CAL_OFF starts flagged (calibration across midnight);
    an unclosed CAL_ON stays flagged to the end of the day.
    """
    flag = np.zeros(86400, dtype=np.int32)
    state = 1 if marks and marks[0][1] == 0 else 0
    since = day0 if state else None
    intervals = []
    for when, on in marks + [(day0 + 86400, None)]:
        if state and (on == 0 or on is None):
            a, b = int(np.ceil(since - day0)), int(np.ceil(when - day0))
            flag[a:b] = 1
            intervals.append({"start_utc": iso(since), "end_utc": iso(when) if on is not None else None})
            state = 0
        elif not state and on == 1:
            state, since = 1, when
    return flag, intervals


# ── output ──────────────────────────────────────────────────────────────────
def iso(seconds: float) -> str:
    return UTCDateTime(seconds).isoformat() + "Z"


def make_traces(loc: str, cha: str, data: np.ndarray, seg: dict) -> list[Trace]:
    traces = []
    for start, piece in zip(seg["starts"], np.split(data, seg["splits"])):
        if len(piece):
            traces.append(Trace(data=np.ascontiguousarray(piece), header={
                "network": NETWORK, "station": STATION, "location": loc, "channel": cha,
                "starttime": UTCDateTime(start), "sampling_rate": seg["sr"],
                "mseed": {"dataquality": DATA_QUALITY}}))
    return traces


def sds_path(root: Path, loc: str, cha: str, day: date) -> Path:
    doy = day.timetuple().tm_yday
    return root / f"{day.year}" / NETWORK / STATION / f"{cha}.D" / \
        f"{NETWORK}.{STATION}.{loc}.{cha}.D.{day.year}.{doy:03d}"


def _samples(traces) -> tuple[np.ndarray, np.ndarray]:
    """Every sample as (time, value), sorted by value and then time."""
    t = np.concatenate([tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate
                        for tr in traces])
    v = np.concatenate([tr.data for tr in traces])
    order = np.lexsort((t, v))
    return t[order], v[order]


def write_channel(path: Path, traces: list[Trace], expected: np.ndarray) -> dict:
    """Write one channel-day atomically, then read it back and check every sample.

    Readers (ObsPy, libmseed) join or split traces that are contiguous within half a
    sample, so on days with jittery arrivals the trace layout read back differs from
    the one written and joined samples move by up to half a sample. The check is
    therefore layout-free: the same multiset of values must come back, and pairing
    equal values in time order (which minimises the largest pairwise shift) must
    keep every sample within READBACK_TOL_SAMPLES of where it was written.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    Stream(traces).write(str(tmp), format="MSEED", reclen=REC_LEN)
    back = read(str(tmp), format="MSEED")
    t_in, v_in = _samples(traces)
    t_out, v_out = _samples(back)
    limit = READBACK_TOL_SAMPLES / min(tr.stats.sampling_rate for tr in traces)
    ok = len(v_in) == len(expected) == len(v_out) and np.array_equal(v_in, v_out)
    shift = float(np.abs(t_in - t_out).max(initial=0.0)) if ok else math.inf
    if not ok or shift > limit:
        tmp.unlink()
        raise ValueError(f"{path.name}: read-back check failed (values {'match' if ok else 'differ'}, "
                         f"max time shift {shift:.3f} s)")
    tmp.replace(path)
    return {"file": str(path), "traces": len(back), "samples": int(len(v_out)),
            "encoding": back[0].stats.mseed.encoding, "max_readback_time_shift_s": round(shift, 6)}


def local_sources_complete(manifest: Path, raw_dir: Path) -> bool:
    """True if an earlier conversion's manifest lists source files that are all still on disk.

    A day converted once had its complete piweb file list downloaded, so re-converting it
    (e.g. after a converter change) needs no network access.
    """
    try:
        sources = json.loads(manifest.read_text())["sources"]
    except (OSError, ValueError, KeyError):
        return False
    files = [Path(f["path"]).name for paths in sources.values() for f in paths]
    return bool(files) and all((raw_dir / name).is_file() for name in files)


def raw_files(raw_dir: Path, channel: str, day: date, fetch: bool) -> list[Path]:
    if fetch:
        raw_dir.mkdir(parents=True, exist_ok=True)
        try:
            urls = discover_files(BASE_URL, channel, day)
        except ValueError:                      # piweb has no file for this channel-day
            urls = []
        for url in urls:
            try:
                download(url, raw_dir / url.rsplit("/", 1)[1])
            except ValueError as exc:           # piweb has some empty (0-byte) files
                print(f"skipped: {exc}")
    return sorted(p for p in raw_dir.glob("*.dat") if re.fullmatch(source_pattern(channel, day), p.name))


def convert_day(day: date, raw_dir: Path, out: Path = DEFAULT_MSEED, manifest_dir: Path = DEFAULT_MANIFEST,
                download: bool = False, verbose: bool = True) -> dict:
    """Convert one UTC day; returns the manifest (also written to manifest_dir)."""
    day0 = datetime(day.year, day.month, day.day, tzinfo=timezone.utc).timestamp()
    raw_dir = raw_dir.expanduser().resolve()
    if download and local_sources_complete(manifest_dir / f"scpr_mseed_{day}.json", raw_dir):
        download = False                    # re-conversion from files already on disk: no piweb
    sources = {ch: raw_files(raw_dir, ch, day, download) for ch in ("PARO1", "PARO2", "SCPR")}
    if not any(sources.values()):
        raise ValueError(f"No raw SCPR files for {day} in {raw_dir}")

    manifest = {"date_utc": str(day), "network": NETWORK, "station": STATION,
                "converter_version": CONVERTER_VERSION, "segment_tol_samples": SEGMENT_TOL,
                "timing": "duplicates dropped; exact-fill bursts repaired; gap_algorithms anomaly + timing segments",
                "sources": {ch: [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                                 for p in paths] for ch, paths in sources.items()},
                "streams": {}, "channels": {}, "skipped": []}
    outputs: dict[tuple[str, str], tuple[list[Trace], np.ndarray]] = {}

    # Records logged without <log-ts> stamps are re-timed from the instrument clock
    # (scpr_retime); everything else keeps its arrival-time timing.
    scpr_rec = parse_scpr(sources["SCPR"], day0, allow_unstamped=True) if sources["SCPR"] else None
    gauge_runs = {ch: retime.unstamped_runs(sources[ch], retime.GAUGE_RECORD) for ch in ("PARO1", "PARO2")}
    scpr_runs = (retime.unstamped_runs(sources["SCPR"], retime.SCPR_RECORD, (CAL_ON, CAL_ON_FALLBACK, CAL_OFF))
                 if sources["SCPR"] else [])
    recon, recon_report = {"PARO1": [], "PARO2": [], "SCPR": []}, []
    if scpr_runs or any(gauge_runs.values()):
        stamped = (scpr_rec["t"], scpr_rec["values"][:, 0]) if scpr_rec else (np.array([]), np.array([]))
        offset, offset_how = retime.predict_offset(retime.offset_history(manifest_dir), day)
        recon, recon_report = retime.reconstruct(gauge_runs, scpr_runs, day0, stamped, offset)
        manifest["timing_reconstructed"] = {"method": "scpr_retime (instrument clock; see README)",
                                            "clock_offset": {"value_s": offset, "how": offset_how},
                                            "validation": retime.VALIDATION, "runs": recon_report,
                                            "records": {}}

    for channel in ("PARO1", "PARO2"):
        if not sources[channel]:
            manifest["skipped"].append(f"{channel}: no raw file")
            continue
        rec = parse_gauge(sources[channel], channel, day0, allow_unstamped=True)
        added = add_reconstructed(rec, recon[channel], lambda p, ch=channel: gauge_row(p, ch), 3, day0)
        if added:
            manifest["timing_reconstructed"]["records"][channel] = added
        if not len(rec["t"]):
            why = "; ".join(str(r.get(channel) or r.get("result")) for r in recon_report) or "no reconstruction attempted"
            require_stamped(channel, 0, rec["log_files"], why)
        t, bursts = repair_bursts(rec["t"])
        seg = segment(t)
        manifest["streams"][channel] = {**seg["summary"], "duplicates_dropped": rec["duplicates_dropped"],
                                        "malformed_records_skipped": rec["malformed_records"],
                                        "log_files": rec["log_files"], "reconstructed_records": added,
                                        "bursts_repaired": bursts}
        for key, data in gauge_channels(channel, rec["values"]).items():
            outputs[key] = (make_traces(*key, data, seg), data)

    if scpr_rec is not None:
        rec = scpr_rec
        stamped_t, stamped_unix = rec["t"].copy(), rec["values"][:, 0].copy()   # LCE: measured, not reconstructed
        added = add_reconstructed(rec, recon["SCPR"], scpr_row, len(AUX_COLUMNS), day0)
        if added:
            manifest["timing_reconstructed"]["records"]["SCPR"] = added
            # calibration lines inside re-timed runs: time of the record before them + 0.5 s
            for times, _, markers in recon["SCPR"]:
                for i, text in markers:
                    when = float(times[max(i, 0)]) + (0.5 if i >= 0 else -0.5)
                    if day0 <= when < day0 + 86400:
                        rec["cal_marks"].append((when, 0 if CAL_OFF in text else 1))
            rec["cal_marks"].sort()
        if not len(rec["t"]):
            why = "; ".join(str(r.get("SCPR") or r.get("result")) for r in recon_report) or "no reconstruction attempted"
            require_stamped("SCPR", 0, rec["log_files"], why)
        t, bursts = repair_bursts(rec["t"])
        seg = segment(t)
        channels, skipped = aux_channels(rec["values"])
        manifest["skipped"] += skipped
        for key, data in channels.items():
            outputs[key] = (make_traces(*key, data, seg), data)
        # LCE: one value per minute; contiguous minutes share a trace. Only from records with
        # real arrival stamps: on re-timed records it would just echo the predicted offset.
        lce_t, lce = clock_offset(stamped_t, stamped_unix, day0) if len(stamped_t) else (np.array([]), np.array([]))
        if len(lce):
            breaks = (np.flatnonzero(np.diff(lce_t) > 60.5) + 1).tolist()
            lce_seg = {"sr": 1 / 60, "splits": breaks, "starts": [lce_t[i] for i in [0] + breaks]}
            outputs[("40", "LCE")] = (make_traces("40", "LCE", lce, lce_seg), lce)
        else:
            manifest["skipped"].append("40.LCE: no timestamped records (clock offset only measurable from arrival stamps)")
        flag, intervals = calibration_flag(rec["cal_marks"], day0)
        outputs[("40", "LEF")] = (make_traces("40", "LEF", flag, {"sr": 1.0, "splits": [], "starts": [day0]}), flag)
        manifest["streams"]["SCPR"] = {**seg["summary"], "duplicates_dropped": rec["duplicates_dropped"],
                                       "log_files": rec["log_files"], "reconstructed_records": added,
                                       "bursts_repaired": bursts, "calibration_flag_intervals": intervals,
                                       "clock_offset_s": ({"min": float(lce.min()), "median": float(np.median(lce)),
                                                           "max": float(lce.max()),
                                                           # when the offset was measured (stamped part of the day)
                                                           "t_mean_day_fraction": float((np.mean(lce_t) - day0) / 86400)}
                                                          if len(lce) else None)}
    else:
        manifest["skipped"].append("SCPR: no raw file (auxiliary channels, LCE and LEF)")

    lost = {name: sum(f["untimestamped_records"] for f in st.get("log_files", [])) - st.get("reconstructed_records", 0)
            for name, st in manifest["streams"].items()}
    manifest["untimestamped_records_not_converted"] = {k: v for k, v in lost.items() if v}
    for (loc, cha), (traces, data) in sorted(outputs.items()):
        manifest["channels"][f"{loc}.{cha}"] = write_channel(sds_path(out.resolve(), loc, cha, day), traces, data)

    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / f"scpr_mseed_{day}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if verbose:
        for name, s in manifest["streams"].items():
            print(f"{name}: {s['samples']} samples, {s['traces']} traces, dt={s['dt_true_s']:.7f} s, "
                  f"{len(s['bursts_repaired'])} bursts repaired, "
                  f"max |grid - stamp| = {s['max_abs_grid_minus_timestamp_s']} s")
        print(f"{len(manifest['channels'])} channels written to {out}; skipped: {manifest['skipped'] or 'none'}")
        print(f"Manifest: {manifest_path}")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--date", type=date.fromisoformat, required=True, help="UTC day YYYY-MM-DD")
    parser.add_argument("--raw-dir", type=Path, help="Folder with the piweb .dat files (default: sample_data/DATE/raw)")
    parser.add_argument("--out", type=Path, default=DEFAULT_MSEED, help="SDS root for MiniSEED")
    parser.add_argument("--manifest-dir", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--download", action="store_true", help="Download missing raw files from piweb first")
    args = parser.parse_args()
    convert_day(args.date, args.raw_dir or SAMPLE_ROOT / str(args.date) / "raw", args.out, args.manifest_dir,
                args.download)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
