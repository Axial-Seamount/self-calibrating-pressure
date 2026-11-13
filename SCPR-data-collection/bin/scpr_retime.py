"""Reconstruct sample times for SCPR records that piweb logged without timestamps.

Some periods were logged with no <log-ts:...> arrival stamps (2024-04-11 19:22 to
2024-04-17 16:19 and 2026-10-17 20:56 to 2026-11-13 20:15 UTC, plus short stretches next
to network outages). The instrument kept recording: the gauges write one line per second
and every auxiliary (SCPR) record carries the instrument clock (UNIX_SEC).

Method (validated on 186 timestamped days by stripping the stamps and comparing with the
converter's normal times; see README "Periods without timestamps"):
  * instrument-clock offset (clock - UTC arrival): a straight line between clock resets,
    predicted from the daily offsets in the surrounding days' manifests (interpolated when
    there are stamps on both sides of the gap), or measured from the day's own stamps;
  * SCPR (auxiliary) records: time = instrument clock - offset, across outages too;
  * gauge records (no clock): one per second, counted forward from a start anchor
    (midnight file start + 0.5 s, whole seconds corrected by the clock only when it
    disagrees by > CLOCK_OVERRIDE_S; restart files: the clock) and checked against an end
    anchor where there is one (first real stamp after the run - 1 s, or next midnight - 0.5 s
    when the file and the clock run to midnight); short outages (<= MAX_CLOCK_GAP_S) are
    placed from the clock;
  * one longer outage inside a run: lines before it counted forward from the start anchor,
    lines after it counted back from the end anchor (each stream drops at the outage but
    reconnects at its own time); BREAK_MARGIN lines each side of the break are dropped;
  * acceptance: forward count within END_TOL_S of the end anchor, or (no end anchor) line
    counts within MAX_COUNT_DIFF of the SCPR run; and PARO1-PARO2 alignment: the two gauges
    record the same seafloor-pressure noise, so their cross-correlation must peak at lag 0
    (+-XCORR_LAG_TOL) in every informative window. Anything else is reported, not converted.
Reconstructed days are labelled in the manifest ("timing_reconstructed").
"""

from __future__ import annotations

from datetime import date, datetime, timezone
import json
from pathlib import Path
import re

import numpy as np

MIDNIGHT_OFFSET_S = 0.5          # first record of a midnight file arrives within its first second
MAX_COUNT_DIFF = 1               # gauge lines vs SCPR records in the same run
MAX_CLOCK_GAP_S = 5              # instrument-clock jumps longer than this (s) are outages
CLOCK_OVERRIDE_S = 1.5           # clock must disagree with the midnight anchor by more than this
RESET_JUMP_S = 3.0               # a day-to-day offset change larger than this is a clock reset
END_TOL_S = 1.5                  # counted gauge time vs its end anchor (next real stamp / next midnight)
MAX_END_SHIFT_S = 10             # clock-indicated early/late file end beyond this: no midnight end anchor
BREAK_MARGIN = 10                # gauge lines dropped on each side of an outage inside a run
XCORR_WINDOW_S = 1800            # PARO1-PARO2 alignment check: window length,
XCORR_MAX_LAG = 300              # lags searched,
XCORR_MIN_CORR = 0.3             # correlation needed for a window to count,
XCORR_LAG_TOL = 1                # and largest lag (s) accepted
VALIDATION = ("186 timestamped days, stamps removed (all day / until a random time): gauges 97.2% / 97.8% of "
              "gauge-day runs with every sample within 0.75 s of the normal timing (98.3% / 99.4% with >= 99% "
              "of samples); runs split around an outage 12/12 within 0.75 s; SCPR channels 94.0% / 96.8% of "
              "days fully within 0.75 s; clock-offset prediction median error 0.01 s")

STAMPED = re.compile(r"^<log-ts:")
LOG_TS = re.compile(r"<log-ts:([0-9T:.\-]+)>")
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


def _stamp_time(text: str) -> float | None:
    found = LOG_TS.findall(text)
    if not found:
        return None
    try:
        return datetime.fromisoformat(found[-1]).replace(tzinfo=timezone.utc).timestamp()
    except ValueError:
        return None


def unstamped_runs(paths: list[Path], record: re.Pattern, markers: tuple[str, ...] = ()) -> list[dict]:
    """Maximal runs of untimestamped record lines, per file, in order.

    Each run: file name, whether it starts at the first record of a midnight file,
    its record payloads, (index, text) of marker lines found inside it, the arrival
    stamp of the first timestamped record after it in the same file ("next_stamp",
    if the stamps came back), and whether it runs to the end of the day's last file
    ("ends_day_file").
    """
    runs = []
    for n_path, path in enumerate(paths):
        midnight = (m := FILE_TIME.search(path.name)) is not None and m[1] == "0000"
        current, pending, seen_record = None, None, False
        with path.open("rb") as handle:
            for raw in handle:
                text = raw.decode("ascii", "replace").rstrip("\r\n").replace("\0", "")
                if STAMPED.match(text):
                    seen_record = True
                    if current is not None:
                        pending, current = current, None
                    body = text[text.rfind(">") + 1:].lstrip()
                    if pending is not None and record.match(body):
                        pending["next_stamp"], pending["next_payload"] = _stamp_time(text), body
                        pending = None
                    continue
                if record.match(text):
                    if current is None:
                        current = {"file": path.name, "at_midnight_start": midnight and not seen_record,
                                   "payloads": [], "markers": [], "next_stamp": None, "next_payload": None,
                                   "ends_day_file": False}
                        runs.append(current)
                        pending = None
                    current["payloads"].append(text)
                    seen_record = True
                elif current is not None and any(k in text for k in markers):
                    current["markers"].append((len(current["payloads"]) - 1, text))
        if current is not None and n_path == len(paths) - 1:
            current["ends_day_file"] = True
    return runs


def _gauge_pressure(payloads: list[str]) -> np.ndarray:
    out = np.full(len(payloads), np.nan)
    for i, p in enumerate(payloads):
        try:
            out[i] = float(p.split(",")[1])
        except (IndexError, ValueError):
            pass
    return out


def _highpass(x: np.ndarray, width: int = 31) -> np.ndarray:
    ok = np.isfinite(x)
    k = np.ones(width)
    mean = np.convolve(np.where(ok, x, 0.0), k, "same") / np.maximum(np.convolve(ok.astype(float), k, "same"), 1)
    return np.where(ok, x - mean, np.nan)


def gauge_alignment(t1: np.ndarray, p1: np.ndarray, t2: np.ndarray, p2: np.ndarray) -> dict:
    """Relative timing of the two gauges from their common pressure fluctuations.

    Both gauges sit side by side and record the same seafloor-pressure noise (microseisms,
    infragravity waves), so on correctly timed data the cross-correlation of the two
    high-passed series peaks at lag 0 in every window. Returns the windows checked and the
    lags found (seconds; positive = PARO2 late).
    """
    start = np.floor(min(t1.min(), t2.min()))
    n = int(np.ceil(max(t1.max(), t2.max()) - start)) + 1
    if n < XCORR_WINDOW_S:                                     # too short to check
        return {"windows": 0, "lags": [], "ok": True}
    grid = [np.full(n, np.nan), np.full(n, np.nan)]
    for g, t, p in ((grid[0], t1, p1), (grid[1], t2, p2)):
        g[np.clip(np.floor(t - start).astype(int), 0, n - 1)] = p
    a, b = _highpass(grid[0]), _highpass(grid[1])
    lags, informative = [], 0
    for w0 in range(0, n - XCORR_WINDOW_S + 1, XCORR_WINDOW_S):
        x = a[w0:w0 + XCORR_WINDOW_S]
        if np.isfinite(x).sum() < XCORR_WINDOW_S // 2:
            continue
        best = (-1.0, 0)
        for lag in range(-XCORR_MAX_LAG, XCORR_MAX_LAG + 1):
            lo, hi = w0 + lag, w0 + lag + XCORR_WINDOW_S
            if lo < 0 or hi > n:
                continue
            y = b[lo:hi]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < XCORR_WINDOW_S // 2:
                continue
            c = np.corrcoef(x[ok], y[ok])[0, 1]
            if c > best[0]:
                best = (c, lag)
        if best[0] >= XCORR_MIN_CORR:
            informative += 1
            lags.append(best[1])
    return {"windows": informative, "lags": sorted(set(lags)),
            "ok": all(abs(x) <= XCORR_LAG_TOL for x in lags)}


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
        jumps = np.diff(unix) - 1                              # seconds missing between records
        outages = np.flatnonzero(jumps > MAX_CLOCK_GAP_S)      # real outages inside the run
        entry = {"run": k, "file": srun["file"], "scpr_records": len(unix), "clock_gap_s": clock_gap,
                 "outages": [{"after_record": int(i), "missing_s": int(jumps[i])} for i in outages]}
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
        t_scpr = anchor + rel
        how += f"; SCPR: instrument clock, offset {source}"
        entry["anchor"] = how
        # End of the day's last file: next midnight - 0.5 s, corrected by whole seconds when the
        # clock says the logger switched files early/late (as for the start anchor); unusable if
        # the clock shows the stream stopped well before midnight.
        midnight_end, midnight_how = None, ""
        if srun.get("ends_day_file"):
            midnight_end, midnight_how = day0 + 86400 - MIDNIGHT_OFFSET_S, "next midnight - 0.5 s"
            d_end = t_scpr[-1] - midnight_end
            if abs(d_end) > MAX_END_SHIFT_S:
                midnight_end = None
            elif abs(d_end) > CLOCK_OVERRIDE_S:
                midnight_end += round(d_end)
                midnight_how += f", shifted {round(d_end):+d} s by the instrument clock"
        # seconds missing after each SCPR record, counted back from the run's last record
        missing_back = (rel[-1] - rel) - (len(rel) - 1 - np.arange(len(rel)))
        # Stamps resume right after the run (no outage in between, by the SCPR clock): only then
        # is a gauge's first real stamp an end anchor for its unstamped lines.
        try:
            stamps_resume = (srun.get("next_payload") is not None
                             and int(srun["next_payload"].split()[1]) - unix[-1] <= 1 + MAX_CLOCK_GAP_S)
        except (IndexError, ValueError):
            stamps_resume = False
        gauges_ok, staged = {}, {}
        for ch in ("PARO1", "PARO2"):
            runs = gauge_runs.get(ch, [])
            if k >= len(runs):
                gauges_ok[ch] = "no matching gauge run"
                continue
            grun = runs[k]
            n = len(grun["payloads"])
            diff = n - len(unix)
            entry[f"{ch}_records"] = n
            if grun.get("next_stamp") is not None and stamps_resume:
                end_anchor, end_how = grun["next_stamp"] - 1.0, "first real stamp after the run - 1 s"
            elif grun.get("ends_day_file") and midnight_end is not None:
                end_anchor, end_how = midnight_end, midnight_how
            else:
                end_anchor, end_how = None, ""
            idx = np.arange(n)
            forward = gauge_anchor + idx + missing[np.minimum(idx, len(missing) - 1)]
            if (not srun["at_midnight_start"] and not len(outages) and grun.get("next_stamp") is not None
                    and stamps_resume):
                # Restart file: the lines flushed from the terminal server's buffer on reconnection
                # run straight into the stamped lines, but each port reconnects at its own moment,
                # so count back from this gauge's own first stamp, not from the SCPR clock.
                staged[ch] = [(end_anchor - (n - 1 - idx), grun["payloads"])]
                gauges_ok[ch] = f"reconstructed (counted back from {end_how})"
                continue
            if not srun["at_midnight_start"] and len(outages):
                gauges_ok[ch] = "not reconstructed: outage inside a run that does not start at midnight"
                continue
            if len(outages) == 0:
                end_err = None if end_anchor is None else float(forward[-1] - end_anchor)
                if end_err is not None and abs(end_err) > END_TOL_S:
                    gauges_ok[ch] = (f"not reconstructed: {ch} counted to {end_err:+.1f} s from the end anchor "
                                     f"({end_how}); {diff:+d} lines vs SCPR")
                    continue
                if abs(diff) > MAX_COUNT_DIFF and end_err is None:
                    gauges_ok[ch] = f"not reconstructed: {ch} {diff:+d} lines vs SCPR, no end anchor"
                    continue
                staged[ch] = [(forward, grun["payloads"])]
                gauges_ok[ch] = "reconstructed" + ("" if end_err is None else
                                                   f" (end within {end_err:+.2f} s of {end_how})")
                continue
            if len(outages) > 1:
                gauges_ok[ch] = f"not reconstructed: {len(outages)} outages inside the run"
                continue
            if end_anchor is None:
                gauges_ok[ch] = "not reconstructed: outage inside the run and no end anchor"
                continue
            # One outage: lines before it counted forward from the start anchor, lines after it
            # counted back from the end anchor; the clock says where the outage starts.
            g = int(outages[0])
            n_before = int(round(t_scpr[g] - gauge_anchor - missing[g])) + 1
            j = n - 1 - idx                                    # lines before the run's last line
            backward = end_anchor - j - missing_back[np.maximum(len(rel) - 1 - j, g + 1)]
            keep_f = idx < n_before - BREAK_MARGIN
            keep_b = idx >= n_before + BREAK_MARGIN
            if not keep_f.any() or not keep_b.any():
                gauges_ok[ch] = "not reconstructed: outage too close to the run's ends"
                continue
            gauge_outage = float(backward[n_before] - forward[n_before - 1] - 1) if 0 < n_before < n else -1.0
            if gauge_outage < -END_TOL_S:
                gauges_ok[ch] = (f"not reconstructed: {ch} has {-gauge_outage:.0f} more lines than the time "
                                 f"between the anchors allows")
                continue
            p = grun["payloads"]
            staged[ch] = [(forward[keep_f], [p[i] for i in np.flatnonzero(keep_f)]),
                          (backward[keep_b], [p[i] for i in np.flatnonzero(keep_b)])]
            gauges_ok[ch] = (f"reconstructed in two parts around the {int(jumps[g])}-s outage (gauge outage "
                             f"{gauge_outage:.0f} s; {2 * BREAK_MARGIN} lines at the break dropped; "
                             f"end anchor: {end_how})")
        if len(staged) == 2:
            t1 = np.concatenate([t for t, _ in staged["PARO1"]]); t2 = np.concatenate([t for t, _ in staged["PARO2"]])
            p1 = _gauge_pressure([x for _, ps in staged["PARO1"] for x in ps])
            p2 = _gauge_pressure([x for _, ps in staged["PARO2"] for x in ps])
            check = gauge_alignment(t1, p1, t2, p2)
            entry["gauge_alignment"] = check
            if not check["ok"]:
                for ch in staged:
                    gauges_ok[ch] = f"not reconstructed: PARO1 and PARO2 disagree by {check['lags']} s"
                staged = {}
        for ch, pieces in staged.items():
            out[ch].extend(pieces)
        out["SCPR"].append((t_scpr, srun["payloads"], srun["markers"]))
        entry["SCPR"] = "reconstructed" + (f" (instrument clock across {len(outages)} outage(s))" if len(outages) else "")
        entry.update(gauges_ok)
        report.append(entry)
    return out, report
