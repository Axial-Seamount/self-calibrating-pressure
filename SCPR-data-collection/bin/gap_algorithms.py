# Vendored unchanged from coszo-hub/absolute-seafloor-pressure
# PREST-data-collection/bin/gap_algorithms.py @ e2c34f3, so SCPR MiniSEED uses the
# same gap detection and trace segmentation as PREST/VEL3D. Keep in sync by copying.
"""
gap_algorithms.py — shared gap-detection interface for the COSZO PREST pipeline.

Two pure functions return the same GapResult shape so the cron pipeline,
testk smoke-test, and offline plotting tools can swap algorithms via a single
config flag (`gap_algo` in run_prest.txt).

    detect_gaps_legacy  : robust median Δt + adaptive multiplier × sp threshold
                          (mirrors OOI_data_request_and_convert_mseed.py
                          ~lines 460–700 prior to the migration)
    detect_gaps_anomaly : OLS Δt_true + integer-step reconstruction +
                          true_missing = n_ideal − n_points (Option A:
                          split only when true_missing > 0)

Lightweight — pure numpy. No matplotlib import, safe for the resource-
constrained VM.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class GapResult:
    sr: float                       # 1 / sp — lands in MiniSEED header
    sp: float                       # day-global sample interval (s)
    gap_idx: np.ndarray             # indices into dt_all where a gap starts
    segment_splits: list[int]       # split points for np.split on the timestamp array
    is_full: bool                   # wall-clock completeness flag
    n_gaps: int                     # post-correction gap count (anomaly: 0 if true_missing==0)
    n_segments: int                 # number of MiniSEED segments produced
    diagnostics: dict[str, Any] = field(default_factory=dict)


# ════════════════════════════════════════════════════════════════════════════
# Legacy — robust median + adaptive threshold
# ════════════════════════════════════════════════════════════════════════════
def detect_gaps_legacy(
    t_sec: np.ndarray,
    sp_nominal: float,
    req_duration: float | None = None,
) -> GapResult:
    """
    Pipeline-faithful legacy detector.

    Sample period:  median(Δt) after stripping the upper 10 % gap tail.
    Threshold:      multiplier × sp, where multiplier ∈ {2.0…4.0} depending on
                    sample-period regime and wall-clock completeness.

    `req_duration` is the requested window length in seconds (e.g. 86400 for a
    24 h pipeline run). If omitted, falls back to the actual t_sec span — only
    appropriate for offline tools where the request window equals the data span.
    """
    n = int(len(t_sec))
    if n < 2:
        raise ValueError(f"Need ≥2 samples, got {n}")

    dt_all = np.diff(t_sec)
    dt_pos = dt_all[np.isfinite(dt_all) & (dt_all > 0.0)]

    if dt_pos.size == 0:
        sp = float("nan")
    else:
        gap_cut = np.percentile(dt_pos, 90.0)
        dt_clean = dt_pos[dt_pos <= gap_cut]
        if dt_clean.size == 0:
            dt_clean = dt_pos
        sp = float(np.median(dt_clean))

    if not np.isfinite(sp) or sp <= 0:
        raise ValueError(f"Non-positive sample period: {sp}")

    sr = 1.0 / sp

    # Wall-clock completeness — pipeline uses requested window duration
    if req_duration is None:
        req_duration = float(t_sec[-1] - t_sec[0])
    expected_npts = int(round(req_duration / sp)) + 1
    tol = max(5, int(0.001 * expected_npts))
    is_full = abs(n - expected_npts) <= tol

    # Adaptive multiplier — same breakpoints as the pre-migration pipeline
    if sp >= 10.0:
        multiplier = 4.0 if is_full else 3.0
    elif sp >= 0.5:
        multiplier = 3.5 if is_full else 2.5
    else:
        multiplier = 2.5 if is_full else 2.0
    gap_threshold = multiplier * sp

    gap_idx = np.where(dt_all > gap_threshold)[0]
    segment_splits = (gap_idx + 1).tolist()

    if sp > 0 and len(gap_idx) > 0:
        gap_total_missing_est = int(sum(
            round(float(dt_all[gi]) / sp) - 1 for gi in gap_idx
        ))
    else:
        gap_total_missing_est = 0

    return GapResult(
        sr=sr,
        sp=sp,
        gap_idx=gap_idx,
        segment_splits=segment_splits,
        is_full=is_full,
        n_gaps=int(len(gap_idx)),
        n_segments=int(len(segment_splits) + 1),
        diagnostics={
            "algorithm": "legacy",
            "multiplier": multiplier,
            "gap_threshold": gap_threshold,
            "expected_npts": expected_npts,
            "gap_total_missing_est": gap_total_missing_est,
            "dt_all": dt_all,
        },
    )


# ════════════════════════════════════════════════════════════════════════════
# Anomaly — OLS Δt_true + integer-step + true_missing (Option A)
# ════════════════════════════════════════════════════════════════════════════
def detect_gaps_anomaly(
    t_sec: np.ndarray,
    sp_nominal: float,
    req_duration: float | None = None,    # accepted for interface parity; unused
) -> GapResult:
    """
    Integer-step reconstruction with OLS rate fit.

    Steps:
        Δt'        = diff(t_sec)
        Δt_FG      = median(Δt')
        Δi_int     = round(Δt' / Δt_FG)
        i_j        = cumsum(Δi_int)  (with i_0 = 0)
        n_ideal    = i_j[-1] + 1
        Δt_true    = polyfit(i_j, t_sec, 1).slope
        true_missing = max(0, n_ideal − n_points)

    Option A splitting: split MiniSEED only when true_missing > 0. Jitter-only
    days produce a single segment spanning the full window.
    """
    n = int(len(t_sec))
    if n < 2:
        raise ValueError(f"Need ≥2 samples, got {n}")

    dt_prime = np.diff(t_sec)
    dt_FG = float(np.median(dt_prime))
    if dt_FG <= 0:
        raise ValueError(f"Non-positive first-guess interval: {dt_FG}")

    delta_i_float = dt_prime / dt_FG
    delta_i_int   = np.round(delta_i_float).astype(np.int64)

    gap_positions_raw = np.where(delta_i_int > 1)[0]
    n_gaps_raw = int(gap_positions_raw.size)

    i_j = np.empty(n, dtype=np.int64)
    i_j[0] = 0
    np.cumsum(delta_i_int, out=i_j[1:])
    n_ideal = int(i_j[-1]) + 1

    true_missing = max(0, n_ideal - n)

    # OLS slope = Δt_true; intercept ignored (segment starttime stays as raw t_sec[0])
    slope, _ = np.polyfit(i_j.astype(float), t_sec, 1)
    dt_true = float(slope)
    if dt_true <= 0:
        raise ValueError(f"Non-positive OLS Δt_true: {dt_true}")

    sp = dt_true
    sr = 1.0 / dt_true

    # Jitter / instability diagnostics
    t_fit = i_j.astype(float) * dt_true   # intercept ≈ 0 in the t_sec frame
    e     = t_sec - (t_sec[0] + (t_fit - t_fit[0]))   # residuals re. observed start
    f     = e / dt_true
    fmax  = float(np.max(np.abs(f)))
    sigma_ms = float(np.std(1000.0 * e, ddof=0))
    jitter_unstable = bool(fmax > 0.4)

    # Wall-clock completeness — under the new algorithm "full" means n_ideal==n
    is_full = (true_missing == 0)

    # Option A — split only when true_missing > 0
    if true_missing == 0:
        gap_idx = np.array([], dtype=np.int64)
        segment_splits: list[int] = []
        n_gaps = 0
    else:
        gap_idx = gap_positions_raw
        segment_splits = (gap_idx + 1).tolist()
        n_gaps = n_gaps_raw

    return GapResult(
        sr=sr,
        sp=sp,
        gap_idx=gap_idx,
        segment_splits=segment_splits,
        is_full=is_full,
        n_gaps=n_gaps,
        n_segments=int(len(segment_splits) + 1),
        diagnostics={
            "algorithm": "anomaly",
            "dt_true": dt_true,
            "dt_FG": dt_FG,
            "n_ideal": n_ideal,
            "true_missing": true_missing,
            "n_gaps_raw": n_gaps_raw,
            "frac_maxabs": fmax,
            "jitter_unstable": jitter_unstable,
            "sigma_ms": sigma_ms,
            "dt_all": dt_prime,
        },
    )


# ════════════════════════════════════════════════════════════════════════════
# Duplicate records — OOI double ingestion
# ════════════════════════════════════════════════════════════════════════════
# Some OOI VEL3D-C periods hold every record twice: identical values (and
# ensemble counter) with the second copy offset by a drifting ~0.6–0.8 s, so a
# day has up to 2× the samples and the copies interleave. Left in, they shred
# the timing analysis and the MiniSEED. Only applied when a day has clearly
# more samples than its nominal rate allows, so normal days are untouched.
def excess_samples(t_sec, sp_nominal, excess=1.02):
    """True if the day holds clearly more samples than its nominal rate allows."""
    t_sec = np.asarray(t_sec, dtype=float)
    return bool(len(t_sec) > 1 and sp_nominal and
                len(t_sec) > excess * ((t_sec[-1] - t_sec[0]) / sp_nominal + 1))


def invalid_clock_mask(t_abs, internal_ts, sp_nominal, tol_s=1.0, min_good=0.4,
                       phase_tol=0.2):
    """Keep-mask dropping records whose instrument clock disagrees with time.

    In the OOI VEL3D-C excess-sample periods the extra records (duplicate
    copies, or foreign records interleaved from another time) carry a garbage
    internal_timestamp (unset, ~1e9 s off) and sit OFF the regular sample grid,
    while genuine records have internal_timestamp == time. A few genuine
    records also have a garbage clock but sit ON the grid in an otherwise
    empty slot — those are kept. Apply only on excess-sample days. Returns
    (keep, n_dropped); no-op if fewer than `min_good` of the records agree
    (then the instrument clock isn't a usable reference for this day).
    """
    t_abs = np.asarray(t_abs, dtype=float)
    d = np.abs(t_abs - np.asarray(internal_ts, dtype=float))
    good = np.isfinite(d) & (d <= tol_s)
    if good.sum() < min_good * len(d) or good.all():
        return np.ones(len(d), dtype=bool), 0
    sp = float(sp_nominal)
    g = t_abs[good]
    phase0 = np.angle(np.mean(np.exp(2j * np.pi * (g % sp) / sp)))   # grid phase
    keep = good.copy()
    bad = np.flatnonzero(~good)
    tb = t_abs[bad]
    ph = np.angle(np.exp(1j * (2 * np.pi * (tb % sp) / sp - phase0))) / (2 * np.pi)
    j = np.clip(np.searchsorted(g, tb), 1, len(g) - 1)
    nearest = np.minimum(np.abs(tb - g[j - 1]), np.abs(g[j] - tb))
    rescue = (np.abs(ph) <= phase_tol) & (nearest >= 0.5 * sp)   # on grid, empty slot
    keep[bad[rescue]] = True
    return keep, int((~keep).sum())


def duplicate_mask(t_sec, cols, sp_nominal, window_s=30.0, excess=1.02, prefer=None):
    """Boolean keep-mask dropping duplicated records.

    A sample duplicates another if ALL its values in `cols` (list of equal-
    length arrays, e.g. the three velocity components) are identical and the
    two are within `window_s`. Of each pair, keep the original: if `prefer`
    is given (e.g. |time − instrument internal_timestamp|; the duplicate copy
    carries a garbage internal timestamp) keep the copy with the smaller
    value; otherwise keep the copy whose time phase (t mod sp) matches the
    day's regular stream (circular-mean phase of non-duplicated samples). Returns
    (keep, n_dropped). No-op unless n > excess · (span / sp_nominal + 1).
    """
    t_sec = np.asarray(t_sec, dtype=float)
    n = len(t_sec)
    keep = np.ones(n, dtype=bool)
    if n < 2 or not cols or not sp_nominal:
        return keep, 0
    if n <= excess * ((t_sec[-1] - t_sec[0]) / sp_nominal + 1):
        return keep, 0
    vals = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    _, inv = np.unique(vals, axis=0, return_inverse=True)
    inv = inv.ravel()
    order = np.lexsort((t_sec, inv))                 # group by key, then time
    same = (inv[order][1:] == inv[order][:-1]) & \
           (np.diff(t_sec[order]) <= window_s)
    pairs = [(order[i], order[i + 1]) for i in np.flatnonzero(same)]
    if not pairs:
        return keep, 0
    in_pair = np.zeros(n, dtype=bool)
    for a, b in pairs:
        in_pair[a] = in_pair[b] = True
    ref = t_sec[~in_pair] if (~in_pair).sum() >= 10 else t_sec
    ang = 2 * np.pi * (ref % sp_nominal) / sp_nominal
    phase0 = np.angle(np.mean(np.exp(1j * ang)))    # dominant stream phase

    def off(i):
        d = 2 * np.pi * (t_sec[i] % sp_nominal) / sp_nominal - phase0
        return abs(np.angle(np.exp(1j * d)))

    for a, b in pairs:                              # a is the earlier copy
        if not (keep[a] and keep[b]):
            continue
        if prefer is not None and prefer[a] != prefer[b]:
            keep[b if prefer[a] < prefer[b] else a] = False
        else:
            keep[b if off(a) <= off(b) else a] = False
    return keep, int((~keep).sum())


# ════════════════════════════════════════════════════════════════════════════
# MiniSEED segmentation — where to start a new trace
# ════════════════════════════════════════════════════════════════════════════
# A MiniSEED trace stores ONE start time + a rate; sample k sits at
# start + k·sp. "gaps" (the original behaviour, kept selectable) starts a new
# trace only at detected gaps, anchored on the segment's first timestamp. A
# burst of late-delivered samples (stamped on arrival) or a clock step inside a
# segment then shifts every later sample — measured up to ~2 s at 1 Hz on real
# RCA days. "timing" also starts a new trace wherever the recorded timestamps
# leave the regular grid by more than `tol`·sp, so every written sample stays
# within ½ sample of the time OOI recorded. Exception: a lone glitched
# timestamp whose neighbours stay on the grid is kept in its grid slot (the
# neighbours fix its position better than its own bad stamp).
SEGMENTING_MODES = ("timing", "gaps")


def timing_segments(t_sec, sp, base_splits, tol=0.5, level_window=25, glitch_run=3):
    """Refine gap-based segments so each trace follows its timestamps.

    Returns (splits, starts): `splits` is a superset of `base_splits` for
    np.split on the timestamp array; `starts[i]` is the start time of segment
    i in the t_sec frame — the segment's robust grid level (median of
    t_k − k·sp over its leading samples), not just its first timestamp.
    """
    t_sec = np.asarray(t_sec, dtype=float)
    bounds = [0] + list(base_splits) + [len(t_sec)]
    thr = tol * sp
    splits, starts = [], []
    for a, b in zip(bounds[:-1], bounds[1:]):
        r = t_sec[a:b] - np.arange(b - a) * sp          # constant (± jitter) on a clean run
        s = 0
        while s < b - a:
            # Level of the run starting at s: median over the leading samples
            # that agree with r[s] (within thr). A short shifted run (e.g. 8
            # samples stamped one period early) gets its own level instead of
            # being judged against whatever follows it.
            w = r[s:s + level_window]
            off = np.flatnonzero(np.abs(w - w[0]) > thr)
            lead = w[:off[0]] if len(off) else w
            level = float(np.median(lead))
            if len(lead) == 1 and len(w) > 1 and abs(w[1] - level) > thr:
                # isolated leading sample (e.g. one burst sample after an
                # outage): own one-sample trace at its timestamp
                e = s + 1
                seg_start = float(r[s])
            else:
                # Scan forward in growing chunks (not the whole remainder each
                # time — that is quadratic on days with many traces).
                e = b - a
                lo, chunk, found = s + 1, 256, False
                while lo < b - a and not found:
                    hi = min(b - a, lo + chunk)
                    bad = np.flatnonzero(np.abs(r[lo:hi] - level) > thr) + lo
                    for k in bad:
                        nxt = np.abs(r[k + 1:k + 1 + glitch_run] - level) > thr
                        if len(nxt) == glitch_run and not nxt.any():
                            continue                     # lone glitch: keep on grid
                        e, found = int(k), True
                        break
                    lo, chunk = hi, min(chunk * 4, 1 << 20)
                seg_start = level
            if a + s > 0:
                splits.append(a + s)
            starts.append(seg_start + s * sp)            # back to t_sec frame
            s = e
    return splits, starts


def mseed_segments(mode, t_sec, gap_result, tol=0.5):
    """(splits, starts) for the MiniSEED writer. mode='gaps' reproduces the
    original behaviour exactly (gap splits, start = first timestamp)."""
    if mode == "gaps":
        splits = list(gap_result.segment_splits)
        return splits, [float(t_sec[i]) for i in [0] + splits]
    if mode == "timing":
        return timing_segments(t_sec, gap_result.sp, gap_result.segment_splits, tol=tol)
    raise ValueError(f"Unknown mseed segmenting {mode!r}; expected one of {SEGMENTING_MODES}")


# ════════════════════════════════════════════════════════════════════════════
# Selector
# ════════════════════════════════════════════════════════════════════════════
ALGORITHMS = {
    "legacy":  detect_gaps_legacy,
    "anomaly": detect_gaps_anomaly,
}

# Below this sample count the OLS fit and integer-step reconstruction become
# numerically suspect; anomaly mode falls back to legacy automatically. 100 is
# 25 minutes of 15 s data and 100 seconds of 1 Hz data — both safely larger
# than any backfill chunk we'd ever care to run through anomaly.
MIN_N_FOR_ANOMALY = 100


def detect_gaps(algo: str, t_sec: np.ndarray, sp_nominal: float, **kwargs) -> GapResult:
    """Dispatch helper. `algo` ∈ {"legacy", "anomaly"}.

    Short-window guardrail: if `algo == "anomaly"` and len(t_sec) < MIN_N_FOR_ANOMALY,
    silently fall back to legacy with a warning. The actual algorithm used is
    reflected in `result.diagnostics["algorithm"]` so callers (e.g. the CSV
    writer) can record what really ran.

    Extra kwargs (e.g. `req_duration`) are forwarded to the algorithm
    function. Both detectors accept `req_duration`; anomaly ignores it.
    """
    if algo == "anomaly" and len(t_sec) < MIN_N_FOR_ANOMALY:
        print(
            f"WARNING: n={len(t_sec)} below MIN_N_FOR_ANOMALY={MIN_N_FOR_ANOMALY}; "
            f"falling back to legacy for this window"
        )
        algo = "legacy"

    fn = ALGORITHMS.get(algo)
    if fn is None:
        raise ValueError(
            f"Unknown gap_algo {algo!r}; expected one of {list(ALGORITHMS)}"
        )
    return fn(t_sec, sp_nominal, **kwargs)
