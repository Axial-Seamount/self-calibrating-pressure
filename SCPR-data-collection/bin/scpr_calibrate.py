#!/usr/bin/env python3
"""SCPR calibration from MiniSEED: extract calibration offsets, then drift-correct.

Works on the SDS tree written by scpr_mseed.py (OO.AXCC1, loc 40 = gauge 1,
loc 41 = gauge 2).

  extract  Find calibrations in LDO.40/LDO.41 and estimate each gauge's offset
           (gauge minus internal barometer, LDI.40) with the same algorithm as
           extract_scpr_calibration.py: t0 = minimum within 100 s of the downward
           threshold crossing, 12 x 100-s bins, one-sigma trim, mean of bin means.
           Results are merged into calibration/calibrations.csv (committed).

  apply    Build each gauge's drift curve from those offsets and write
           drift-corrected LDO MiniSEED (data quality Q) to a separate SDS tree.
           Drift is relative (SCPR practice): delta_k = d_k - d_ref, and
           corrected = P_obs - delta(t). Before differencing, each offset gets the
           dead-weight-tester thermal term d_k + d_k * 9e-6 * (T_PRT,k - T_ref)
           (piston + cylinder area expansion; T_ref = median event PRT). Default curve: weighted fit of
           a + b*t + c*exp(-t/tau) with tau grid-searched (Wilcock et al. 2021);
           a linear trend if the exponential is not significant; plain linear
           interpolation between events if fewer than 5 events or --method linear.
           The curve is held before the first event. After the last event both methods
           extrapolate a weighted line through the last 4 events (best in a
           leave-future-out hindcast), with uncertainty sqrt(40^2 + (170 Pa/yr * dt)^2);
           correction stops when that exceeds --max-drift-sigma (100 Pa, ~200 days).
           Drift is referenced to the smooth curve at the first trusted event.
           Calibrations (LEF = 1 or pressure below the threshold) are masked from
           60 s before to --settle-seconds after. Raw LDO is never modified.
           Each corrected day gets output/manifest/scpr_corrected_DATE.json with the
           table version (hash) used and whether it is provisional (after the last
           calibration, so the curve there is extrapolated).

  model    Write calibration/drift_model.json and drift.png from the table alone.

Requires numpy and ObsPy (conda env ooi_env).
    python scpr_calibrate.py extract --start 2020-05-04 --end 2020-05-04
    python scpr_calibrate.py apply   --start 2020-05-01 --end 2020-05-31
    python scpr_calibrate.py model
Most users want scpr.py, which wraps download, conversion and correction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date, datetime, timedelta, timezone
import math
from pathlib import Path
import sys

import numpy as np
from obspy import Stream, UTCDateTime, read
from obspy.clients.filesystem.sds import Client

import extract_scpr_calibration as cal
import scpr_paths as paths


NETWORK, STATION = "OO", "AXCC1"
GAUGE_LOC = {"PARO1": "40", "PARO2": "41"}
DEFAULT_MSEED = paths.MSEED
# Committed with the code: the table and the drift model are all a user needs to correct a day.
MODEL_DIR = paths.CALIBRATION
DEFAULT_TABLE = MODEL_DIR / "calibrations.csv"
DEFAULT_CORRECTED = paths.CORRECTED
DEFAULT_REPORT = paths.REPORT
DEFAULT_MANIFEST = paths.MANIFEST
# Piston + cylinder thermal expansion of the DWT effective area (Supplementary Materials).
DWT_AREA_ALPHA_PER_C = 4.5e-6 + 4.5e-6
# Deployment-day calibrations (2018-07-06) were not in thermal equilibrium; excluded by default.
DEFAULT_EXCLUDE_BEFORE = "2018-07-07"
MIN_EVENTS_FOR_FIT = 5
# Event QC: good calibrations have SE <= ~85 Pa (2018-2026); larger means the
# dead-weight tester was not floating steadily during [t0, t0 + 1200 s).
MAX_SE_PA = 200.0
# Event-to-event scatter not captured by an event's own SE (which only measures noise
# inside its 20-min window; 6-20 Pa). Calibrations scatter 30-45 Pa about smooth drift,
# so each event's uncertainty is sqrt(SE^2 + EVENT_SCATTER_PA^2) for weighting.
EVENT_SCATTER_PA = 30.0
# After the last calibration (chosen by a leave-future-out hindcast, 2026-11-12: predict
# each calibration from the earlier ones; rms error PARO1/PARO2: hold flat 112/103 Pa,
# continue the fitted curve 116/58 Pa, weighted line through the last 4 events 66/54 Pa):
# extrapolate a weighted straight line through the last TAIL_EVENTS events, tapering the
# step from the curve to that line over the last gap so the drift stays continuous at the
# last event (costs ~6 Pa in the hindcast; a step would look like a pressure signal).
# Reaching back further (6+ events, or >= 1 year in 2018-19) pulls in the fast early drift
# and roughly doubles the PARO1 error.
TAIL_EVENTS = 4
# Drift uncertainty: DRIFT_SIGMA0_PA inside the record (event scatter about a smooth
# curve), growing as sqrt(s0^2 + (rate * days since last event)^2) afterwards (rate
# fitted to the hindcast errors). Correction stops when it exceeds MAX_DRIFT_SIGMA_PA
# (~1 cm of water; ~200 days) or at --max-extrapolation-days, whichever comes first.
DRIFT_SIGMA0_PA, DRIFT_SIGMA_RATE_PA_PER_YR, MAX_DRIFT_SIGMA_PA = 40.0, 170.0, 100.0
DRIFT_MODEL_VERSION = 3          # bump when the drift logic changes: local corrections get redone
TAU_GRID_DAYS = np.geomspace(2.0, 3000.0, 240)
MASK_BEFORE_S = 60.0
PA_PER_PSI = cal.KPA_PER_PSI * 1000.0
THRESHOLD_PA = 15400e3                      # same 15400 kPa threshold as the .dat pipeline
TABLE_COLUMNS = ["gauge", "location", "t0_utc", "offset_pa", "offset_psi", "bin_means_sd_pa", "se_pa",
                 "native_samples", "retained_samples", "max_abs_barometer_dt_s", "barometer_mean_pa",
                 "prt_mean_degc", "flag_fraction", "status",
                 # second ("masked") method, see masked_offset()
                 "offset2_pa", "se2_pa", "offset2_trim_pa", "bins2", "window2_s", "masked2_s",
                 "slope2_pa_per_min", "split_half2_pa", "status2"]

# Second offset method (masked): mechanical activity logged by the controller during
# the calibration is masked, the first SKIP2_S of settling is skipped, and each
# 100-s bin is summarised by its median (bins weighted by retained samples).
RAW_ROOT = paths.RAW
DISTURBANCE_TEXT = ("Re-Spin", "Spin motor", "Belt", "Charge Valve Open", "Charge Valve Close",
                    "Opening Valve", "Closing Valve", "Valve 1 closed")
DISTURB_GROUP_S, DISTURB_BEFORE_S, DISTURB_AFTER_S = 30.0, 5.0, 60.0
SKIP2_S, END_MARGIN2_S, MAX_END2_S, BIN2_S, MIN_BIN_FRACTION = 300.0, 60.0, 1500.0, 100.0, 0.5


# ── MiniSEED access ─────────────────────────────────────────────────────────
def load(client: Client, loc: str, cha: str, start: UTCDateTime, end: UTCDateTime) -> tuple[np.ndarray, np.ndarray]:
    """Sample times (epoch s) and values of one channel, time-sorted across traces/days."""
    # merge=False: each day has its own fitted rate, so ObsPy cannot merge adjacent days.
    st: Stream = client.get_waveforms(NETWORK, STATION, loc, cha, start, end, merge=False)
    if not len(st):
        return np.array([]), np.array([])
    st.sort(keys=["starttime"])
    t = np.concatenate([tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate for tr in st])
    v = np.concatenate([tr.data.astype(float) for tr in st])
    order = np.argsort(t, kind="stable")
    t, v = t[order], v[order]
    keep = np.ones(len(t), dtype=bool)
    keep[1:] = np.diff(t) > 1e-6            # traces from adjacent day files may share an edge sample
    return t[keep], v[keep]


def to_samples(t: np.ndarray, pressure_pa: np.ndarray, temp_t: np.ndarray, temp: np.ndarray) -> list[cal.Sample]:
    """extract_scpr_calibration.Sample list (pressure in psi, nearest temperature)."""
    if len(temp_t):
        j = np.clip(np.searchsorted(temp_t, t), 1, len(temp_t) - 1)
        j -= (t - temp_t[j - 1]) <= (temp_t[j] - t)
        tt = temp[j]
    else:
        tt = np.full(len(t), math.nan)
    return [cal.Sample(datetime.fromtimestamp(a, timezone.utc), p / PA_PER_PSI, b)
            for a, p, b in zip(t.tolist(), pressure_pa.tolist(), tt.tolist())]


def window_mean(t: np.ndarray, v: np.ndarray, start: float, seconds: float) -> float:
    sel = (t >= start) & (t < start + seconds)
    return float(np.mean(v[sel])) if sel.any() else math.nan


# ── second offset method ────────────────────────────────────────────────────
def log_disturbances(day: date) -> np.ndarray | None:
    """(start, end) epochs of mechanical activity in the controller log, padded.

    Lines naming a re-spin, the belt or the charge valve are grouped when less
    than DISTURB_GROUP_S apart, then padded by DISTURB_BEFORE_S / DISTURB_AFTER_S.
    None if the raw SCPR log for the day is not on disk.
    """
    from scpr_mseed import read_log, to_epoch
    paths = sorted(p for d in (day, day + timedelta(days=1))
                   for p in (RAW_ROOT / str(d)).glob("SCPRAA301_SCPR_*.dat"))
    if not any(str(day) in p.parent.name for p in paths):
        return None
    stamps = [stamp for stamp, payload in read_log(paths) if any(k in payload for k in DISTURBANCE_TEXT)]
    if not stamps:
        return np.empty((0, 2))
    t = np.sort(to_epoch(stamps))
    breaks = np.flatnonzero(np.diff(t) > DISTURB_GROUP_S)
    starts = np.r_[t[0], t[breaks + 1]] - DISTURB_BEFORE_S
    ends = np.r_[t[breaks], t[-1]] + DISTURB_AFTER_S
    return np.column_stack([starts, ends])


def masked_offset(t: np.ndarray, d: np.ndarray, t0: float, cal_end: float, disturb: np.ndarray) -> dict:
    """Offset from gauge-minus-barometer samples d(t) with disturbances masked.

    Window [t0 + SKIP2_S, min(cal_end - END_MARGIN2_S, t0 + MAX_END2_S)); 100-s bins;
    a bin with < MIN_BIN_FRACTION of its seconds left is dropped; offset = mean of
    bin medians weighted by retained samples. QC: slope of bin medians (Pa/min) and
    the difference between the second and first half of the bins.
    """
    lo, hi = t0 + SKIP2_S, min(cal_end - END_MARGIN2_S, t0 + MAX_END2_S)
    sel = (t >= lo) & (t < hi)
    t, d = t[sel], d[sel]
    clean = ~in_intervals(t, disturb)
    out = {"window2_s": round(hi - lo, 1), "masked2_s": int((~clean).sum())}
    meds, trims, mids, counts = [], [], [], []
    for a in np.arange(lo, hi - BIN2_S / 2, BIN2_S):
        x = d[clean & (t >= a) & (t < a + BIN2_S)]
        if len(x) < MIN_BIN_FRACTION * min(BIN2_S, hi - a):
            continue
        keep = np.abs(x - x.mean()) <= x.std()
        meds.append(np.median(x))
        trims.append(x[keep].mean())
        mids.append((a + BIN2_S / 2 - t0) / 60.0)
        counts.append(len(x))
    if len(meds) < 3:
        return {**out, "bins2": len(meds), "status2": f"failed: only {len(meds)} usable bins"}
    w = np.array(counts, dtype=float)
    meds, trims, mids = map(np.array, (meds, trims, mids))
    mean = float(np.average(meds, weights=w))
    sd = float(np.sqrt(np.average((meds - mean) ** 2, weights=w)))
    half = len(meds) // 2
    return {**out, "bins2": len(meds), "offset2_pa": mean, "se2_pa": sd / math.sqrt(len(meds)),
            "offset2_trim_pa": float(np.average(trims, weights=w)),
            "slope2_pa_per_min": float(np.polyfit(mids, meds, 1)[0]),
            "split_half2_pa": float(meds[half:].mean() - meds[:half].mean()), "status2": "ok"}


def nearest(t_ref: np.ndarray, v_ref: np.ndarray, t: np.ndarray) -> np.ndarray:
    j = np.clip(np.searchsorted(t_ref, t), 1, len(t_ref) - 1)
    j -= (t - t_ref[j - 1]) <= (t_ref[j] - t)
    return v_ref[j]


def calibration_end(ft: np.ndarray, flag: np.ndarray, gt: np.ndarray, p: np.ndarray, t0: float) -> float:
    """End of the calibration containing t0: LEF back to 0, else pressure back above threshold."""
    after = ft > t0
    off = after & (flag == 0)
    if off.any() and (flag[(ft >= t0) & (ft < ft[off][0])] > 0).any():
        return float(ft[off][0])
    up = (gt > t0) & (p >= THRESHOLD_PA)
    return float(gt[up][0]) if up.any() else float(gt[-1])


# ── extract ─────────────────────────────────────────────────────────────────
def extract(args: argparse.Namespace) -> None:
    client = Client(str(args.mseed))
    day = args.start
    found: list[dict] = []
    empty_days = 0
    while day <= args.end:
        lo = UTCDateTime(day.isoformat())
        # Neighbouring days let an event near midnight be measured in full.
        start, end = lo - 86400, lo + 2 * 86400
        bt, baro = load(client, "40", "LDI", start, end)
        pt, prt = load(client, "40", "LKI", start, end)
        ft, flag = load(client, "40", "LEF", start, end)
        if not len(bt) or not (bt >= lo.timestamp).any() or not (bt < lo.timestamp + 86400).any():
            empty_days += 1
            day += timedelta(days=1)
            continue
        auxiliary = to_samples(bt, baro, pt, prt)
        disturb = log_disturbances(day)
        for gauge, loc in GAUGE_LOC.items():
            gt, pressure = load(client, loc, "LDO", start, end)
            if not len(gt):
                print(f"{day} {gauge}: no LDO.{loc} data")
                continue
            # LK1 only rides along as the sample temperature; offsets do not use it.
            samples = to_samples(gt, pressure, *load(client, loc, "LK1", start, end))
            origins, rejected = cal.detect_events(samples, THRESHOLD_PA / PA_PER_PSI)
            for origin in origins:
                t0 = samples[origin].time
                if t0.date() != day:
                    continue
                try:
                    _, _, s = cal.summarize_event(gauge, 1, samples, origin, auxiliary, None)
                    status = "ok"
                except ValueError as exc:   # coverage failure: keep a record, no offset
                    s, status = {"t0_utc": cal.iso(t0)}, f"failed: {exc}"
                t0s = t0.timestamp()
                row = {"gauge": gauge, "location": loc, "t0_utc": s["t0_utc"], "status": status,
                       "barometer_mean_pa": window_mean(bt, baro, t0s, cal.DURATION),
                       "prt_mean_degc": window_mean(pt, prt, t0s, cal.DURATION),
                       "flag_fraction": window_mean(ft, flag, t0s, cal.DURATION)}
                if status == "ok" and s["sasagawa_se_psi"] * PA_PER_PSI > MAX_SE_PA:
                    # e.g. 2021-05-12: the DWT piston did not float for 20 minutes and the
                    # gauge wandered by 50 kPa; SE ~3900 Pa against <= 85 Pa for good events.
                    status = (f"rejected: unstable, SE {s['sasagawa_se_psi'] * PA_PER_PSI:.0f} Pa "
                              f"> {MAX_SE_PA:.0f} Pa")
                    row["status"] = status
                if "offset_pa" not in row and "calibration_pressure_psi" in s:
                    row.update({"offset_psi": s["calibration_pressure_psi"],
                                "offset_pa": s["calibration_pressure_psi"] * PA_PER_PSI,
                                "bin_means_sd_pa": s["bin_means_sd_psi"] * PA_PER_PSI,
                                "se_pa": s["sasagawa_se_psi"] * PA_PER_PSI,
                                "native_samples": s["native_samples"], "retained_samples": s["retained_samples"],
                                "max_abs_barometer_dt_s": s["max_abs_barometer_dt_s"]})
                if disturb is None:
                    row["status2"] = "failed: raw SCPR log not in output/raw"
                else:
                    row.update(masked_offset(gt, pressure - nearest(bt, baro, gt), t0s,
                                             calibration_end(ft, flag, gt, pressure, t0s), disturb))
                found.append(row)
                print(f"{gauge} t0={row['t0_utc']}: {status}"
                      + (f", offset {row['offset_pa']:.1f} Pa ({row['offset_psi']:.9f} psi), "
                         f"flag {row['flag_fraction']:.2f}" if "offset_pa" in row else ""))
            for r in rejected:
                if r["t0_utc"][:10] == str(day):
                    print(f"{gauge} candidate {r['t0_utc']} rejected: {r['reason']}")
        day += timedelta(days=1)
    if empty_days:
        print(f"{empty_days} day(s) without LDI.40 barometer data were skipped")
    merge_table(args.table, found, args.start, args.end)


def read_table(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def merge_table(path: Path, rows: list[dict], first: date, last: date) -> None:
    """Replace this date range's rows in the table; keep everything else."""
    def in_range(r: dict) -> bool:
        return first <= date.fromisoformat(r["t0_utc"][:10]) <= last
    kept = [r for r in read_table(path) if not in_range(r)]
    table = sorted(kept + rows, key=lambda r: (r["t0_utc"], r["gauge"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".csv.part")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, TABLE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in table:
            writer.writerow({k: (format(v, ".10g") if isinstance(v, float) else v) for k, v in r.items()})
    tmp.replace(path)
    print(f"{len(rows)} calibration rows for {first}..{last}; table now {len(table)} rows: {path}")


# ── drift model ─────────────────────────────────────────────────────────────
class DriftModel:
    """Drift of one gauge (Pa), relative to the smooth curve at its first trusted calibration.

    Inside the record: the weighted parametric fit (method "fit") or interpolation through
    every event (method "linear"). Before the first event: held. After the last event:
    a weighted line through the last TAIL_EVENTS events (see the constants), with a
    per-time uncertainty sigma_drift(t); NaN once that exceeds max_sigma or past
    max_extrapolation_days.
    """

    def __init__(self, gauge: str, rows: list[dict], method: str, max_extrapolation_days: float,
                 thermal: bool, offsets: str = "glenn", event_scatter: float = EVENT_SCATTER_PA,
                 tail_events: int = TAIL_EVENTS, max_sigma: float = MAX_DRIFT_SIGMA_PA) -> None:
        self.gauge = gauge
        self.t = np.array([UTCDateTime(r["t0_utc"]).timestamp for r in rows])
        d = np.array([float(r[OFFSET_COLUMNS[offsets][0]]) for r in rows])
        self.se = np.array([float(r[OFFSET_COLUMNS[offsets][1]]) for r in rows])
        self.sigma = np.sqrt(self.se**2 + event_scatter**2)       # used for weights and error bars
        prt = np.array([float(r["prt_mean_degc"]) for r in rows])
        # DWT area grows with temperature, lowering the generated pressure; add it back.
        # Referenced to the median PRT: the first events (days after deployment) were not
        # yet at the ambient ~7.2 degC (2018-07-09: 7.86 degC).
        self.thermal = d * DWT_AREA_ALPHA_PER_C * (prt - np.median(prt)) if thermal else np.zeros(len(d))
        self.d = d
        self.y = d + self.thermal - d[0]                      # adjusted offsets, relative to event 1
        self.days = (self.t - self.t[0]) / 86400.0
        n = len(rows)
        smooth = self._fit() if n >= MIN_EVENTS_FOR_FIT else None
        self.kind, self.coef, self.tau = "linear interpolation", None, None
        if method == "fit" and smooth:
            self.kind, self.coef, self.tau = smooth

        # Tail: weighted line through the last events, x = days since the last event.
        i0 = max(0, n - tail_events)
        self.tail_events = n - i0
        if self.tail_events >= 2:
            w = 1.0 / self.sigma[i0:]
            X = np.column_stack([np.ones(self.tail_events), self.days[i0:] - self.days[-1]])
            self.tail = tuple(float(v) for v in np.linalg.lstsq(X * w[:, None], self.y[i0:] * w, rcond=None)[0])
        else:
            self.tail = (float(self.y[-1]), 0.0)
        self.taper_days = float(self.days[-1] - self.days[-2]) if n >= 2 else 0.0

        # Reference: the smooth curve at the first event, the same for both methods
        # (the raw first event, 2018-07-09, sits ~250 Pa off the curve).
        if smooth:
            self.ref = float(self._param(smooth[1], smooth[2], np.array([0.0]))[0])
        elif n >= 2:
            w = 1.0 / self.sigma
            X = np.column_stack([np.ones(n), self.days])
            self.ref = float(np.linalg.lstsq(X * w[:, None], self.y * w, rcond=None)[0][0])
        else:
            self.ref = float(self.y[0])
        self.residual = self.y - self._interior(self.t)
        self.rms = float(np.sqrt(np.mean(self.residual**2)))
        dt_sigma = (math.sqrt(max_sigma**2 - DRIFT_SIGMA0_PA**2) / DRIFT_SIGMA_RATE_PA_PER_YR * 365.25
                    if max_sigma > DRIFT_SIGMA0_PA else 0.0)
        self.t_last_valid = self.t[-1] + min(max_extrapolation_days, dt_sigma) * 86400.0

    def _fit(self) -> tuple[str, np.ndarray, float | None]:
        """Weighted a + b*t + c*exp(-t/tau), tau by grid search; a + b*t if c is not significant."""
        w = 1.0 / self.sigma
        best = None
        for tau in TAU_GRID_DAYS:
            X = np.column_stack([np.ones_like(self.days), self.days, np.exp(-self.days / tau)])
            coef, *_ = np.linalg.lstsq(X * w[:, None], self.y * w, rcond=None)
            sse = float(np.sum(((X @ coef - self.y) * w) ** 2))
            if best is None or sse < best[0]:
                best = (sse, tau, coef, X)
        sse, tau, coef, X = best
        n = len(self.y)
        cov = np.linalg.pinv((X * w[:, None]).T @ (X * w[:, None])) * sse / max(n - 3, 1)
        if abs(coef[2]) >= 2.0 * math.sqrt(max(cov[2, 2], 0.0)):
            return "a + b*t + c*exp(-t/tau)", coef, float(tau)
        X = np.column_stack([np.ones_like(self.days), self.days])      # exponential not significant
        coef, *_ = np.linalg.lstsq(X * w[:, None], self.y * w, rcond=None)
        return "a + b*t (exponential not significant)", np.append(coef, 0.0), None

    @staticmethod
    def _param(coef: np.ndarray, tau: float | None, days: np.ndarray) -> np.ndarray:
        a, b, c = coef
        return a + b * days + (c * np.exp(-days / tau) if tau else 0.0)

    def _interior(self, t: np.ndarray) -> np.ndarray:
        """The curve inside the record (interp or fit), relative to event 1."""
        if self.coef is None:
            return np.interp(t, self.t, self.y)
        return self._param(self.coef, self.tau, (t - self.t[0]) / 86400.0)

    def _curve(self, t: np.ndarray) -> np.ndarray:
        out = self._interior(t)
        late = t > self.t[-1]
        if late.any():
            dt = (t[late] - self.t[-1]) / 86400.0
            a, b = self.tail
            step = float(self._interior(self.t[-1:])[0]) - a          # curve minus line at the last event
            taper = np.clip(1.0 - dt / self.taper_days, 0.0, 1.0) if self.taper_days > 0 else 0.0
            out[late] = a + b * dt + step * taper
        return out

    def delta(self, t: np.ndarray) -> np.ndarray:
        """Drift at t; held before the first event, tail line after the last, NaN past the limit."""
        out = self._curve(np.maximum(t, self.t[0])) - self.ref
        out[t > self.t_last_valid] = np.nan
        return out

    def sigma_drift(self, t: np.ndarray) -> np.ndarray:
        """1-sigma uncertainty of delta(t) (Pa)."""
        dt_yr = np.maximum(t - self.t[-1], 0.0) / 86400.0 / 365.25
        return np.sqrt(DRIFT_SIGMA0_PA**2 + (DRIFT_SIGMA_RATE_PA_PER_YR * dt_yr) ** 2)

    def report(self, rows: list[dict]) -> dict:
        chi2 = float(np.mean((self.residual / self.sigma) ** 2))
        return {"gauge": self.gauge, "model": self.kind, "tau_days": self.tau,
                "coefficients_pa_per_day": None if self.coef is None else [float(c) for c in self.coef],
                "rms_residual_pa": self.rms, "mean_squared_normalised_residual": chi2,
                "reference": {"epoch_utc": rows[0]["t0_utc"], "relative_pa": self.ref,
                              "definition": "smooth curve at the first trusted event"},
                "tail": {"events": self.tail_events, "slope_pa_per_year": self.tail[1] * 365.25,
                         "level_at_last_event_pa": self.tail[0] - self.ref, "taper_days": self.taper_days},
                "uncertainty": {"sigma0_pa": DRIFT_SIGMA0_PA, "rate_pa_per_year": DRIFT_SIGMA_RATE_PA_PER_YR,
                                "formula": "sqrt(sigma0^2 + (rate * years since last event)^2)"},
                "last_event_utc": rows[-1]["t0_utc"], "valid_until_utc": iso_s(self.t_last_valid),
                "events": [{"t0_utc": r["t0_utc"], "offset_pa": float(d), "se_pa": float(se), "sigma_pa": float(sg),
                            "prt_mean_degc": float(r["prt_mean_degc"]), "thermal_adjust_pa": float(th),
                            "drift_pa": float(y - self.ref), "residual_pa": float(res)}
                           for r, d, se, sg, th, y, res in zip(rows, self.d, self.se, self.sigma, self.thermal,
                                                               self.y, self.residual)]}


def hindcast(gauge: str, rows: list[dict], args: argparse.Namespace) -> dict:
    """Predict each event from the earlier ones (leave-future-out): the tail rule's track record."""
    errs, leads = [], []
    for k in range(MIN_EVENTS_FOR_FIT + 1, len(rows)):
        m = DriftModel(gauge, rows[:k], args.method, 1e9, not args.no_thermal,
                       getattr(args, "offsets", "glenn"), getattr(args, "event_scatter", EVENT_SCATTER_PA),
                       getattr(args, "tail_events", TAIL_EVENTS), math.inf)
        t_k = UTCDateTime(rows[k]["t0_utc"]).timestamp
        full = DriftModel(gauge, rows[:k + 1], args.method, 1e9, not args.no_thermal,
                          getattr(args, "offsets", "glenn"), getattr(args, "event_scatter", EVENT_SCATTER_PA))
        errs.append(float(m._curve(np.array([t_k]))[0] - full.y[k]))
        leads.append((t_k - m.t[-1]) / 86400.0)
    e, lead = np.array(errs), np.array(leads)
    if not len(e):
        return {}
    pred = np.sqrt(DRIFT_SIGMA0_PA**2 + (DRIFT_SIGMA_RATE_PA_PER_YR * lead / 365.25) ** 2)
    return {"trials": len(e), "lead_days": [float(lead.min()), float(lead.max())],
            "rms_error_pa": float(np.sqrt(np.mean(e**2))), "mean_error_pa": float(e.mean()),
            "mean_squared_error_over_sigma": float(np.mean((e / pred) ** 2))}


def iso_s(seconds: float) -> str:
    return UTCDateTime(seconds).strftime("%Y-%m-%dT%H:%M:%SZ")


# Offset method -> (offset column, SE column, status column)
OFFSET_COLUMNS = {"glenn": ("offset_pa", "se_pa", "status"), "masked": ("offset2_pa", "se2_pa", "status2")}


def trusted_rows(table: Path, gauge: str, exclude_before: date, offsets: str = "glenn") -> list[dict]:
    """Usable events: Glenn's status ok (QC), plus the chosen method's own status ok."""
    col, se, status = OFFSET_COLUMNS[offsets]
    rows = [r for r in read_table(table) if r["gauge"] == gauge and r["status"] == "ok"
            and r.get(status) == "ok" and r.get(col) and float(r["se_pa"]) <= MAX_SE_PA
            and date.fromisoformat(r["t0_utc"][:10]) >= exclude_before]
    return sorted(rows, key=lambda r: r["t0_utc"])


# ── apply ───────────────────────────────────────────────────────────────────
def mask_intervals(ft: np.ndarray, flag: np.ndarray, pt: np.ndarray, p: np.ndarray, settle: float) -> np.ndarray:
    """(start, end) epochs to drop: LEF = 1 or pressure below threshold, padded."""
    bad = np.sort(np.concatenate([ft[flag > 0], pt[p < THRESHOLD_PA]]))
    if not len(bad):
        return np.empty((0, 2))
    breaks = np.flatnonzero(np.diff(bad) > 60.0)
    starts = np.r_[bad[0], bad[breaks + 1]] - MASK_BEFORE_S
    ends = np.r_[bad[breaks], bad[-1]] + settle
    merged = [[starts[0], ends[0]]]
    for a, b in zip(starts[1:], ends[1:]):
        if a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return np.array(merged)


def in_intervals(t: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    if not len(intervals):
        return np.zeros(len(t), dtype=bool)
    i = np.searchsorted(intervals[:, 0], t, side="right") - 1
    return (i >= 0) & (t <= intervals[np.maximum(i, 0), 1])


def table_version(table: Path) -> str:
    """Short content hash of the calibration table; recorded with every corrected day."""
    return hashlib.sha256(table.read_bytes()).hexdigest()[:12] if table.exists() else "none"


def settings(args: argparse.Namespace) -> dict:
    return {"table_version": table_version(args.table), "method": args.method,
            "offsets": getattr(args, "offsets", "glenn"),
            "event_scatter_pa": getattr(args, "event_scatter", EVENT_SCATTER_PA),
            "tail_events": getattr(args, "tail_events", TAIL_EVENTS),
            "max_drift_sigma_pa": getattr(args, "max_drift_sigma", MAX_DRIFT_SIGMA_PA),
            "drift_model_version": DRIFT_MODEL_VERSION,
            "thermal_adjustment": not args.no_thermal, "exclude_before": str(args.exclude_before),
            "max_extrapolation_days": args.max_extrapolation_days, "settle_seconds": args.settle_seconds}


def build_models(args: argparse.Namespace, verbose: bool = True) -> tuple[dict[str, DriftModel], dict]:
    """Drift model per gauge from the whole table, plus the model report."""
    models, reports = {}, []
    for gauge in GAUGE_LOC:
        rows = trusted_rows(args.table, gauge, args.exclude_before, getattr(args, "offsets", "glenn"))
        if not rows:
            print(f"{gauge}: no usable calibrations in {args.table}; not corrected")
            continue
        model = DriftModel(gauge, rows, args.method, args.max_extrapolation_days, not args.no_thermal,
                           getattr(args, "offsets", "glenn"), getattr(args, "event_scatter", EVENT_SCATTER_PA),
                           getattr(args, "tail_events", TAIL_EVENTS),
                           getattr(args, "max_drift_sigma", MAX_DRIFT_SIGMA_PA))
        models[gauge] = model
        reports.append(model.report(rows))
        if verbose:
            print(f"{gauge}: {len(rows)} calibrations ({rows[0]['t0_utc'][:10]} .. {rows[-1]['t0_utc'][:10]}), "
                  f"model {model.kind}" + (f", tau {model.tau:.0f} d" if model.tau else "")
                  + f", rms residual {model.rms:.1f} Pa")
    if not models:
        raise ValueError(f"No calibrations to build a drift curve from in {args.table}")
    report = {**settings(args), "dwt_area_alpha_per_c": DWT_AREA_ALPHA_PER_C,
              "mask": {"before_s": MASK_BEFORE_S, "after_s": args.settle_seconds, "threshold_pa": THRESHOLD_PA},
              "gauges": reports, "paro1_minus_paro2_drift_at_events_pa": paro_difference(models)}
    return models, report


def correct_day(client: Client, models: dict[str, DriftModel], day: date, out: Path, settle: float) -> dict:
    """Write drift-corrected LDO for one day; returns per-gauge statistics."""
    from scpr_mseed import sds_path, write_channel    # same SDS layout and read-back check
    lo = UTCDateTime(day.isoformat())
    ft, flag = load(client, "40", "LEF", lo - 86400, lo + 2 * 86400)
    entry: dict = {"date_utc": str(day)}
    for gauge, model in models.items():
        loc = GAUGE_LOC[gauge]
        # Correct exactly the traces of this day's raw file. A time-window query would also
        # return the previous day's last sample when its grid time falls just after midnight,
        # and readers would then join that 1-sample trace onto this day's trace (1-3 s shifts).
        raw_file = sds_path(Path(client.sds_root), loc, "LDO", day)
        if not raw_file.exists():
            continue
        st = read(str(raw_file), format="MSEED")
        pt, p = load(client, loc, "LDO", lo - 86400, lo + 2 * 86400)
        masks = mask_intervals(ft, flag, pt, p, settle)
        traces, kept, n_in, n_mask, n_late = [], [], 0, 0, 0
        st.sort(keys=["starttime"])
        for tr in st:
            t = tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate
            drift = model.delta(t)
            masked = in_intervals(t, masks)
            keep = ~masked & np.isfinite(drift)
            n_in += len(t)
            n_mask += int(masked.sum())
            n_late += int((~masked & ~np.isfinite(drift)).sum())
            edges = np.flatnonzero(np.diff(np.r_[0, keep.astype(int), 0]))
            for a, b in zip(edges[::2], edges[1::2]):
                piece = tr.data[a:b].astype(float) - drift[a:b]
                new = tr.copy()
                new.data = piece
                new.stats.starttime = tr.stats.starttime + a / tr.stats.sampling_rate
                new.stats.mseed = {"dataquality": "Q"}
                traces.append(new)
                kept.append(piece)
        last_event = model.t[-1]
        stats = {"samples_in": n_in, "masked_calibration": n_mask, "past_extrapolation_limit": n_late,
                 "written": int(sum(len(k) for k in kept)),
                 "last_calibration_utc": iso_s(last_event),
                 # After the last calibration the curve is extrapolated: re-run when the next one arrives.
                 "provisional": bool(lo.timestamp + 86400 > last_event),
                 "valid_until_utc": iso_s(model.t_last_valid),
                 "drift_sigma_pa_at_day_start_end": [float(x) for x in
                                                     model.sigma_drift(np.array([lo.timestamp, lo.timestamp + 86399.0]))],
                 "drift_at_day_start_end_pa": [float(x) if np.isfinite(x) else None
                                               for x in model.delta(np.array([lo.timestamp, lo.timestamp + 86399.0]))]}
        if traces:
            stats.update(write_channel(sds_path(out.resolve(), loc, "LDO", day), traces, np.concatenate(kept)))
        entry[gauge] = stats
    return entry


def write_day_manifest(args: argparse.Namespace, entry: dict) -> Path:
    args.manifest_dir.mkdir(parents=True, exist_ok=True)
    path = args.manifest_dir / f"scpr_corrected_{entry['date_utc']}.json"
    path.write_text(json.dumps({**settings(args), **entry}, indent=2) + "\n", encoding="utf-8")
    return path


def apply(args: argparse.Namespace) -> list[dict]:
    client = Client(str(args.mseed))
    models, report = build_models(args)
    days_out = []
    day = args.start
    while day <= args.end:
        entry = correct_day(client, models, day, args.out, args.settle_seconds)
        if any(g in entry for g in models):
            write_day_manifest(args, entry)
        days_out.append(entry)
        day += timedelta(days=1)

    args.report_dir.mkdir(parents=True, exist_ok=True)
    path = args.report_dir / f"drift_correction_{args.start}_{args.end}.json"
    path.write_text(json.dumps({**report, "days": days_out}, indent=2) + "\n", encoding="utf-8")
    written = sum(e[g]["written"] for e in days_out for g in models if g in e)
    print(f"{written} corrected samples written to {args.out}; report {path}")
    return days_out


def model(args: argparse.Namespace) -> None:
    """Publish the drift curve: calibration/drift_model.json and drift.png."""
    models, report = build_models(args)
    for g in report["gauges"]:
        g["hindcast"] = hindcast(g["gauge"], trusted_rows(args.table, g["gauge"], args.exclude_before,
                                                          getattr(args, "offsets", "glenn")), args)
        h = g["hindcast"]
        if h:
            print(f"{g['gauge']} hindcast: {h['trials']} predictions {h['lead_days'][0]:.0f}-{h['lead_days'][1]:.0f} d "
                  f"ahead, rms error {h['rms_error_pa']:.0f} Pa, mean {h['mean_error_pa']:+.0f} Pa, "
                  f"(error/sigma)^2 {h['mean_squared_error_over_sigma']:.2f}; tail slope "
                  f"{g['tail']['slope_pa_per_year']:.0f} Pa/yr; corrected until {g['valid_until_utc'][:10]}")
    args.model_dir.mkdir(parents=True, exist_ok=True)
    path = args.model_dir / "drift_model.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    plot_drift(models, args.model_dir / "drift.png")
    print(f"Model: {path}; plot: {args.model_dir / 'drift.png'}")


def paro_difference(models: dict[str, DriftModel]) -> list[dict]:
    """PARO1 minus PARO2 drift at each PARO1 event: a gauge-only jump shows up here."""
    if len(models) < 2:
        return []
    m1, m2 = models["PARO1"], models["PARO2"]
    return [{"t_utc": iso_s(t), "difference_pa": float(a - b)}
            for t, a, b in zip(m1.t, m1.delta(m1.t), m2.delta(m1.t))]


def plot_drift(models: dict[str, DriftModel], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.dates import date2num
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    t_all = np.concatenate([m.t for m in models.values()])
    grid = np.linspace(t_all.min(), max(m.t_last_valid for m in models.values()), 2000)
    when = lambda t: date2num([datetime.fromtimestamp(x, timezone.utc) for x in t])
    for (gauge, m), color in zip(models.items(), ("tab:blue", "tab:orange")):
        axes[0].plot(when(grid), m.delta(grid), color=color, label=f"{gauge}: {m.kind}")
        late = grid > m.t[-1]
        axes[0].fill_between(when(grid[late]), (m.delta(grid) - m.sigma_drift(grid))[late],
                             (m.delta(grid) + m.sigma_drift(grid))[late], color=color, alpha=0.2, lw=0)
        axes[0].errorbar(when(m.t), m.y - m.ref, yerr=m.sigma, fmt="o", color=color, ms=4, capsize=2)
    axes[0].set_ylabel("drift (Pa)")
    axes[0].legend(fontsize=8)
    if len(models) == 2:
        axes[1].plot(when(grid), models["PARO1"].delta(grid) - models["PARO2"].delta(grid), color="k")
    axes[1].set_ylabel("PARO1 - PARO2 drift (Pa)")
    axes[1].xaxis_date()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def add_model_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    p.add_argument("--method", choices=("fit", "linear"), default="fit")
    p.add_argument("--max-extrapolation-days", type=float, default=365.0,
                   help="Hard limit after the last calibration (default 365; the sigma limit usually stops first)")
    p.add_argument("--max-drift-sigma", type=float, default=MAX_DRIFT_SIGMA_PA,
                   help="Stop correcting when the drift uncertainty exceeds this (Pa, default 100)")
    p.add_argument("--tail-events", type=int, default=TAIL_EVENTS,
                   help="Calibrations in the line extrapolated after the last one (default 4)")
    p.add_argument("--settle-seconds", type=float, default=1800.0,
                   help="Mask after each calibration ends (default 1800)")
    p.add_argument("--exclude-before", type=date.fromisoformat, default=date.fromisoformat(DEFAULT_EXCLUDE_BEFORE))
    p.add_argument("--no-thermal", action="store_true", help="Skip the DWT thermal adjustment")
    p.add_argument("--offsets", choices=tuple(OFFSET_COLUMNS), default="glenn",
                   help="Offset method: glenn (Sasagawa 1-sigma trim, default) or masked")
    p.add_argument("--event-scatter", type=float, default=EVENT_SCATTER_PA,
                   help="Pa added in quadrature to each event's SE for the fit weights (default 30)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, func, text in (("extract", extract, "Find calibrations and estimate gauge offsets"),
                             ("apply", apply, "Write drift-corrected LDO MiniSEED"),
                             ("model", model, "Write the drift model JSON and plot")):
        p = sub.add_parser(name, help=text)
        p.set_defaults(func=func)
        if name == "model":
            add_model_options(p)
            p.add_argument("--model-dir", type=Path, default=MODEL_DIR)
            continue
        p.add_argument("--start", type=date.fromisoformat, required=True)
        p.add_argument("--end", type=date.fromisoformat, required=True)
        p.add_argument("--mseed", type=Path, default=DEFAULT_MSEED, help="SDS root written by scpr_mseed.py")
        if name == "extract":
            p.add_argument("--table", type=Path, default=DEFAULT_TABLE)
        else:
            add_model_options(p)
            p.add_argument("--out", type=Path, default=DEFAULT_CORRECTED, help="SDS root for corrected LDO")
            p.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
            p.add_argument("--manifest-dir", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    if getattr(args, "start", None) and args.end < args.start:
        parser.error("--end is before --start")
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
