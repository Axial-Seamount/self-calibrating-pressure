#!/usr/bin/env python3
"""Export the SCPR pressure record as one time series file (hourly means, CSV).

Reads the raw MiniSEED (output/mseed) and the committed calibration record, and writes
one row per hour with, for each gauge:
  raw pressure, drift-corrected pressure with both calibration averagings (glenn = default,
  masked = OOI practice), the drift uncertainty, and the number of 1-Hz samples averaged.
Calibrations are left out exactly as in the corrected MiniSEED (60 s before LEF = 1 or
pressure < 15 400 kPa, to 30 min after). Corrected values are empty where the drift
uncertainty exceeds the limit (no correction), i.e. after about 2026-08-21.

    python scpr_export.py          # 2018-07-06 .. yesterday -> timeseries/scpr_pressure_hourly.csv,
                                   #   scpr_pressure.png and scpr_uplift.png
    (scpr.py update rebuilds both every day, after the day's correction, and --git-push commits them)
    python scpr_export.py --interval 86400 --out ../output/scpr_pressure_daily.csv
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read
from obspy.clients.filesystem.sds import Client

import scpr_calibrate as sc
import scpr_paths as paths
from scpr import calibration_options
from scpr_mseed import sds_path

FIRST_DAY = date(2018, 7, 6)
REF_DAYS = 30                                   # uplift zero: mean of the first 30 days
GAUGES = ("PARO1", "PARO2")
_STATE: dict = {}


def _models() -> dict:
    if not _STATE:
        _STATE["glenn"] = sc.build_models(calibration_options(offsets="glenn"), verbose=False)[0]
        _STATE["masked"] = sc.build_models(calibration_options(offsets="masked"), verbose=False)[0]
        _STATE["client"] = Client(str(paths.MSEED))
    return _STATE


def godin(x: np.ndarray, interval: int) -> np.ndarray:
    """Godin tidal filter (running means of 24, 24 and 25 h) on an evenly spaced series;
    NaN where a window is less than 90 % covered. For daily or longer intervals: unchanged."""
    if interval >= 86400:
        return x.copy()
    steps = [round(h * 3600 / interval) for h in (24, 24, 25)]
    y = x.astype(float)
    for w in steps:
        ok = np.isfinite(y)
        k = np.ones(w)
        n = np.convolve(ok.astype(float), k, "same")
        with np.errstate(invalid="ignore", divide="ignore"):
            y = np.where(n >= 0.9 * w, np.convolve(np.where(ok, y, 0.0), k, "same") / n, np.nan)
    return y


def _epoch(d: date) -> float:
    return datetime(d.year, d.month, d.day, tzinfo=timezone.utc).timestamp()


def day_rows(job: tuple[date, int]) -> list[tuple]:
    """Bin one day: [(bin start epoch, {gauge: (raw, corr_glenn, corr_masked, sigma, n)})]."""
    day, interval = job
    st = _models()
    client, settle = st["client"], 1800.0
    lo = UTCDateTime(day.isoformat())
    day0 = lo.timestamp
    nbin = 86400 // interval
    try:
        ft, flag = sc.load(client, "40", "LEF", lo - 86400, lo + 2 * 86400)
    except Exception:
        ft, flag = np.array([]), np.array([])
    out = {}
    for gauge in GAUGES:
        loc = sc.GAUGE_LOC[gauge]
        f = sds_path(paths.MSEED, loc, "LDO", day)
        if not f.exists():
            continue
        tr_list = read(str(f), format="MSEED")
        t = np.concatenate([tr.stats.starttime.timestamp + np.arange(tr.stats.npts) / tr.stats.sampling_rate
                            for tr in tr_list])
        p = np.concatenate([tr.data.astype(float) for tr in tr_list])
        pt, pp = sc.load(client, loc, "LDO", lo - 86400, lo + 2 * 86400)
        keep = ~sc.in_intervals(t, sc.mask_intervals(ft, flag, pt, pp, settle)) & (t >= day0) & (t < day0 + 86400)
        t, p = t[keep], p[keep]
        b = ((t - day0) // interval).astype(int)
        n = np.bincount(b, minlength=nbin).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            raw = np.bincount(b, p, nbin) / n
            cols = [raw]
            for method in ("glenn", "masked"):
                corr = p - st[method][gauge].delta(t)               # NaN past the correction limit
                ok = np.isfinite(corr)
                nc = np.bincount(b[ok], minlength=nbin).astype(float)
                c = np.bincount(b[ok], corr[ok], nbin) / nc
                c[nc < n] = np.nan                                  # partly past the limit: leave empty
                cols.append(c)
        mid = day0 + (np.arange(nbin) + 0.5) * interval
        cols.append(st["glenn"][gauge].sigma_drift(mid))
        cols.append(n)
        out[gauge] = np.column_stack(cols)
    return [(day0 + i * interval, {g: out[g][i] for g in out if out[g][i, 4] > 0}) for i in range(nbin)]


DEFAULT_OUT = paths.TIMESERIES / "scpr_pressure_hourly.csv"
DEFAULT_PNG = paths.TIMESERIES / "scpr_pressure.png"
UPLIFT_PNG = paths.TIMESERIES / "scpr_uplift.png"
# Seafloor pressure change -> vertical displacement: dP = -rho * g * dz, with the in-situ
# seawater density at the SCPR (~1530 m, ~2.3 degC, S ~34.6: ~1035 kg/m3) and g at 46 N.
RHO_SEAFLOOR, G_LOCAL = 1035.0, 9.808


def export(start: date = FIRST_DAY, end: date | None = None, interval: int = 3600,
           out: Path = DEFAULT_OUT, workers: int = 6) -> Path:
    """Write the whole time series file (rebuilt from scratch: a new calibration changes it all)."""
    from scpr import yesterday
    args = argparse.Namespace(start=start, end=end or yesterday(), interval=interval,
                              out=out, workers=workers)
    if 86400 % args.interval:
        raise ValueError("interval must divide 86400")
    days = [args.start + timedelta(days=i) for i in range((args.end - args.start).days + 1)]
    with ProcessPoolExecutor(args.workers) as ex:
        rows = [r for chunk in ex.map(day_rows, [(d, args.interval) for d in days], chunksize=8) for r in chunk]
    models = _models()
    version = sc.table_version(sc.DEFAULT_TABLE)
    last = {g: sc.iso_s(models["glenn"][g].t[-1]) for g in GAUGES}
    valid = {g: sc.iso_s(models["glenn"][g].t_last_valid) for g in GAUGES}
    # Uplift: tides removed from the corrected pressure (Godin 24-24-25 h filter on the hourly
    # series), relative to its mean over the first REF_DAYS days of the record (a fixed zero,
    # whatever --start is), divided by rho*g.
    rho_g = RHO_SEAFLOOR * G_LOCAL
    detided, ref = {}, {}
    for g in GAUGES:
        corr = np.array([vals[g][1] if g in vals else np.nan for _, vals in rows])
        detided[g] = godin(corr, args.interval)
    t0s = np.array([r[0] for r in rows])
    in_ref = t0s < _epoch(FIRST_DAY + timedelta(days=REF_DAYS))
    for g in GAUGES:
        v = detided[g][in_ref]
        ref[g] = float(np.nanmean(v)) if args.start == FIRST_DAY and np.isfinite(v).any() else np.nan
    if args.start != FIRST_DAY:
        print("uplift columns left empty: the reference period needs --start", FIRST_DAY)
    names = ("raw", "corr", "corr_masked", "uplift", "drift_sigma", "n")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".part")
    written = 0
    with tmp.open("w") as fh:
        fh.write(f"# SCPR OO.AXCC1 seafloor pressure, {args.interval}-s means, {args.start} .. {args.end} (UTC)\n"
                 "# Gauges: PARO1 = loc 40 (Paroscientific SN 136576), PARO2 = loc 41 (SN 132123). Pressure in Pa (absolute).\n"
                 "# time_utc: start of the averaging interval. Calibrations removed (60 s before to 30 min after).\n"
                 "# *_raw: gauge pressure as recorded. *_corr: raw minus drift, calibration averaging 'glenn'\n"
                 "#   (Sasagawa et al. 2016; default). *_corr_masked: the same with the 'masked' averaging (spin-up\n"
                 "#   skipped, re-spins removed). Both are relative to each gauge's first calibration (2018-07-09).\n"
                 f"# *_uplift_m: estimated seafloor uplift (m) = -(detided corr - reference) / (rho g), rho = {RHO_SEAFLOOR:.0f} kg/m3,\n"
                 "#   tides removed with a Godin 24-24-25 h filter (empty within ~1.5 days of a data gap),\n"
                 f"#   g = {G_LOCAL} m/s2 (1 kPa = {1000 / rho_g * 100:.2f} cm); reference = mean corr of {FIRST_DAY} .. "
                 f"{FIRST_DAY + timedelta(days=REF_DAYS - 1)}\n"
                 f"#   (PARO1 {ref['PARO1']:.1f} Pa, PARO2 {ref['PARO2']:.1f} Pa). No reference station: ocean pressure changes remain.\n"
                 "# *_drift_sigma: 1-sigma uncertainty of the drift removed. *_n: 1-Hz samples averaged.\n"
                 f"# Empty corr: drift uncertainty > 100 Pa (corrected until {valid['PARO1'][:10]}); after the last calibration\n"
                 f"#   ({last['PARO1'][:10]}) the drift is extrapolated (provisional).\n"
                 f"# Calibration table {version}; generated {datetime.now(timezone.utc):%Y-%m-%dT%H:%M}Z by scpr_export.py\n")
        unit = {"n": "", "uplift": "_m"}
        fh.write("time_utc," + ",".join(f"{g.lower()}_{n}{unit.get(n, '_pa')}" for g in GAUGES for n in names) + "\n")
        for i, (t0, vals) in enumerate(rows):
            if not vals:
                continue
            cells = [datetime.fromtimestamp(t0, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")]
            for g in GAUGES:
                v = vals.get(g)
                if v is None:
                    cells += [""] * 6
                    continue
                up = -(detided[g][i] - ref[g]) / rho_g
                f1 = lambda x: "" if not np.isfinite(x) else f"{x:.1f}"
                cells += [f1(v[0]), f1(v[1]), f1(v[2]), "" if not np.isfinite(up) else f"{up:.4f}", f1(v[3]),
                          str(int(v[4]))]
            fh.write(",".join(cells) + "\n")
            written += 1
    tmp.replace(args.out)
    print(f"Time series: {written} rows -> {args.out}")
    if args.out == DEFAULT_OUT:
        plot(args.out, DEFAULT_PNG)
        plot_uplift(args.out, UPLIFT_PNG)
    return args.out


# Categorical slots 1-2 of the reference palette (validated light-mode pair); text in ink tokens.
COLORS = {"PARO1": "#2a78d6", "PARO2": "#eb6834"}
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def label_extrapolated(ax, last_cal, valid, y_frac: float = 0.03, va: str = "bottom") -> None:
    """Say what the shaded span after the last calibration is, inside it."""
    import matplotlib.transforms as mtransforms
    tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(last_cal, y_frac, f" after the last\n calibration\n ({str(last_cal)[:10]}):\n drift extrapolated;\n"
                              f" not corrected\n after {str(valid)[:10]}",
            transform=tr, color=INK2, fontsize=7, va=va, ha="left", linespacing=1.15)


def read_csv(path: Path) -> tuple[np.ndarray, dict, list[str]]:
    lines = path.read_text().splitlines()
    head = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    cols = body[0].split(",")
    t = np.array([np.datetime64(r.split(",", 1)[0].rstrip("Z")) for r in body[1:]])
    vals = np.array([[float(x) if x else np.nan for x in r.split(",")[1:]] for r in body[1:]])
    return t, {c: vals[:, i] for i, c in enumerate(cols[1:])}, head


def daily(t: np.ndarray, v: np.ndarray, min_hours: int = 20) -> tuple[np.ndarray, np.ndarray]:
    d = t.astype("datetime64[D]")
    days, inv = np.unique(d, return_inverse=True)
    ok = np.isfinite(v)
    n = np.bincount(inv[ok], minlength=len(days))
    s = np.bincount(inv[ok], v[ok], minlength=len(days))
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.where(n >= min_hours, s / n, np.nan)
    return days.astype("datetime64[h]") + 12, m


def plot(csv: Path, png: Path) -> Path:
    """Two panels from the time series file: corrected pressure change, and PARO1 - PARO2."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    t, c, _ = read_csv(csv)
    models = _models()["glenn"]
    cal_days = np.array(sorted({sc.iso_s(x)[:10] for g in GAUGES for x in models[g].t}), dtype="datetime64[D]")
    valid = np.datetime64(sc.iso_s(models["PARO1"].t_last_valid)[:19])
    last_cal = np.datetime64(sc.iso_s(models["PARO1"].t[-1])[:19])
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True, height_ratios=(2.2, 1),
                                 facecolor=SURFACE, constrained_layout=True)
    for ax in (a1, a2):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.6)
        ax.axvspan(last_cal, t[-1], color=GRID, alpha=0.45, lw=0)
    for g in GAUGES:
        key = g.lower()
        x, y = daily(t, c[f"{key}_corr_pa"])
        ref = y[np.isfinite(y)][:30].mean()
        y = (y - ref) / 1000.0
        _, sig = daily(t, c[f"{key}_drift_sigma_pa"])
        a1.fill_between(x, y - sig / 1000, y + sig / 1000, color=COLORS[g], alpha=0.12, lw=0)
        a1.plot(x, y, color=COLORS[g], lw=1.2, label=f"{g} (loc {sc.GAUGE_LOC[g]})")
    lo = a1.get_ylim()[0]
    a1.plot(cal_days.astype("datetime64[h]"), np.full(len(cal_days), lo), "|", color=INK2, ms=7, mew=1,
            label="calibration")
    a1.set_ylabel("Pressure change (kPa)")
    a1.set_title("SCPR seafloor pressure, drift-corrected (daily means; shading = drift uncertainty, 1\u03c3)",
                 loc="left", color=INK, fontsize=10)
    label_extrapolated(a1, last_cal, valid, y_frac=0.97, va="top")
    a1.legend(loc="upper right", bbox_to_anchor=(0.84, 1.0), frameon=False, fontsize=8, ncol=1)
    for kind, colr, lw, lab in (("raw", INK2, 0.9, "raw"), ("corr", COLORS["PARO1"], 1.2, "drift-corrected")):
        x, d1 = daily(t, c[f"paro1_{kind}_pa"]); _, d2 = daily(t, c[f"paro2_{kind}_pa"])
        dd = (d1 - d2)
        dd = (dd - np.nanmean(dd[np.isfinite(dd)][:30])) / 1000.0
        a2.plot(x, dd, color=colr, lw=lw, label=lab)
    a2.set_ylabel("PARO1 \u2212 PARO2 (kPa)")
    a2.set_title("Difference between the two gauges (relative to the first month)", loc="left", color=INK,
                 fontsize=10)
    a2.legend(loc="upper left", frameon=False, fontsize=8, ncol=2)
    a2.xaxis.set_major_locator(mdates.YearLocator())
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    a2.set_xlabel(f"OO.AXCC1 SCPR, Axial Seamount. Updated {str(t[-1])[:10]} from timeseries/{csv.name}; "
                  f"corrected until {str(valid)[:10]}.", fontsize=7, color=INK2, loc="right")
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Figure -> {png}")
    return png


def running_mean(y: np.ndarray, width: int) -> np.ndarray:
    ok = np.isfinite(y)
    k = np.ones(width)
    n = np.convolve(ok.astype(float), k, "same")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n >= width * 0.6, np.convolve(np.where(ok, y, 0.0), k, "same") / n, np.nan)


def plot_uplift(csv: Path, png: Path) -> Path:
    """Estimated seafloor uplift since the first calibration, and its rate."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    t, c, _ = read_csv(csv)
    models = _models()["glenn"]
    last_cal = np.datetime64(sc.iso_s(models["PARO1"].t[-1])[:19])
    rho_g = RHO_SEAFLOOR * G_LOCAL
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True, height_ratios=(2.2, 1),
                                 facecolor=SURFACE, constrained_layout=True)
    ups = []
    for g in GAUGES:
        x, up = daily(t, c[f"{g.lower()}_uplift_m"])             # m, from the file's uplift column
        ups.append(up)
        a1.plot(x, up, color=COLORS[g], lw=0.7, alpha=0.55, label=f"{g} (daily mean)")
    with np.errstate(invalid="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", RuntimeWarning)       # days with no gauge
        mean_up = np.nanmean(np.vstack(ups), axis=0)
    smooth = running_mean(mean_up, 30)
    a1.plot(x, smooth, color=INK, lw=1.6, label="both gauges, 30-day mean")
    fin = np.flatnonzero(np.isfinite(smooth))
    a1.annotate(f"{smooth[fin[-1]] * 100:+.0f} cm", (x[fin[-1]], smooth[fin[-1]]), xytext=(5, 0),
                textcoords="offset points", va="center", color=INK, fontsize=8.5)
    valid = np.datetime64(sc.iso_s(models["PARO1"].t_last_valid)[:19])
    label_extrapolated(a1, last_cal, valid)
    a1.set_ylabel("Uplift since July 2018 (m)")
    a1.set_title("Estimated seafloor uplift at the central caldera, from drift-corrected SCPR pressure",
                 loc="left", color=INK, fontsize=10)
    a1.legend(loc="upper left", frameon=False, fontsize=8)
    # rate: slope of the 30-day mean over a centred 180-day window, cm/yr
    rate = np.full(len(smooth), np.nan)
    half, days = 90, (x - x[0]).astype("timedelta64[h]").astype(float) / 24.0
    for i in range(half, len(smooth) - half):
        w = slice(i - half, i + half + 1)
        ok = np.isfinite(smooth[w])
        if ok.sum() > 0.8 * (2 * half + 1):
            rate[i] = np.polyfit(days[w][ok], smooth[w][ok], 1)[0] * 365.25 * 100
    a2.plot(x, rate, color=INK, lw=1.3)
    a2.axhline(0, color=INK2, lw=0.6)
    a2.set_ylabel("Uplift rate (cm/yr)")
    a2.set_title("Uplift rate (180-day window)", loc="left", color=INK, fontsize=10)
    for ax in (a1, a2):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.6)
        ax.axvspan(last_cal, t[-1], color=GRID, alpha=0.45, lw=0)
    a2.xaxis.set_major_locator(mdates.YearLocator())
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    a2.set_xlabel(f"Uplift = -\u0394P / (\u03c1g), \u03c1 = {RHO_SEAFLOOR:.0f} kg/m\u00b3, g = {G_LOCAL} m/s\u00b2 "
                  f"(1 kPa \u2248 {1000 / rho_g * 100:.1f} cm). No reference station: ocean pressure changes "
                  f"(~cm) remain. Shaded: after the last calibration. Updated {str(t[-1])[:10]}.",
                  fontsize=7, color=INK2, loc="right")
    png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f"Figure -> {png}")
    return png


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--start", type=date.fromisoformat, default=FIRST_DAY)
    ap.add_argument("--end", type=date.fromisoformat, default=None, help="Last day (default: yesterday, UTC)")
    ap.add_argument("--interval", type=int, default=3600, help="Averaging interval in s; must divide 86400")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--plot-only", action="store_true", help="Redraw the figure from the existing file")
    a = ap.parse_args()
    if a.plot_only:
        plot(a.out, DEFAULT_PNG)
        plot_uplift(a.out, UPLIFT_PNG)
        return
    export(a.start, a.end, a.interval, a.out, a.workers)


if __name__ == "__main__":
    main()
