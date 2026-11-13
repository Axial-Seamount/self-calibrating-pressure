#!/usr/bin/env python3
"""SCPR at Axial (OO.AXCC1): one command from piweb to drift-corrected MiniSEED.

Users
  day DATE [--end DATE]
      Download the day(s) from piweb, convert to raw MiniSEED (data quality D,
      output/mseed/) and write drift-corrected LDO (quality Q, output/mseed_corrected/)
      using the committed calibration/calibrations.csv. No archive needed.
      Days after the last calibration are marked provisional and carry a drift
      uncertainty; once that exceeds 100 Pa (~200 days) they are converted but
      not corrected.

Maintainers (the calibration table is shared through git)
  build [--start DATE] [--end DATE]
      Find calibration days on piweb (scpr_events.py), convert only those days,
      extract every calibration into calibration/calibrations.csv and rewrite
      calibration/drift_model.json + drift.png. Commit the calibration/ folder.
  update [--no-yesterday]
      Daily (cron): scan new days, add any new calibration to the table, refit the
      drift model, re-correct EVERY local corrected day whose correction used an older
      table or drift logic (a new calibration changes the whole curve), and convert +
      correct yesterday.

Run with the conda env ooi_env (numpy, ObsPy, matplotlib):
    python scpr.py day 2020-05-04
    python scpr.py day 2020-05-01 --end 2020-05-31
    python scpr.py build
    python scpr.py update
"""

from __future__ import annotations

import argparse
from datetime import date, datetime, timedelta, timezone
import json
from pathlib import Path
import sys

from obspy.clients.filesystem.sds import Client

import scpr_calibrate as sc
import scpr_events as ev
from scpr_mseed import CONVERTER_VERSION, DEFAULT_MANIFEST, DEFAULT_MSEED, convert_day, sds_path

ROOT = Path(__file__).resolve().parent
RAW_ROOT = ev.RAW_ROOT


def yesterday() -> date:
    return datetime.now(timezone.utc).date() - timedelta(days=1)


def calibration_options(**overrides) -> argparse.Namespace:
    """The scpr_calibrate settings, defaults unless overridden."""
    opts = dict(table=sc.DEFAULT_TABLE, method="fit", max_extrapolation_days=365.0, settle_seconds=1800.0,
                exclude_before=date.fromisoformat(sc.DEFAULT_EXCLUDE_BEFORE), no_thermal=False,
                mseed=DEFAULT_MSEED, out=sc.DEFAULT_CORRECTED, manifest_dir=DEFAULT_MANIFEST,
                report_dir=sc.DEFAULT_REPORT, model_dir=sc.MODEL_DIR)
    opts.update(overrides)
    return argparse.Namespace(**opts)


def ensure_converted(day: date, force: bool = False) -> bool:
    """Raw MiniSEED for the day exists (downloading/converting if needed). False if piweb has nothing."""
    manifest = DEFAULT_MANIFEST / f"scpr_mseed_{day}.json"
    if (not force and manifest.exists() and sds_path(DEFAULT_MSEED, "40", "LDO", day).exists()
            and json.loads(manifest.read_text()).get("converter_version") == CONVERTER_VERSION):
        return True
    try:
        convert_day(day, RAW_ROOT / str(day), download=True, verbose=False)
    except ValueError as exc:
        if "No raw SCPR files" in str(exc):
            return False
        raise
    return True


def late_calibration_before(day: date) -> bool:
    """A calibration late on the previous day can reach past midnight (mask needs that day)."""
    prev = str(day - timedelta(days=1))
    return any(r["date"] == prev and r["first_marker_utc"][11:13] >= "21" for r in ev.read_days())


def describe(entry: dict) -> str:
    parts = []
    for gauge in sc.GAUGE_LOC:
        s = entry.get(gauge)
        if s is None:
            parts.append(f"{gauge} no data")
        elif not s["written"]:
            parts.append(f"{gauge} NOT corrected (after {s['valid_until_utc'][:10]}: drift uncertainty "
                         f"> {sc.MAX_DRIFT_SIGMA_PA:.0f} Pa since the last calibration {s['last_calibration_utc'][:10]})")
        else:
            note = " PROVISIONAL" if s["provisional"] else ""
            sigma = max(s["drift_sigma_pa_at_day_start_end"])
            parts.append(f"{gauge} {s['written']} samples, {s['masked_calibration']} masked, "
                         f"drift +/-{sigma:.0f} Pa{note}")
    return "; ".join(parts)


# ── users ───────────────────────────────────────────────────────────────────
def _convert_task(job: tuple[date, bool]) -> tuple[date, str | None, bool]:
    """Worker: (day, error or None, has data)."""
    day, reconvert = job
    try:
        return day, None, ensure_converted(day, reconvert)
    except (OSError, ValueError) as exc:
        return day, str(exc), False


_WORKER_STATE: dict = {}


def _correct_task(job: tuple[date, str]) -> tuple[date, str | None, dict | None]:
    """Worker: correct one already-converted day; (day, error or None, manifest entry)."""
    day, method = job
    try:
        if _WORKER_STATE.get("method") != method:      # models and client once per worker process
            opts = calibration_options(method=method)
            _WORKER_STATE.update(method=method, opts=opts, models=sc.build_models(opts, verbose=False)[0],
                                 client=Client(str(DEFAULT_MSEED)))
        st = _WORKER_STATE
        entry = sc.correct_day(st["client"], st["models"], day, st["opts"].out, st["opts"].settle_seconds)
        if any(g in entry for g in st["models"]):
            sc.write_day_manifest(st["opts"], entry)
        return day, None, entry
    except (OSError, ValueError) as exc:
        return day, str(exc), None


def cmd_day(args: argparse.Namespace) -> None:
    end = args.end or args.date
    if end < args.date:
        raise ValueError("--end is before DATE")
    if not sc.DEFAULT_TABLE.exists():
        raise ValueError(f"{sc.DEFAULT_TABLE} is missing; pull the repo or run 'python scpr.py build'")
    opts = calibration_options(method=args.method)
    models, _ = sc.build_models(opts, verbose=False)
    print(f"Calibration table {sc.table_version(opts.table)}: "
          + ", ".join(f"{g} {len(m.t)} events to {sc.iso_s(m.t[-1])[:10]}" for g, m in models.items()))
    DEFAULT_MSEED.mkdir(parents=True, exist_ok=True)     # fresh machine: the SDS root must exist
    days = [args.date + timedelta(days=i) for i in range((end - args.date).days + 1)]
    # A calibration late on the day before the range can reach past midnight (mask needs it).
    extra = [args.date - timedelta(days=1)] if late_calibration_before(args.date) else []
    provisional, failed, nodata = False, [], []
    workers = max(1, getattr(args, "workers", 1))

    # Phase 1: download + convert, each day exactly once (also the calibration-adjacent days
    # inside the range, since every day of the range is converted here).
    from concurrent.futures import ProcessPoolExecutor
    pool = ProcessPoolExecutor(workers) if workers > 1 else None
    mapper = pool.map if pool else map
    converted = []
    try:
        for day, error, has_data in mapper(_convert_task, [(d, args.reconvert and d not in extra)
                                                           for d in extra + days]):
            if day in extra:
                continue
            if error:
                print(f"{day}: ERROR {error}", flush=True)
                failed.append(day)
            elif not has_data:
                print(f"{day}: no SCPR data on piweb", flush=True)
                nodata.append(day)
            else:
                converted.append(day)
        # A day whose previous day had a late calibration needs that day converted too.
        for day in [d for d in converted if late_calibration_before(d)]:
            prev = day - timedelta(days=1)
            if prev not in days:
                _convert_task((prev, False))

        # Phase 2: drift-correct (reads the converted days, writes only its own day).
        for day, error, entry in mapper(_correct_task, [(d, args.method) for d in converted]):
            if error:
                print(f"{day}: ERROR {error}", flush=True)
                failed.append(day)
                continue
            provisional |= any(entry.get(g, {}).get("provisional") and entry[g]["written"] for g in models)
            print(f"{day}: {describe(entry)}", flush=True)
    finally:
        if pool:
            pool.shutdown()
    print(f"{len(converted)} day(s) converted and corrected, {len(nodata)} without data, {len(failed)} failed")
    print(f"Raw MiniSEED (D): {DEFAULT_MSEED}\nCorrected LDO (Q): {opts.out}")
    if provisional:
        print("PROVISIONAL: these days are after the last calibration; the drift there is extrapolated.\n"
              "  Run this again after the next calibration is added to the table (git pull).")
    if failed:
        raise ValueError(f"{len(failed)} day(s) failed: {', '.join(map(str, failed))}")


# ── maintainers ─────────────────────────────────────────────────────────────
def build(start: date, end: date) -> tuple[str, str]:
    """Discover, convert and extract calibrations in [start, end]; refresh the model."""
    before = sc.table_version(sc.DEFAULT_TABLE)
    days = [r for r in ev.find(start, end) if start <= date.fromisoformat(r["date"]) <= end]
    print(f"Converting and extracting {len(days)} calibration day(s)")
    for r in days:
        day = date.fromisoformat(r["date"])
        if not ensure_converted(day):
            print(f"{day}: no SCPR data on piweb")
            continue
        if r["first_marker_utc"][11:13] >= "21":       # event may run past midnight
            ensure_converted(day + timedelta(days=1))
        sc.extract(argparse.Namespace(mseed=DEFAULT_MSEED, table=sc.DEFAULT_TABLE, start=day, end=day))
    if sc.DEFAULT_TABLE.exists():
        sc.model(calibration_options())
    after = sc.table_version(sc.DEFAULT_TABLE)
    print(f"Calibration table: {before} -> {after}" if after != before else "Calibration table unchanged")
    return before, after


def stale_days(version: str) -> list[date]:
    out = []
    for path in sorted(DEFAULT_MANIFEST.glob("scpr_corrected_*.json")):
        info = json.loads(path.read_text())
        # Re-correct when the table or the model weighting behind a day's correction changed.
        if info.get("table_version") != version or info.get("drift_model_version") != sc.DRIFT_MODEL_VERSION:
            out.append(date.fromisoformat(info["date_utc"]))
    return out


def cmd_build(args: argparse.Namespace) -> None:
    build(args.start, args.end)


def publish_calibration(message: str) -> bool:
    """Commit and push calibration/ if it changed. Uses the machine's own git identity."""
    import subprocess
    run = lambda *cmd: subprocess.run(["git", "-C", str(ROOT), *cmd], check=True, capture_output=True, text=True)
    target = str(sc.MODEL_DIR)
    if not run("status", "--porcelain", "--", target).stdout.strip():
        print("calibration/ unchanged; nothing to push")
        return False
    run("add", "--", target)
    run("commit", "-m", message, "--", target)
    run("push")
    print(f"Pushed calibration/: {message}")
    return True


def cmd_update(args: argparse.Namespace) -> None:
    cache = ev.load_cache()
    settled = [d for d, e in cache.items() if not e.get("unsettled")]
    start = date.fromisoformat(max(settled)) - timedelta(days=3) if settled else ev.FIRST_DAY
    before, version = build(start, yesterday())
    redo = stale_days(version)
    if not args.no_yesterday and yesterday() not in redo:
        redo.append(yesterday())
    failed = []
    if redo:
        print(f"Correcting {len(redo)} day(s) with table {version}")
        for day in redo:
            try:
                cmd_day(argparse.Namespace(date=day, end=None, method="fit", reconvert=False, workers=1))
            except (OSError, ValueError) as exc:     # e.g. yesterday not timestamped; keep going
                print(f"{day}: {exc}")
                failed.append(day)
    if args.git_push:
        rows = [r for r in sc.read_table(sc.DEFAULT_TABLE) if r["status"] == "ok"]
        last = max((r["t0_utc"][:10] for r in rows), default="none")
        try:
            publish_calibration(f"Update SCPR calibration table ({len(rows) // 2} events, last {last})")
        except Exception as exc:                      # report, but keep the processing result
            stderr = getattr(exc, "stderr", "") or ""
            raise ValueError(f"git push of calibration/ failed: {exc} {stderr.strip()}") from exc
    elif before != version:
        print("calibration/ changed: commit and push it (or run update --git-push) so users get it")
    if failed:
        # Non-zero exit so cron/monitoring notices; everything else above still ran.
        raise ValueError(f"{len(failed)} day(s) could not be processed: {', '.join(map(str, failed))}")
    print("Done.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("day", help="Download, convert and drift-correct one day or a range")
    p.add_argument("date", type=date.fromisoformat, help="UTC day YYYY-MM-DD")
    p.add_argument("--end", type=date.fromisoformat, help="Last day of a range (inclusive)")
    p.add_argument("--method", choices=("fit", "linear"), default="fit",
                   help="Drift curve: fit (default) or linear interpolation between calibrations")
    p.add_argument("--reconvert", action="store_true", help="Convert again even if MiniSEED exists")
    p.add_argument("--workers", type=int, default=1,
                   help="Days processed in parallel (e.g. 4 for a backfill; each downloads ~22 MB/day)")
    p.set_defaults(func=cmd_day)
    p = sub.add_parser("build", help="Maintainers: build the calibration table from piweb")
    p.add_argument("--start", type=date.fromisoformat, default=ev.FIRST_DAY)
    p.add_argument("--end", type=date.fromisoformat, default=yesterday())
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("update", help="Maintainers: incremental build + re-correct stale days (cron)")
    p.add_argument("--no-yesterday", action="store_true", help="Do not convert and correct yesterday")
    p.add_argument("--git-push", action="store_true",
                   help="Commit and push calibration/ when it changed (machine's git identity)")
    p.set_defaults(func=cmd_update)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
