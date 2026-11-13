#!/usr/bin/env python3
"""Find SCPR calibration days on piweb without downloading the whole archive.

piweb serves one SCPRAA301_SCPR log per day. A normal day's log has a fixed
size (86400 type-1 records); a calibration day's log does not, because records
pause during the calibration and the controller's messages are added. So:

  1. scan     HEAD every daily SCPR log (a few hundred bytes each; cached in
              $SCPR_OUTPUT/cache/piweb_scpr_sizes.json, so later scans only see new days);
  2. screen   days whose size is not one of the common "normal" sizes become
              candidates (calibrations, but also outages and restarts);
  3. confirm  download each candidate's SCPR log and look for the controller's
              "Switching SV to DWT" / "Starting Calibration Run" lines; if none,
              download PARO1 and look for pressure below the 15400 kPa threshold.

Confirmed days go to calibration/calibration_days.csv (committed with the code),
cross-checked against the event list in Glenn Sasagawa's CSCPR FAQ.

    python scpr_events.py --start 2018-07-01 --end 2020-12-31
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import date, datetime, timedelta, timezone
import json
from pathlib import Path
import re
import sys
from urllib.request import Request, urlopen

from download_scpr_day import BASE_URL, DirectoryLinks, download
import scpr_paths as paths

CACHE = paths.CACHE / "piweb_scpr_sizes.json"
DAYS_CSV = paths.CALIBRATION / "calibration_days.csv"
RAW_ROOT = paths.RAW
FIRST_DAY = date(2018, 7, 6)                 # deployment (FAQ)
THRESHOLD_KPA = 15400.0
MIN_LOW_SAMPLES = 600                        # 10 min; a calibration holds > 1200 s
NORMAL_SHARE = 0.02                          # a size seen on >= 2% of days counts as normal
FILE_NAME = re.compile(r"SCPRAA301_(PARO1|PARO2|SCPR)_[^/]+_(\d{8})T\d{4}_UTC\.dat")
CAL_TEXT = ("Switching SV to DWT", "Starting Calibration Run")
# 2: also reads untimestamped log lines. Older "no calibration" results for days whose raw
# SCPR log is untimestamped are checked again (from the raw files on disk; no download).
CHECK_VERSION = 2
DAYS_COLUMNS = ["date", "first_marker_utc", "evidence", "faq_event"]

# CSCPR FAQS.docx, "Calibration Events" (complete to 2020-10-06).
FAQ_EVENTS = {
    1: "2018-07-06T08:03:35", 2: "2018-07-06T21:40:06", 3: "2018-07-09T22:13:52",
    4: "2018-07-24T18:11:09", 5: "2018-07-31T17:38:32", 6: "2018-08-07T17:31:10",
    7: "2018-08-14T17:33:30", 8: "2018-09-07T17:35:00", 9: "2018-09-21T17:45:00",
    10: "2018-10-05T16:10:20", 11: "2018-11-01T16:28:30", 12: "2018-12-03T17:38:30",
    13: "2019-01-04T18:10:00", 14: "2019-03-01T17:07:25", 15: "2019-05-01T16:53:30",
    16: "2019-07-01T16:29:15", 17: "2019-08-12T16:12:30", 18: "2019-11-04T17:13:00",
    19: "2020-02-03T17:11:00", 20: "2020-05-04T16:08:00", 21: "2020-07-07T16:02:00",
    22: "2020-10-06T15:56:00",
}


def month_listing(channel: str, year: int, month: int) -> dict[str, list[str]]:
    """{YYYY-MM-DD: [file URLs]} for one channel-month directory (empty if absent)."""
    url = f"{BASE_URL}SCPRAA301_{channel}/{year:04d}/{month:02d}/"
    try:
        with urlopen(url, timeout=60) as response:
            parser = DirectoryLinks()
            parser.feed(response.read().decode("utf-8"))
    except OSError as exc:
        if "404" in str(exc):
            return {}
        raise
    out: dict[str, list[str]] = {}
    for name in sorted(set(parser.links)):
        m = FILE_NAME.fullmatch(name)
        if m and m[1] == channel:
            d = m[2]
            out.setdefault(f"{d[:4]}-{d[4:6]}-{d[6:]}", []).append(url + name)
    return out


def head_size(url: str) -> int:
    with urlopen(Request(url, method="HEAD"), timeout=60) as response:
        return int(response.headers["Content-Length"])


def load_cache() -> dict:
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def save_cache(cache: dict) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    tmp = CACHE.with_suffix(".json.part")
    tmp.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n")
    tmp.replace(CACHE)


def months(start: date, end: date):
    y, m = start.year, start.month
    while (y, m) <= (end.year, end.month):
        yield y, m
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)


def scan(start: date, end: date, cache: dict, workers: int = 8) -> None:
    """Fill cache[day] = {"files": {url: size}} for days not cached yet.

    Days from yesterday on are marked unsettled (piweb may still be writing them)
    and are scanned again next time.
    """
    for day in [d for d, e in cache.items() if e.get("unsettled")]:
        del cache[day]
    settled = datetime.now(timezone.utc).date() - timedelta(days=1)
    todo: list[tuple[str, str]] = []
    fresh: dict[str, dict] = {}
    for y, m in months(start, end):
        for day, urls in month_listing("SCPR", y, m).items():
            if start <= date.fromisoformat(day) <= end and day not in cache:
                fresh[day] = {"files": {}}
                todo += [(day, u) for u in urls]
    with ThreadPoolExecutor(workers) as pool:
        for (day, url), size in zip(todo, pool.map(lambda du: head_size(du[1]), todo)):
            fresh[day]["files"][url] = size
    for day, entry in fresh.items():
        if date.fromisoformat(day) >= settled:
            entry["unsettled"] = True
        cache[day] = entry
    print(f"scan {start}..{end}: {len(fresh)} new day(s), {len(todo)} HEAD request(s)")


def normal_sizes(cache: dict) -> set[int]:
    single = Counter(next(iter(e["files"].values())) for e in cache.values() if len(e["files"]) == 1)
    return {size for size, n in single.items() if n >= NORMAL_SHARE * len(cache)}


def raw_files(day: str, channel: str, urls: list[str]) -> list[Path]:
    folder = RAW_ROOT / day
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for url in urls:
        path = folder / url.rsplit("/", 1)[1]
        try:
            download(url, path)
        except ValueError as exc:            # piweb has some empty (0-byte) files
            print(f"  skipped: {exc}")
            continue
        paths.append(path)
    return paths


def confirm(day: str, entry: dict) -> dict:
    """Look inside a candidate day's logs. Returns the check result stored in the cache."""
    stamp = re.compile(r"^<log-ts:([^>]+)>\s?(.*)$")
    for path in raw_files(day, "SCPR", sorted(entry["files"])):
        with path.open("rb") as handle:
            for raw in handle:
                line = raw.decode("ascii", "replace").replace("\0", "")
                if any(text in line for text in CAL_TEXT):
                    m = stamp.match(line.strip())
                    if m and m[1].startswith(day):
                        return {"calibration": True, "first_marker_utc": m[1][:19], "check_version": CHECK_VERSION,
                                "evidence": "controller log: " + next(t for t in CAL_TEXT if t in m[2])}
                    if not m:
                        # Logged without <log-ts> (untimestamped periods): the file is this day's;
                        # the time comes from the re-timed MiniSEED when the event is extracted.
                        return {"calibration": True, "first_marker_utc": f"{day}T00:00:00",
                                "check_version": CHECK_VERSION,
                                "evidence": "controller log (untimestamped): "
                                            + next(t for t in CAL_TEXT if t in line)}
    # No controller text: fall back to the gauge pressure itself. A calibration holds
    # the gauge below the threshold for > 1200 s; single corrupt records (e.g. a
    # truncated "*0001,5934,5.85") must not count, so require MIN_LOW_SAMPLES.
    y, mo = int(day[:4]), int(day[5:7])
    paro = month_listing("PARO1", y, mo).get(day, [])
    low: list[str] = []
    for path in raw_files(day, "PARO1", paro):
        with path.open("rb") as handle:
            for raw in handle:
                text = raw.decode("ascii", "replace").strip().replace("\0", "")
                m = stamp.match(text)
                payload = m[2] if (m and m[1].startswith(day)) else (text if text.startswith("*") else "")
                if payload.startswith("*"):
                    parts = payload[1:].split(",")
                    try:
                        if len(parts) == 4 and float(parts[1]) < THRESHOLD_KPA:
                            low.append(m[1][:19] if m else f"{day}T00:00:00")
                    except ValueError:
                        continue
    if len(low) >= MIN_LOW_SAMPLES:
        return {"calibration": True, "first_marker_utc": low[0], "check_version": CHECK_VERSION,
                "evidence": f"PARO1 below 15400 kPa for {len(low)} samples (no controller text)"}
    return {"calibration": False, "low_pressure_samples": len(low), "check_version": CHECK_VERSION}


def prune_raw(day: str) -> None:
    """Drop raw files downloaded only to screen a non-calibration day (unless it was converted)."""
    folder = RAW_ROOT / day
    if (paths.MANIFEST / f"scpr_mseed_{day}.json").exists() or not folder.is_dir():
        return
    for path in folder.glob("SCPRAA301_*.dat"):
        path.unlink()
    if not any(folder.iterdir()):
        folder.rmdir()


def untimestamped_raw(day: str) -> bool:
    """The day's raw SCPR log is on disk and its first line has no <log-ts> stamp."""
    for path in sorted((RAW_ROOT / day).glob("SCPRAA301_SCPR_*.dat")):
        with path.open("rb") as handle:
            if not handle.readline().lstrip(b"\0 ").startswith(b"<log-ts:"):
                return True
    return False


def faq_number(day: str) -> str:
    return ";".join(str(n) for n, t in FAQ_EVENTS.items() if t.startswith(day))


def find(start: date, end: date, workers: int = 8) -> list[dict]:
    """Scan, screen and confirm; rewrite calibration_days.csv. Returns all calibration days."""
    cache = load_cache()
    scan(start, end, cache, workers)
    normal = normal_sizes(cache)
    candidates = [d for d, e in sorted(cache.items())
                  if start <= date.fromisoformat(d) <= end and "check" not in e
                  and (len(e["files"]) != 1 or next(iter(e["files"].values())) not in normal)]
    candidates += [d for d, e in sorted(cache.items())
                   if start <= date.fromisoformat(d) <= end and d not in candidates
                   and e.get("check", {}).get("calibration") is False
                   and e["check"].get("check_version", 1) < CHECK_VERSION and untimestamped_raw(d)]
    print(f"normal SCPR log sizes: {sorted(normal)}; {len(candidates)} new candidate day(s) to check")
    for day in candidates:
        result = cache[day]["check"] = confirm(day, cache[day])
        if not result["calibration"]:
            prune_raw(day)
        print(f"  {day}: {'calibration ' + result['first_marker_utc'] + ' (' + result['evidence'] + ')' if result['calibration'] else 'no calibration'}")
        save_cache(cache)
    save_cache(cache)

    days = []
    for d, e in sorted(cache.items()):
        result = e.get("check")
        if result and result["calibration"]:
            days.append({"date": d, "first_marker_utc": result["first_marker_utc"],
                         "evidence": result["evidence"], "faq_event": faq_number(d)})
    DAYS_CSV.parent.mkdir(parents=True, exist_ok=True)
    tmp = DAYS_CSV.with_suffix(".csv.part")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, DAYS_COLUMNS)
        writer.writeheader()
        writer.writerows(days)
    tmp.replace(DAYS_CSV)

    found = {r["date"] for r in days}
    scanned = {d for d in cache if start <= date.fromisoformat(d) <= end}
    missing = [f"{n} ({t})" for n, t in FAQ_EVENTS.items() if t[:10] in scanned and t[:10] not in found]
    print(f"{len(days)} calibration day(s) in {DAYS_CSV}")
    if missing:
        print("WARNING: FAQ events not found: " + ", ".join(missing))
    return days


def read_days() -> list[dict]:
    if not DAYS_CSV.exists():
        return []
    with DAYS_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start", type=date.fromisoformat, default=FIRST_DAY)
    parser.add_argument("--end", type=date.fromisoformat,
                        default=datetime.now(timezone.utc).date() - timedelta(days=1))
    parser.add_argument("--workers", type=int, default=8, help="Parallel HEAD requests")
    args = parser.parse_args()
    find(args.start, args.end, args.workers)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
