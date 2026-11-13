#!/usr/bin/env python3
"""Double-check the AXCC1 StationXML before it goes anywhere.

The StationXML is generated from `COSZO Instrument Data Channels.xlsx`, and the sheet
can be wrong. This tests the XML against sources that do not come from the sheet:

  1. the converter (scpr_mseed.CHANNEL_UNITS): the unit each channel is really written in;
  2. the MiniSEED on disk: same channel set, matching sample rates, data inside epochs;
  3. EarthScope's live station record: the station element must be identical (we add
     channels to an existing station), channels must lie near it;
  4. FDSN StationXML conventions: SI-style ASCII unit names, singular "count", no blank or
     padded descriptions, sensitivity 1 for data already in physical units.

Exit status 1 if anything FAILs (WARN does not fail).
    python check_metadata.py                      # output/xml/OO_AXCC1_40.xml and _41.xml
    python check_metadata.py --days 2020-05-04 2024-08-28 --offline
"""

from __future__ import annotations

import argparse
import json
from datetime import date
import math
from pathlib import Path
import sys

from obspy import UTCDateTime, read, read_inventory

import scpr_paths as paths
from scpr_mseed import CHANNEL_UNITS, NETWORK, STATION

XML = [paths.ROOT / "output" / "xml" / f"OO_AXCC1_{loc}.xml" for loc in ("40", "41")]
DEFAULT_DAYS = ["2018-08-14", "2020-05-04", "2024-08-28"]   # early, mid, late (calibration days)
MAX_CHANNEL_OFFSET_M = 200.0
FDSN = "https://service.earthscope.org/fdsnws/station/1/query"

problems: list[str] = []


def report(level: str, message: str) -> None:
    print(f"{level:4s} {message}")
    if level == "FAIL":
        problems.append(message)


def distance_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    dn = (lat2 - lat1) * 111132.9
    de = (lon2 - lon1) * 111319.5 * math.cos(math.radians((lat1 + lat2) / 2))
    return math.hypot(dn, de)


def earthscope_station():
    from urllib.request import urlopen
    url = f"{FDSN}?net={NETWORK}&sta={STATION}&level=station&format=text"
    with urlopen(url, timeout=60) as response:
        lines = [l for l in response.read().decode().splitlines() if l and not l.startswith("#")]
    f = lines[-1].split("|")                       # latest epoch
    return {"lat": float(f[2]), "lon": float(f[3]), "elev": float(f[4]), "site": f[5],
            "start": UTCDateTime(f[6]), "end": UTCDateTime(f[7]) if f[7].strip() else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xml", type=Path, nargs="+", default=XML)
    parser.add_argument("--days", nargs="+", default=DEFAULT_DAYS, help="Days of MiniSEED to compare against")
    parser.add_argument("--offline", action="store_true", help="Skip the EarthScope comparison")
    args = parser.parse_args()

    inv = read_inventory(str(args.xml[0]))
    for f in args.xml[1:]:
        inv += read_inventory(str(f))
    stations = [s for n in inv for s in n]
    channels = {(c.location_code, c.code): (s, c) for s in stations for c in s}
    print(f"{len(channels)} channels in {', '.join(p.name for p in args.xml)}\n")

    # 3. Station element vs EarthScope
    if not args.offline:
        es = earthscope_station()
        for s in stations:
            for name, ours, theirs in (("start", s.start_date, es["start"]), ("end", s.end_date, es["end"]),
                                       ("latitude", s.latitude, es["lat"]), ("longitude", s.longitude, es["lon"]),
                                       ("elevation", s.elevation, es["elev"]), ("site", s.site.name, es["site"])):
                same = abs(ours - theirs) < 1e-6 if isinstance(theirs, float) else ours == theirs
                report("OK" if same else "FAIL", f"station {name}: ours {ours} / EarthScope {theirs}")
            break                                  # both files carry the same station element
        if len({(str(s.start_date), float(s.latitude), float(s.longitude), float(s.elevation), s.site.name)
                for s in stations}) > 1:
            report("FAIL", "the two XML files disagree on the station element")

    # 1 + 4. Units and conventions, per channel
    print()
    for (loc, cha), (sta, c) in sorted(channels.items()):
        sens = c.response.instrument_sensitivity
        tag = f"{loc}.{cha}"
        expected = CHANNEL_UNITS.get(cha)
        if expected is None:
            report("FAIL", f"{tag}: not a channel the converter writes")
            continue
        for side, unit, desc in (("input", sens.input_units, sens.input_units_description),
                                 ("output", sens.output_units, sens.output_units_description)):
            if unit != expected:
                report("FAIL", f"{tag}: {side} unit {unit!r}, converter writes {expected!r}")
            if not unit.isascii():
                report("FAIL", f"{tag}: {side} unit {unit!r} is not ASCII")
            if not desc or desc != desc.strip():
                report("FAIL", f"{tag}: {side} unit description {desc!r} is blank or padded")
        if abs(sens.value - 1.0) > 1e-12:
            report("FAIL", f"{tag}: sensitivity {sens.value}, data are already in {expected} (expected 1)")
        if not c.description or c.description != c.description.strip():
            report("WARN", f"{tag}: channel description {c.description!r}")
        off = distance_m(sta.latitude, sta.longitude, c.latitude, c.longitude)
        if off > MAX_CHANNEL_OFFSET_M:
            report("FAIL", f"{tag}: {off:.0f} m from the station point")
        if c.start_date < sta.start_date or (sta.end_date and c.end_date and c.end_date > sta.end_date):
            report("FAIL", f"{tag}: channel epoch outside the station epoch")
    offsets = {round(distance_m(s.latitude, s.longitude, c.latitude, c.longitude)) for s, c in channels.values()}
    report("OK", f"units, descriptions, sensitivities and epochs checked for {len(channels)} channels; "
                 f"channel offset from the station point: {sorted(offsets)} m")

    # 2. Against MiniSEED on disk
    print()
    for d in args.days:
        day = date.fromisoformat(d)
        written = {}
        for f in sorted((paths.MSEED / str(day.year) / NETWORK / STATION).glob(f"*/{NETWORK}.{STATION}.*.{day.year}.{day.timetuple().tm_yday:03d}")):
            st = read(str(f), headonly=True)
            written[(st[0].stats.location, st[0].stats.channel)] = st
        if not written:
            report("WARN", f"{d}: no MiniSEED on disk to compare with")
            continue
        missing_meta = sorted(set(written) - set(channels))
        no_data = sorted(set(channels) - set(written))
        manifest = paths.MANIFEST / f"scpr_mseed_{day}.json"
        if ("40", "LCE") in no_data and manifest.exists() and "timing_reconstructed" in json.loads(manifest.read_text()):
            no_data.remove(("40", "LCE"))           # re-timed day: no arrival stamps, so no LCE by design
            report("OK", f"{d}: re-timed day (no timestamps on piweb); LCE absent by design")
        report("FAIL" if missing_meta else "OK", f"{d}: channels with data but no metadata: {missing_meta or 'none'}")
        report("FAIL" if no_data else "OK", f"{d}: channels in metadata without data: {no_data or 'none'}")
        for key, st in written.items():
            if key not in channels:
                continue
            c = channels[key][1]
            sr = st[0].stats.sampling_rate
            if abs(sr - c.sample_rate) / c.sample_rate > 1e-3:
                report("FAIL", f"{d} {key[0]}.{key[1]}: data rate {sr} Hz vs metadata {c.sample_rate} Hz")
            first = min(tr.stats.starttime for tr in st)
            if first < c.start_date or (c.end_date and first > c.end_date):
                report("FAIL", f"{d} {key[0]}.{key[1]}: data at {first} outside the channel epoch")

    print(f"\n{'FAILED: ' + str(len(problems)) + ' problem(s)' if problems else 'All checks passed.'}")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
