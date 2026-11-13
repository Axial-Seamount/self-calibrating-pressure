#!/usr/bin/env python3
"""Download one UTC day of SCPR data and write unaveraged physical-unit .dat files.

Uses only the Python standard library. Default: calibration event 20,
2020-05-04 (FAQ event time 16:08:00 UTC). No calibration extraction or correction.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import date, datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import urljoin, urlsplit
from urllib.request import urlopen


BASE_URL = "http://piweb.ooirsn.uw.edu/scpr/data/"
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "sample_data"   # SCPR-data-collection/sample_data
KPA_PER_PSI = 6.894757293168361
CHANNELS = ("PARO1", "PARO2", "SCPR")
# Supplementary Materials Axial Cabled.docx: original SCPR unit 001 gauges.
GAUGES = {
    "PARO1": {"id": "0001", "serial": "136576", "U0": 5.856127,
              "Y1": -3969.862, "Y2": -9294.779, "Y3": 0.0},
    "PARO2": {"id": "0002", "serial": "132123", "U0": 5.853962,
              "Y1": -3917.198, "Y2": -9346.861, "Y3": 0.0},
}
LOG_RECORD = re.compile(rb"^<log-ts:([^>]+)>\s*(.*)$")


class DirectoryLinks(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            self.links.extend(value for key, value in attrs if key == "href" and value)


def source_pattern(channel: str, day: date) -> str:
    return rf"SCPRAA301_{channel}_[^/]+_{day:%Y%m%d}T\d{{4}}_UTC\.dat"


def discover_files(base_url: str, channel: str, day: date) -> list[str]:
    directory = f"{base_url.rstrip('/')}/SCPRAA301_{channel}/{day:%Y/%m}/"
    with urlopen(directory, timeout=60) as response:
        parser = DirectoryLinks()
        parser.feed(response.read().decode("utf-8"))
    # Only direct data-file links; the logger IP/port need not be hard-coded.
    urls = sorted({urljoin(directory, name) for name in parser.links
                   if re.fullmatch(source_pattern(channel, day), name)})
    if not urls:
        raise ValueError(f"No {channel} files for {day} in {directory}")
    return urls


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        print(f"Reuse: {path.name}", flush=True)
        return
    temporary = path.with_suffix(path.suffix + ".part")
    try:
        with urlopen(url, timeout=60) as response, temporary.open("wb") as target:
            expected_size = response.headers.get("Content-Length")
            shutil.copyfileobj(response, target)
        actual_size = temporary.stat().st_size
        if not actual_size or (expected_size and actual_size != int(expected_size)):
            raise ValueError(f"Incomplete download: {url} ({actual_size} bytes)")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Downloaded: {path.name} ({actual_size:,} bytes)", flush=True)


def finite_numbers(values: list[str]) -> list[float]:
    numbers = [float(value) for value in values]
    if not all(math.isfinite(value) for value in numbers):
        raise ValueError("non-finite measurement")
    return numbers


def convert_payload(channel: str, payload: str) -> list[float | int] | None:
    """None denotes a header/command/diagnostic, never a calibration selection."""
    if channel in GAUGES:
        if not payload.startswith("*"):
            return None
        parts = payload[1:].split(",")
        coeff = GAUGES[channel]
        if len(parts) != 4 or parts[0] != coeff["id"]:
            raise ValueError(f"unexpected {channel} gauge identifier or field count")
        pressure_kpa, _pressure_period, temperature_period = finite_numbers(parts[1:])
        u = temperature_period - coeff["U0"]
        temperature = coeff["Y1"] * u + coeff["Y2"] * u**2 + coeff["Y3"] * u**3
        return [pressure_kpa / KPA_PER_PSI, temperature]

    parts = payload.split()
    if not parts or parts[0] != "1":
        return None
    # Operator menus are echoed into the log too (e.g. "1    Auto Calibration ...");
    # only 9 numeric fields make a type-1 measurement.
    if len(parts) != 9:
        return None
    try:
        system_unix = int(parts[1])
        lvdt1, lvdt2, baro, prt, rpm, xtilt, ytilt = finite_numbers(parts[2:])
    except ValueError:
        return None
    # Readme - Cabled SCPR Data files.docx; LVDT a (mm/V) and b (V) from Glenn Sasagawa.
    return [system_unix, (23.9836 * baro + 49.7607) / KPA_PER_PSI,
            1000.0 * prt - 245.8, rpm, 10.377 * (xtilt - 2.479),
            10.312 * (ytilt - 2.551), 1.222 * (lvdt1 - 0.25), 1.222 * (lvdt2 - 0.25)]


def write_header(handle, channel: str, day: date, sources: list[dict]) -> None:
    handle.write(f"# SCPR unit 001 | {channel} | UTC day {day}\n")
    handle.write("# time_utc: OOI data-center arrival time, not instrument integration time.\n")
    handle.write("# Native samples/order retained, including repeated timestamps; no averaging,\n")
    handle.write("# interpolation, calibration extraction, drift correction or barometer subtraction.\n")
    handle.write(f"# Absolute pressure in psi; 1 psi = {KPA_PER_PSI} kPa.\n")
    for source in sources:
        handle.write(f"# source: {source['url']} | sha256: {source['sha256']}\n")
    if channel in GAUGES:
        coeff = GAUGES[channel]
        handle.write(f"# Gauge S/N {coeff['serial']}; pressure uses the logged controller kPa value.\n")
        handle.write("# Temperature: U=raw_temperature_period_us-U0; T=Y1*U+Y2*U^2+Y3*U^3.\n")
        handle.write(f"# U0={coeff['U0']}; Y1={coeff['Y1']}; Y2={coeff['Y2']}; Y3={coeff['Y3']}.\n")
        columns = "time_utc pressure_psi temperature_degC"
    else:
        handle.write("# Type 1 sensor records only; header/command/diagnostic text remains in raw files.\n")
        handle.write("# barometer_psi=(23.9836*Vbaro+49.7607)/6.894757293168361\n")
        handle.write("# prt_temperature_degC=1000*Vprt-245.8\n")
        handle.write("# tilt_x_deg=10.377*(Vx-2.479); tilt_y_deg=10.312*(Vy-2.551)\n")
        handle.write("# lvdt1_mm=1.222*(V1-0.25); lvdt2_mm=1.222*(V2-0.25) (piston height)\n")
        handle.write("# LVDTs are powered only during calibrations; otherwise (<0.1 V) not piston height.\n")
        handle.write("# system_unix_s is the drifting instrument clock.\n")
        columns = ("time_utc system_unix_s barometer_psi prt_temperature_degC "
                   "rotation_rpm tilt_x_deg tilt_y_deg lvdt1_mm lvdt2_mm")
    handle.write("# " + "\t".join(columns.split()) + "\n")


def convert_files(channel: str, day: date, sources: list[dict], output: Path) -> dict:
    counts: Counter = Counter()
    seen_times: set[str] = set()
    first_time = last_time = None
    temporary = output.with_suffix(".dat.part")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            write_header(handle, channel, day, sources)
            for source in sources:
                path = Path(source["path"])
                with path.open("rb") as raw:
                    for line_number, line in enumerate(raw, 1):
                        counts["input_lines"] += 1
                        match = LOG_RECORD.match(line.rstrip(b"\r\n"))
                        if match is None:
                            counts["non_record_lines"] += 1
                            continue
                        try:
                            values = convert_payload(channel, match[2].decode("ascii"))
                            if values is None:
                                counts["header_command_or_diagnostic_lines"] += 1
                                continue
                            timestamp = match[1].decode("ascii")
                            instant = datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
                            if instant.tzinfo is None:
                                instant = instant.replace(tzinfo=timezone.utc)
                            instant = instant.astimezone(timezone.utc)
                            if instant.date() != day:
                                counts["outside_requested_utc_day"] += 1
                                continue
                        except (ValueError, OverflowError) as exc:
                            raise ValueError(f"{path.name}:{line_number}: {exc}") from exc
                        timestamp = instant.isoformat().replace("+00:00", "Z")
                        if timestamp in seen_times:
                            counts["repeated_timestamps_kept"] += 1
                        seen_times.add(timestamp)
                        first_time = timestamp if first_time is None else min(first_time, timestamp)
                        last_time = timestamp if last_time is None else max(last_time, timestamp)
                        # Python float remains double precision, including the temperature period.
                        handle.write(timestamp + "\t" + "\t".join(format(v, ".17g") for v in values) + "\n")
                        counts["output_rows"] += 1
        if not counts["output_rows"]:
            raise ValueError(f"No valid {channel} measurements on {day}")
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    return {"output": str(output), "first_time_utc": first_time,
            "last_time_utc": last_time, "counts": dict(counts), "sources": sources}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", type=date.fromisoformat, default=date(2020, 5, 4),
                        help="UTC day YYYY-MM-DD (default: 2020-05-04, FAQ event 20)")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT,
                        help="Output root; a YYYY-MM-DD subdirectory is created")
    parser.add_argument("--base-url", default=BASE_URL)
    parser.add_argument("--offline", action="store_true",
                        help="Convert the existing day/raw/*.dat files without network access")
    args = parser.parse_args()
    day_dir = args.out_dir.expanduser().resolve() / args.date.isoformat()
    raw_dir = day_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"date_utc": args.date.isoformat(), "kpa_per_psi": KPA_PER_PSI,
                "processing": "native samples; no calibration extraction or correction",
                "documentation": ["Readme - Cabled SCPR Data files.docx",
                                  "Supplementary Materials Axial Cabled.docx", "CSCPR FAQS.docx"],
                "channels": {}}
    for channel in CHANNELS:
        if args.offline:
            urls = [f"{args.base_url.rstrip('/')}/SCPRAA301_{channel}/{args.date:%Y/%m}/{p.name}"
                    for p in sorted(raw_dir.glob("*.dat"))
                    if re.fullmatch(source_pattern(channel, args.date), p.name)]
            if not urls:
                raise ValueError(f"No cached {channel} data for {args.date} in {raw_dir}")
        else:
            urls = discover_files(args.base_url, channel, args.date)
        sources = []
        for url in urls:
            path = raw_dir / Path(urlsplit(url).path).name
            if not args.offline:
                download(url, path)
            sources.append({"url": url, "path": str(path), "bytes": path.stat().st_size,
                            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        output = day_dir / f"SCPR_{channel}_{args.date:%Y%m%d}_physical.dat"
        report = convert_files(channel, args.date, sources, output)
        manifest["channels"][channel] = report
        print(f"{channel}: {report['counts']['output_rows']:,} rows -> {output}", flush=True)
    manifest_path = day_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Processing record: {manifest_path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
