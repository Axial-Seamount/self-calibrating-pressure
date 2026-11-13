"""Where the SCPR pipeline reads and writes.

Layout: <repo>/SCPR-data-collection/{bin, calibration, output, run}; this file is in bin/.
CALIBRATION is part of the git checkout (the table and drift model shared with users).
OUTPUT holds everything generated and large (raw downloads, MiniSEED, manifests, cache);
it defaults to SCPR-data-collection/output and can be moved with the SCPR_OUTPUT environment variable,
e.g. on the processing VM:  export SCPR_OUTPUT=/data/scpr
"""

import os
from pathlib import Path

BIN = Path(__file__).resolve().parent
ROOT = BIN.parent                         # SCPR-data-collection
CALIBRATION = ROOT / "calibration"
OUTPUT = Path(os.environ.get("SCPR_OUTPUT") or ROOT / "output").expanduser().resolve()

RAW = OUTPUT / "raw"                      # piweb .dat files, one folder per day
MSEED = OUTPUT / "mseed"                  # raw MiniSEED (quality D), SDS
CORRECTED = OUTPUT / "mseed_corrected"    # drift-corrected LDO (quality Q), SDS
MANIFEST = OUTPUT / "manifest"            # per-day JSON records
REPORT = OUTPUT / "calibration"           # per-run apply reports
CACHE = OUTPUT / "cache"                  # piweb file-size scan
