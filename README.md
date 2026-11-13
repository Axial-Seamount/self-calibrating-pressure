# Self-Calibrating Pressure Recorder (SCPR) — Axial Seamount

Code that turns the cabled **Self-Calibrating Pressure Recorder** at Axial Seamount central
caldera (OOI RSN, reference designator `SCPRAA301`) into MiniSEED for **OO.AXCC1**, and
corrects its two Paroscientific pressure gauges for long-term drift using the instrument's
own in-situ calibrations against a dead-weight tester.

* **Raw MiniSEED** (data quality `D`): all 23 SCPR channels, location codes 40 and 41. This,
  with the StationXML, is what goes to EarthScope.
* **Drift-corrected pressure** (data quality `Q`): `LDO` for both gauges, with a per-day drift
  uncertainty. A product users make with this repository; it is **not** sent to EarthScope.
* **StationXML** for the new channels, checked against EarthScope and the data.
* **Calibration record** (`SCPR-data-collection/calibration/`): every calibration since
  deployment (2018-07-06), the drift model and its plot. It ships with the code, so correcting
  any single day needs no archive.

All times are UTC.

---

## Contents

1. [Repository layout](#repository-layout)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Data source](#data-source)
5. [Output channels](#output-channels)
6. [Timing](#timing)
7. [Calibration and drift correction](#calibration-and-drift-correction)
8. [Daily operation (maintainers / processing VM)](#daily-operation-maintainers--processing-vm)
9. [Station metadata (StationXML)](#station-metadata-stationxml)
10. [Script reference](#script-reference)
11. [Stand-alone calibration check (steps 1–2)](#stand-alone-calibration-check-steps-12)
12. [Known issues and open questions](#known-issues-and-open-questions)
13. [References](#references)

---

## Repository layout

```
README.md                       this file
environment.yml                 conda environment (python 3.11, numpy, obspy, matplotlib)
LICENSE                         MIT
SCPR-data-collection/
  bin/                          all code (run scripts from here)
    scpr.py                     one command for users and maintainers: day / build / update
    scpr_mseed.py               piweb .dat files -> MiniSEED (one UTC day)
    scpr_retime.py              re-timing of records piweb logged without timestamps
    scpr_calibrate.py           calibration offsets, drift model, drift correction
    scpr_events.py              finds calibration days on piweb without downloading everything
    scpr_paths.py               where data is read and written (SCPR_OUTPUT)
    gap_algorithms.py           PREST gap / segment algorithms (vendored, unchanged)
    check_metadata.py           independent check of the StationXML
    create_coszo_metadata.py    StationXML from the COSZO channel spreadsheet
    download_scpr_day.py        stand-alone: one day of piweb data -> physical-unit .dat files
    extract_scpr_calibration.py stand-alone: calibration values from those .dat files
  calibration/                  committed: the calibration record used by every correction
    calibrations.csv            one row per gauge and calibration event
    calibration_days.csv        calibration days found on piweb
    drift_model.json            fitted drift model, hindcast, uncertainty model
    drift.png                   drift curves and PARO1 - PARO2 difference
  output/
    xml/                        committed: OO_AXCC1_40.xml, OO_AXCC1_41.xml
    raw/ mseed/ mseed_corrected/ manifest/ cache/ logs/   generated, not committed
  run/
    run_daily_update.sh         cron wrapper for the processing VM
```

Not in this public repository (git-ignored, kept locally): `param/` (the COSZO team channel
spreadsheet `COSZO Instrument Data Channels.xlsx`), `docs/` (instrument documentation from
G. Sasagawa) and `sample_data/` (test days for the stand-alone scripts).

## Installation

```bash
git clone https://github.com/Axial-Seamount/self-calibrating-pressure.git
cd self-calibrating-pressure
conda env create -f environment.yml        # creates env "scpr"
conda activate scpr
cd SCPR-data-collection/bin
```

By default data goes to `SCPR-data-collection/output/`. To put it elsewhere (e.g. a large disk):

```bash
export SCPR_OUTPUT=/data/scpr     # raw/, mseed/, mseed_corrected/, manifest/, cache/ go here
```

The calibration record (`SCPR-data-collection/calibration/`) always stays in the git checkout.

## Quick start

Raw and drift-corrected MiniSEED for one day or a range, straight from piweb:

```bash
python scpr.py day 2020-05-04
python scpr.py day 2020-05-01 --end 2020-05-31 --workers 4
```

Each day is downloaded (~22 MB), converted, and corrected with the committed calibration
record. Output (SDS layout, one file per channel-day):

| Folder | Content |
|---|---|
| `output/mseed/YEAR/OO/AXCC1/CHA.D/OO.AXCC1.LOC.CHA.D.YEAR.DOY` | raw, all channels, quality `D` |
| `output/mseed_corrected/...` | drift-corrected `LDO` (loc 40 and 41), quality `Q` |
| `output/mseed_corrected_masked/...` | the same with the masked calibration averaging (`--offsets masked`) |
| `output/manifest/scpr_mseed_DATE.json` | inputs (with checksums), timing, per-stream statistics |
| `output/manifest/scpr_corrected_DATE.json` | calibration table version, drift, uncertainty, provisional flag |

* **Provisional:** days after the last calibration use an extrapolated drift. They are marked
  provisional; after `git pull` brings a newer calibration record, run the same command again.
* **Not corrected:** once the drift uncertainty exceeds 100 Pa (about 200 days after the last
  calibration) the day is converted but not corrected.
* `--method linear` interpolates the drift between calibrations instead of using the fitted curve.
* `--offsets masked` uses the second calibration averaging (see
  [Calibration offset](#calibration-offset-per-event-per-gauge)) and writes to
  `output/mseed_corrected_masked/`.

## Data source

piweb (`http://piweb.ooirsn.uw.edu/scpr/data/`), **only these three folders**:

| Folder | Content | Output |
|---|---|---|
| `SCPRAA301_PARO1/` | gauge 1 (Paroscientific 44.4K, SN 136576) | loc 40 |
| `SCPRAA301_PARO2/` | gauge 2 (Paroscientific 46K, SN 132123) | loc 41 |
| `SCPRAA301_SCPR/`  | controller log: auxiliary sensors (type-1 records), calibration messages | loc 40 (LVDT2 at 41) |

The other folders on piweb (`SCPRAA301/`, `SCPR_Data/`) are not used. Files are
`YEAR/MONTH/SCPRAA301_<CH>_..._<YYYYMMDD>T<HHMM>_UTC.dat`, normally one per day; a logger
restart starts a new file (e.g. `_T2056`). Each line is normally
`<log-ts:2020-05-04T00:00:00.659> *0001,15537.993,28.767,5.854`: the arrival time at the OOI
data system, then the instrument output.

## Output channels

Network `OO`, station `AXCC1` (the existing EarthScope station; the SCPR sits 26 m from the
station point and carries its own channel coordinates: 45.9548288 N, 130.0092400 W, −1530.1 m,
from OOI's Regional Cabled Array; Glenn Sasagawa's FAQ gives 45.954907 N, 130.00928 W, 1525 m,
9 m away, which is not used).

| Loc | Channel | Content | Unit | Rate |
|---|---|---|---|---|
| 40 / 41 | `LDO` | gauge 1 / gauge 2 absolute pressure | Pa | 1 Hz |
| 40 / 41 | `LE1` | gauge pressure period | s | 1 Hz |
| 40 / 41 | `LE2` | gauge temperature period | s | 1 Hz |
| 40 / 41 | `LK1` | gauge temperature | C | 1 Hz |
| 40 / 41 | `LE3` | LVDT1 / LVDT2 voltage (piston height; powered only during calibrations) | V | 1 Hz |
| 40 / 41 | `LE9` | LVDT1 / LVDT2 piston height, 1.222 mm/V · (V − 0.25 V) (meaningful only during calibrations) | m | 1 Hz |
| 40 | `LDI` | internal barometer (Vaisala PTB110) | Pa | 1 Hz |
| 40 | `LKI` | PRT temperature (dead-weight tester) | C | 1 Hz |
| 40 | `LAX`, `LAY` | tilt X, Y | degree | 1 Hz |
| 40 | `LE4`, `LE5`, `LE7`, `LE8` | barometer, PRT, tilt X, tilt Y voltages | V | 1 Hz |
| 40 | `LE6` | mass-bell rotation | rpm | 1 Hz |
| 40 | `LCE` | instrument clock minus arrival time (per-minute median) | s | 1/60 Hz |
| 40 | `LEF` | calibration flag: 1 from "Switching SV to DWT" to "Exiting: SV to Sea" | count | 1 Hz |

Values are written as measured (no averaging, no interpolation): FLOAT64, `LEF` STEIM2,
512-byte records. Conversions follow the instrument documentation: barometer
P = 23.9836·V + 49.7607 kPa, PRT T = 1000·V − 245.8 °C, tilt = a·(V − b) with
(10.377, 2.479) and (10.312, 2.551), gauge temperature from the period with
U0/Y1–Y3 per gauge. The gauge coefficients are unchanged since July 2018.

## Timing

* **Sample times are OOI arrival times** (`<log-ts>`), UTC. The instrument's own clock runs
  1–2 minutes behind them; `LCE` records that offset.
* Exact duplicate records are dropped. After a network outage the buffered samples arrive in a
  burst; when the burst exactly fills the outage it is put back into its one-second slots.
* The daily sample interval is fitted and traces are split with PREST's `gap_algorithms`
  (`anomaly` gaps, `timing` segmenting) with a split threshold of **0.75 sample** (PREST uses
  0.5). At 0.5 every split leaves a gap of almost exactly half a sample, which MiniSEED readers
  (ObsPy/libmseed) re-join and re-time — by up to 1 s on days with clock drift. At 0.75 readers
  get back exactly the times written; every sample is within 0.75 s of its arrival time.
* Every written file is read back and checked: identical values, times within half a sample.
* Logger lines partly overwritten with NUL bytes (seen 2024-12-12) are repaired to their last
  complete record; truncated records are skipped and counted in the manifest.

### Periods without timestamps

piweb logged some data **without `<log-ts>` stamps**, so there is no arrival time to use:

| Period (UTC) | What happened |
|---|---|
| **2024-04-11 19:22 → 2024-04-17 16:19** | 5-min stream dropout; the logger came back without timestamps |
| **2026-10-17 20:56 → 2026-11-13 20:15** | 2.9-h dropout (18:04 → 20:56); the logger came back without timestamps until OOI fixed it |
| seconds to minutes on 9 other days (2018-07-06, 2021-06-30, 2022-09-08, 2023-11-03, 2023-12-21, 2024-01-30, 2024-02-05, 2024-02-16, 2025-05-02) | lines flushed from the terminal server's buffer at the start of a restart file (e.g. `_T1933`), before stamping resumed |

The instrument itself was unaffected, and these records are **re-timed** by `scpr_retime.py`
from what the data itself carries. Every re-timed day says so in its manifest
(`timing_reconstructed`: clock offset and how it was obtained, every run with its anchors, checks
and result, records re-timed and records left out).

**1. The instrument clock.** Every SCPR (auxiliary) record carries the instrument's own clock,
`UNIX_SEC`. Its offset from UTC (clock − arrival time) changes along a straight line between clock
resets (−0.555 s/day in 2026). For a day without stamps the offset is predicted from the daily
offsets (in the manifests) of the 30 nearest timestamped days on each side: **interpolated** when
there are stamps on both sides of the gap, extrapolated otherwise. On a partly stamped day it is
measured from that day's stamps. Prediction error on timestamped days: median 0.01 s. Check at the
end of the 27-day 2026 gap: −80.06 s extrapolated from before the gap, −80.32 s measured when the
stamps came back (0.26 s, i.e. ≤ 0.01 s/day).

**2. Auxiliary channels** (`LE3`–`LE9`, `LDI`, `LKI`, `LAX`, `LAY`, `LEF`): time = instrument clock −
offset, record by record, including across outages.

**3. Gauge channels** (`LDO`, `LE1`, `LE2`, `LK1`): the gauges have no clock, but write exactly one
line per second. Their lines are counted from an anchor:

| Anchor | Value | Used when |
|---|---|---|
| start of a midnight file | 00:00:00.5 (the first line arrives within the first second) | the day's file starts without stamps |
| next midnight | 23:59:59.5 | the run continues to the end of the day's last file |
| the gauge's own first real stamp after the run | that stamp − 1 s | stamps resume and the SCPR clock runs straight into its first stamped record (no outage in between) |
| the instrument clock | SCPR clock time of the run's first record | a restart file (e.g. `_T2056`) with no stamps after it |

A midnight anchor is moved by whole seconds when the instrument clock shows the logger switched
files more than 1.5 s early or late (−3 s on 2026-10-19, −4 s at the end of 2026-10-18); a file
that, by the clock, ended more than 10 s before midnight gets no end anchor. Short outages
(≤ 5 s) are placed where the instrument clock shows them.

**4. Cases and checks.** A run of gauge lines is re-timed only if it passes:

| Case | Timing | Check |
|---|---|---|
| no outage, midnight start (most days) | counted forward from the start anchor | lands within **1.5 s** of the end anchor, if there is one; else line count within 1 of the SCPR records |
| stamps resume during the day (2024-04-17, 2026-11-13) | counted forward from midnight | lands within 1.5 s of the gauge's first real stamp (2024-04-17: +0.20 / +0.26 s; 2026-11-13: −0.23 / −0.53 s) |
| restart-file buffer before stamps resume (7 of the 9 short cases) | counted **back** from the gauge's own first real stamp; each port reconnects at its own moment, so the SCPR clock is not used for gauge lines here (it was 105 s off on 2023-11-03) | SCPR clock continuous into the stamps |
| one outage > 5 s inside a midnight-start run (2026-11-08) | lines before the outage counted forward from the start anchor, lines after it counted back from the end anchor; the SCPR clock gives where it starts. Each stream drops at the outage but reconnects at its own time (2026-11-08: SCPR 27 s, PARO1 29 s, PARO2 130 s); 10 lines each side of the break are left out | the gauge outage comes out ≥ 0 s |
| more than one outage, or an outage in a restart file | — | not re-timed (2018-07-06 deployment day: 34 min of gauge data; 2022-09-08: 5 min) |

On top of that, **PARO1 and PARO2 must line up**: both gauges record the same seafloor pressure
noise (microseisms), so the cross-correlation of the two re-timed series must peak at lag 0
(± 1 s) in every 30-min window where it is clear (correlation ≥ 0.3). On 2026-11-08 all 48
windows peak at 0. In 2024 the gauges' output steps (~23 Pa) were coarser than the noise, so the
check is inconclusive there and the anchors carry the timing alone. It checks the gauges against
each other only, not their absolute time.

**5. Result.** All untimestamped days are converted. Left out: the 20 lines at the 2026-11-08
outage, 1–4 lines at the 2026-10-17/10-19 file switches, and the gauge lines on 2018-07-06 and
2022-09-08 above. `LCE` is not written for re-timed records (it would only echo the predicted
offset).

**6. Validation.** Stamps were removed from 186 timestamped days (2019–2026) — for the whole day,
and separately only up to a random time between 04:00 and 20:00 (stamps resuming during the day)
— and every re-timed line was compared with its normal time:

| | Stamps removed all day | Stamps removed until mid-day |
|---|---|---|
| Gauge day-runs re-timed / refused | 353 / 17 | 360 / 12 |
| … with **every** sample within 0.75 s | 343 (97.2 %) | 352 (97.8 %) |
| … with ≥ 99 % of samples within 0.75 s | 347 (98.3 %) | 358 (99.4 %) |
| Runs split around an outage, every sample within 0.75 s | 8 / 8 | 4 / 4 |
| Auxiliary (SCPR) days with every sample within 0.75 s | 173 / 184 (94.0 %) | 180 / 186 (96.8 %) |

The misses are a few 2021–2022 days on which the arrival stamps themselves jumped (up to 22 s on
2021-05-06, 0.03 % of samples) or the clock offset was irregular; the actual gaps (2024-04 and
2026-10/11) fall in periods where the clock offset is a clean straight line.

**For data users:** re-timed data is distributed like any other data, with no separate quality
flag in the MiniSEED or StationXML; its times are accurate to about 0.75 s (≤ 0.3 s for the
2026 gap, from the check in 1). Keep that in mind for anything timing-critical, and note that
`LCE` is absent there. If OOI can re-export these periods with real timestamps, those replace the
reconstruction.

## Calibration and drift correction

### The instrument

Every 2–8 months the controller switches both gauges from seawater to an internal dead-weight
tester (DWT; piston-cylinder PC-7300-2 SN 1466 with a 6.6344 kg mass bell) that generates a known pressure,
referenced to the air pressure inside the housing (internal barometer). A gauge's reading of that
known pressure, compared between calibrations, measures its drift.

### Finding calibrations

`scpr_events.py` finds calibration days without downloading the archive: it HEAD-requests every
daily SCPR log (a calibration day's log is not one of the usual sizes), then confirms each
candidate from the controller text (`Switching SV to DWT`, timestamped or not) or, failing
that, ≥ 600 PARO1 samples below 15 400 kPa. 2018–2020 reproduces all 22 events in the
instrument FAQ; 34 calibration days are known up to 2026-02-05
(`calibration/calibration_days.csv`).

### Calibration offset (per event, per gauge)

Method of Sasagawa et al. (2016), as in `extract_scpr_calibration.py`:

1. t0 = minimum gauge pressure within 100 s after the downward crossing of 15 400 kPa;
2. window [t0, t0 + 1200 s) in 12 bins of 100 s;
3. each sample minus the nearest internal-barometer sample;
4. per bin, samples beyond one population standard deviation are dropped, the rest averaged;
5. offset = equal-weight mean of the 12 bin means; standard error from their spread.

Results go to `calibration/calibrations.csv` (`offset_pa`, `se_pa`, barometer and PRT means,
status). The MiniSEED-based offsets agree with the stand-alone `.dat` pipeline to ≤ 0.04 Pa.

**Event QC:** events with SE > 200 Pa are `rejected: unstable` (good events: ≤ 85 Pa; e.g.
2021-05-12, when the DWT piston did not float and the gauges wandered by 50 kPa). Events failing
the coverage rules are `failed` (2018-09-21, jittery data). Deployment-day events (2018-07-06)
are excluded from the drift fit.

**Two averaging methods are kept.** Every event in the table has both:

| | `offset_pa` — **glenn** (default) | `offset2_pa` — **masked** |
|---|---|---|
| Source | Sasagawa et al. (2016), Glenn Sasagawa's extraction | the SCPR team's practice (email from Matt, Nov 2026): drop the first 1–3 min after spin-up, when the bell can wobble or stick, and remove re-spin spikes |
| Window | 1200 s from the pressure minimum | from t0 + 300 s to 60 s before "Exiting: SV to Sea" (at most t0 + 1500 s) |
| Re-spins etc. | reduced by the 1σ trim in each bin | logged re-spins, spin-motor, belt and valve operations masked (5 s before to 60 s after) |
| Bin value | 1σ-trimmed mean | median, bins weighted by samples kept |
| Drift fit rms since 2019 | 44 / 33 Pa | 37 / 29 Pa |

Both drift-corrected products can be made side by side; the masked one goes to its own folder
and never overwrites the default:

```bash
python scpr.py day 2020-05-04                    # glenn  -> output/mseed_corrected/
python scpr.py day 2020-05-04 --offsets masked   # masked -> output/mseed_corrected_masked/
```

The two differ mostly by a constant (each is relative to its own first calibration; e.g. 63 /
78 Pa on 2020-05-05). The daily update and the committed `drift_model.json` use the default.

### Drift model

* Drift of each event relative to the smooth curve at the first trusted event, plus the DWT
  thermal term `offset · 9e-6 · (T_PRT − median T_PRT)` (piston + cylinder area expansion).
* Curve: weighted least squares `a + b·t + c·exp(−t/τ)` with τ grid-searched
  (Wilcock et al. 2021); a straight line if the exponential is not significant; linear
  interpolation between events with `--method linear`. Weights `1/(SE² + 30² Pa²)`: per-event SE
  (6–20 Pa) underestimates the 30–45 Pa scatter between events.
* Current fit (31 events per gauge, 2018-07-09 → 2026-02-05): rms misfit PARO1 71 Pa,
  PARO2 40 Pa (2019 on: 44 / 33 Pa). 2018 misfits reach 70–260 Pa (fast early drift).
* **After the last calibration** the drift follows a weighted straight line through the last
  4 calibrations, tapered from the curve so it is continuous. Chosen by a leave-future-out
  hindcast (predict each calibration from the earlier ones): rms error PARO1 / PARO2 72 / 54 Pa,
  versus 112 / 103 Pa holding the last value and 116 / 58 Pa continuing the fitted curve.
* **Uncertainty** of the drift: 40 Pa inside the record, √(40² + (170 Pa/yr · Δt)²) after the
  last calibration (57 Pa at 90 days, 93 Pa at 6 months). Correction stops at 100 Pa
  (`--max-drift-sigma`), about 200 days after the last calibration.
* Before the first trusted calibration the drift is held at that value.
* `drift_model.json` holds the coefficients, every event with its residual, the hindcast and the
  uncertainty model; `drift.png` plots it.

### Correction

`corrected LDO = raw LDO − drift(t)`, written as quality `Q`. Calibrations are masked: from 60 s
before `LEF` = 1 (or pressure below 15 400 kPa) to 30 min after it ends (`--settle-seconds`).
Each corrected day records the calibration table version and drift model version used, its
drift and uncertainty, and whether it is provisional. Raw MiniSEED is never modified.

## Daily operation (maintainers / processing VM)

```bash
python scpr.py build     # first time: find calibration days, convert only those, extract every
                         #   calibration, fit the drift model -> calibration/
python scpr.py update    # daily: new calibrations, refitted model, re-correction of every local
                         #   day made with an older table/model, yesterday converted + corrected
python scpr.py update --git-push   # ... and commit + push calibration/ when it changed
```

On the processing VM, cron runs `SCPR-data-collection/run/run_daily_update.sh` (settings:
`SCPR_OUTPUT`, `SCPR_PYTHON`, `SCPR_PUSH`):

```
0 6 * * * /path/to/self-calibrating-pressure/SCPR-data-collection/run/run_daily_update.sh
```

* `update` exits non-zero when a day could not be processed and says which, so cron mail or
  monitoring notices; all other days are still processed.
* Pushing needs the VM's own git credentials (e.g. a deploy key).
* A new calibration changes the whole drift curve, so every corrected day is redone; days after
  the previous last calibration become final.
* Processing versions are recorded per day (`converter_version`, `drift_model_version`); after a
  code change that bumps them, `scpr.py` reprocesses the affected days by itself.
* The first run on an empty `SCPR_OUTPUT` scans all of piweb (HEAD requests only) and converts
  the ~34 calibration days (~1 GB). A full backfill of every day is
  `python scpr.py day 2018-07-06 --end <date> --workers 4` (~65 GB raw + ~45 GB MiniSEED;
  days already on disk are re-converted without contacting piweb).

## Station metadata (StationXML)

`output/xml/OO_AXCC1_40.xml` and `OO_AXCC1_41.xml` are generated from the AXCC1 rows of the COSZO
channel spreadsheet (not in this repository):

```bash
python create_coszo_metadata.py --metadata "../param/COSZO Instrument Data Channels.xlsx" \
    --station-name AXCC1 --location-code 40 --output ../output/xml   # and 41
python check_metadata.py          # always
```

* The **station** element is identical to EarthScope's existing AXCC1 (start 2014-07-25,
  45.954683 N, 130.008979 W, −1527 m, "RSN Axial Central Caldera 1"); the SCPR channels carry
  their own position and start (2018-07-06).
* **The spreadsheet can be wrong**, so `check_metadata.py` tests the XML against sources that are
  not the spreadsheet: the units the converter writes (`scpr_mseed.CHANNEL_UNITS`), the MiniSEED
  on disk (same channels, rates, epochs), EarthScope's live station record, and FDSN conventions
  (ASCII SI unit names, singular `count`, no blank descriptions). Non-zero exit on any failure.

## Script reference

All scripts are in `SCPR-data-collection/bin/`; `--help` lists every option.

| Command | What it does |
|---|---|
| `scpr.py day DATE [--end DATE] [--workers N] [--method fit\|linear] [--offsets glenn\|masked] [--reconvert]` | download, convert, drift-correct |
| `scpr.py build [--start] [--end]` | build the calibration record from piweb |
| `scpr.py update [--no-yesterday] [--git-push]` | daily incremental update |
| `scpr_mseed.py --date DATE [--download] [--raw-dir DIR]` | convert one day only |
| `scpr_calibrate.py extract --start --end` | calibration offsets from MiniSEED into the table |
| `scpr_calibrate.py apply --start --end` | drift-correct a date range |
| `scpr_calibrate.py model` | refit and write `drift_model.json` / `drift.png` |
| `scpr_events.py [--start] [--end]` | find calibration days on piweb |
| `check_metadata.py [--days ...] [--offline]` | check the StationXML |
| `create_coszo_metadata.py ...` | StationXML from the spreadsheet |

## Stand-alone calibration check (steps 1–2)

Two dependency-free scripts (Python ≥ 3.10, standard library only) reproduce one day's
calibration values from piweb text files, independent of the MiniSEED pipeline. Pressures in
psi, temperatures in °C. Data goes to `SCPR-data-collection/sample_data/DATE/`.

1. `download_scpr_day.py --date 2020-05-04 [--offline]` — downloads PARO1, PARO2 and the
   auxiliary log (originals in `raw/`) and writes physical-unit `*_physical.dat` files.
2. `extract_scpr_calibration.py --date 2020-05-04` — writes
   `calibration/SCPR_calibration_samples.dat` (time, pressure1_psi, pressure2_psi,
   temperature1_degC, temperature2_degC, barom_psi) and `SCPR_calibration_summary.dat`
   (time, cal_Paro1, cal_Paro2, barom). Rules:
   * threshold 2233.581161 psi (15 400 kPa, `--threshold-psi`); t0 = each gauge's minimum within
     100 s after crossing below it;
   * the first observation at or after t0 − 300 s must lie within 5 s of that time and above the
     threshold; candidates with pressure > threshold in [t0, t0 + 600 s] or ≥ threshold in
     (t0 + 600 s, t0 + 1200 s) are rejected; data must reach t0 + 1199 s;
   * coverage: no interval > 5 s within [t0, t0 + 1200 s); each 100-s bin ≥ 90 native samples;
   * end: the first upward crossing after t0 + 1200 s that stays at or above the threshold ≥ 30 s;
   * export window: 60 s before the earlier gauge t0 to 60 s after the later gauge end, times
     rounded to the second; sample barometer = nearest observation within 1 s; barometer values
     in the first 1 s after a logged SensePowerOn are excluded;
   * calibration value as in [Calibration offset](#calibration-offset-per-event-per-gauge).

## Known issues and open questions

* **Timestamps were missing 2026-10-17 20:56 → 2026-11-13 20:15 UTC** on piweb (logger issue,
  fixed by OOI). That period, 2024-04-11 → 04-17 and short restart fragments on 9 other days
  are re-timed to about 0.75 s (see [Timing](#periods-without-timestamps)); 34 min of gauge data
  on the deployment day (2018-07-06) and 5 min on 2022-09-08 cannot be re-timed.
* **No calibration since 2026-02-05**: data after about 2026-08-21 is converted but not
  drift-corrected (uncertainty > 100 Pa).
* **LVDT coefficients** (a = 1.222 mm/V, b = 0.25 V for both LVDTs) come from Glenn Sasagawa's
  processing and will be re-checked after the 2027 recovery. Days converted before converter
  version 4 have no `LE9`; `scpr.py day START --end END --reconvert` adds it from the raw files.
* The documented DWT mass and piston area (6.6344 kg, 4.9034 mm²) do not reproduce the measured
  calibration pressure (~15 269 kPa); this does not affect the drift, which uses differences.
* The internal barometer reading fell 0.29 kPa from 2018 to 2026 (mostly in 2018); housing
  pressure or sensor drift — it enters the offsets directly.
* Early drift (2018) is faster than the model follows; use `--method linear` to pass exactly
  through every calibration.

## References

* Sasagawa et al. (2016) — the self-calibrating pressure recorder method and the calibration
  offset algorithm used here (12 × 100-s bins, 1σ trim).
* Wilcock et al. (2021) — the exponential-plus-linear form used for the drift curve.
* (Full citations: see the instrument documentation; to be completed.)
* `gap_algorithms.py` is vendored unchanged from the PREST pipeline
  (absolute-seafloor-pressure, commit e2c34f3).
* Instrument documentation by G. Sasagawa (Scripps): data-file readme, FAQ (deployment, serial
  numbers, calibration events), supplementary materials.

## License

MIT — see [LICENSE](LICENSE).
