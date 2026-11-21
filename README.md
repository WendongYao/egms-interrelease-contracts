# EGMS inter-release differences: frozen analysis contracts, code and tile-level results

Pre-registration record for the study of reprocessing differences between releases of the
European Ground Motion Service (EGMS) Ortho product (L3, vertical and east-west components).
Author: Wendong Yao (School of Computer Science, University College Dublin; ADAPT Centre).

Status: private deposit prepared on 2026-11-21, before manuscript submission.

## Rule followed in every contract
Each contract was written and its SHA-256 recorded in `README_HASH.txt` before any analysis
script for that contract was written or run. After results existed, contracts were not edited;
deviations were appended at the end of the contract text and are listed there with their time
and whether they change any number. Findings that failed to replicate are reported as such.
The contract text still reads "draft, not frozen" in its status line because the frozen
bytes are exactly the text the user approved; the freeze time is the one in `README_HASH.txt`.
Contracts and logs are written in Chinese.

## Contents
| Folder | Frozen (UTC) | Contract SHA-256 | What it tests |
|---|---|---|---|
| `contracts/Confirm01_JAG` | 2026-11-17 10:04 (Amendment 1: 13:11) | `61cd70e9…fa672a` (Amendment 1 `8b7a9609…5fbed5`) | Replication of the inter-release difference on new tiles; stability of management decisions; attribution |
| `contracts/Confirm02_Mechanism` | 2026-11-18 10:36 | `a2a07fe1…e2d15b2` | Decomposition of the difference along the processing chain (GNSS reference-frame term) |
| `contracts/Confirm03_Proximity` | 2026-11-20 12:12 | `257ab6ea…e978af1` | Threshold proximity and churn baseline; GNSS offset decomposition; geographic strata; east-west component (Part D) |
| `contracts/Confirm04_Robustness` | 2026-11-21 06:45 | `c57fbc7e…40a58` | Rounding sensitivity; spatial-block bootstrap; threshold curves; uncertainty buffer |

Each contract folder holds:
- `CONTRACT_CN.md` (+ `AMENDMENT_1_CN.md`, `_superseded/` drafts for Confirm01), `README_HASH.txt`
- `RUN_LOG_CN.md` (time-stamped run log with file hashes), `RESULT_CN.md` (results against the pre-set criteria)
- `code/` with `CODE_HASHES_frozen_before_analysis.txt` (hashes recorded before results were read)
- `results/` tile-, station- and set-level outputs with `RESULT_FILE_HASHES*.txt`
- `data/` sample definitions only: tile draws, frozen parameters, file manifests, GNSS station lists
- `logs/` extraction logs (Confirm01, Confirm02)

`tools/egms_download/` holds the download scripts. They read a CLMS service key from the
environment variable `EGMS_SERVICE_KEY_JSON`; no key, token or credential is included.

## Not included
- EGMS product files and the point data extracted from them (several GB). EGMS data are
  distributed by the Copernicus Land Monitoring Service (https://egms.land.copernicus.eu/);
  the manifests in `data/` list every tile, release and file used.
- Point-level result tables of Confirm01 (`part1_Nprime__by_point.csv`, `part1_Oprime__by_point.csv`).
  Their hashes remain listed in `results/RESULT_FILE_HASHES.txt`.
- A local list of downloaded archives with machine paths (`local_egms_zips.txt`).

## Hash check
`python verify_hashes.py` recomputes every SHA-256 listed in the frozen hash files.
Expected output on this deposit: 75 OK, 1 OK_PREFIX, 2 SUPERSEDED, 2 ABSENT.
- OK_PREFIX: Confirm01 `CONTRACT_CN.md`; the frozen hash covers lines 1-121, the text before
  the deviation entries were appended (`head -n 121 CONTRACT_CN.md | sha256sum`).
- SUPERSEDED: Confirm01 `code/s89_part2.py` (deviation 2: windowed raster reading after an
  out-of-memory stop, before any result; identical output on the check tile) and
  `code/s91_manifest_P2.py` (fix logged in `RUN_LOG_CN.md` before analysis). The current
  files match the later hashes recorded in those places.
- ABSENT: the two point-level tables listed above.
- Some hash files name results `*_N'.csv`/`*_O'.csv`; the same bytes are stored as
  `*_Nprime.csv`/`*_Oprime.csv`. Confirm03 lists inputs taken from Confirm01 `data/` and
  Confirm02 `code/rsc/`; `Confirm03_Proximity/data/relcons_by_station.csv` is the station table
  from an earlier exploratory probe and is included because Confirm03 reads it.

`MANIFEST_SHA256.txt` lists the SHA-256 of every file in this deposit.
