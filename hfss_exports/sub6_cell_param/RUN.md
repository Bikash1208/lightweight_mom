> **STATUS 2026-11-10 (supersedes the text below, which is the original BLOCKED draft).** The PyAEDT model now builds natively in HFSS Eigenmode
> (no mesh operations) and the user's calibration **VERIFY passes** (4/4 rows: f <= 0.64 %, Q <= 13.4 %; `verify_reference/`). The floor-artefact test is
> done and read-back verified (`ARTEFACT_TEST.csv/.json`); the 377-ohm impedance floor at -34.01 mm is the chosen absorbing floor. SCREEN is driven by
> `run_queue.py` (main server, QUEUE_A) and `run_queue_b.py` (lab machine, QUEUE_B, see `REMOTE_RUN.md`); results are merged with `merge_lab_b.py` and
> analysed by `report_screen.py`. The BAR90 section at the end of this file is current; every statement about a blocked BUILD below is historical.

# HFSS-native cell sweep — BLOCKED at BUILD (2026-11-10)

**Do not run SCREEN until the boundary conflict in SERVER_RESULT.md is resolved and VERIFY passes.**
AEDT 2025 R1 explicitly rejects the requested Radiation boundary in an Eigenmode design.
The script preserves that requested boundary and stops; it does not substitute PEC/PML.
Current package is a draft, MODEL_ONLY / NOT_HFSS_VALIDATED, not a working solver delivery.

Use CPython 3.10–3.12 and the current PyAEDT distribution, named `pyaedt` on PyPI;
its import namespace is `ansys.aedt.core`. Server: Python 3.12.16, PyAEDT 1.9.0.

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install pyaedt==1.9.0
cd hfss_exports\sub6_cell_param
..\..\.venv\Scripts\python cell_param_pyaedt.py --build
..\..\.venv\Scripts\python cell_param_pyaedt.py --verify
..\..\.venv\Scripts\python cell_param_pyaedt.py --sweep SCREEN_1.csv
..\..\.venv\Scripts\python cell_param_pyaedt.py --analyse
```

Alternatively, `python cell_param_pyaedt.py --all` builds, verifies, then runs all SCREEN CSVs and analyses.
It prints all VERIFY rows and the gate before SCREEN; f errors must be <=1%, Q errors <=30%.
Use a fresh `--project PATH` for a new build; existing projects are never overwritten.
Every solve requests all logical CPUs through documented `Hfss.analyze_setup(cores=os.cpu_count(), use_auto_settings=True)`.
This is a request, not proof of licensed/used cores: inspect license status and solver profiles.
No reduced-core retry is performed. License feature counts on this server are currently unavailable.

SCREEN files: OAT priority 1 island/slot/diode position (76 rows), 2 disc/bias hardware (56),
3 via positions/stub width (44); four state/direction rows per geometry, phase 90 only.
`*` means unique run ID, not an AEDT variable. Values are numeric in design_variables.csv units.
The first geometry is always nominal. All supplied rows pass analytic 0.15-mm fabrication preflight.
`FULL_nominal.csv` is a 12-row phase template; it is not a solved finalist or part of `--all`.
The `island_offset` variable is explicitly the longitudinal y offset. The island's x position is
derived to maintain the fixed 1.00-mm gap; these coordinates equal BASELINE at nominal.

Results go in `runs/`: all eight f/Q pairs, exact row variables, profiles, VERIFY gate and JSON ranking.
SCREEN-only analysis is provisional: it cannot determine flatness, branch identity, beta at 5.8 or anisotropy.
FULL branch tracking and field review are outstanding; no frequency-only branch assignment is accepted.
**The current draft does not yet implement final FULL branch interpolation/ranking.**
The server BUILD failure prevented checking later setup/solve/export calls in AEDT.
No HFSS numbers have been computed by this package.

---
## BAR90-02LRH update (2026-11-10, S2)

Diode: Infineon BAR90-02LRH (TSLP-2-7, Ls 0.4 nH), datasheet 2012-10-09. Footprint (p.7): copper pads
0.6 (y) x 0.35 (x) mm, 0.30 mm gap, overall 1.0 x 0.6 mm. Sheet = 0.30 (x, current) x 0.60 (y), W/L = 2.0.
5.8 GHz sheet values (MODEL_ONLY until HFSS-run): OFF 22.23 / -335.4 ohm/sq (C_eff 0.15 pF || 3 k, + 0.4 nH);
ON 1.60 / +29.2 (0.8 ohm + 0.4 nH). Sensitivity states (nominal geometry only, SENS_nominal.csv, 8 rows):
OFF_C013 Xs -390.9, OFF_C018 Xs -275.0, OFF_RP2K Rs 33.2, OFF_RP5K Rs 13.4 (derived by `bar90_sheet()`).
Datasheet typical C is 0.18-0.19 pF at 1-1.8 GHz; 0.15 pF is the isolation-curve extrapolation to 5.8 GHz.

`--diode bar50|bar90` is now mandatory (except --analyse):
- `bar50` (`design_variables_bar50.csv`): the user's calibration cell. Only this profile may run `--verify`/`--artefact`.
- `bar90` (`design_variables_bar90.csv`, to become `design_variables.csv`): SCREEN cell.
New AEDT variables `pad_w`, `pad_l` (fixed), `diode_gap` 0.30, `diode_w` 0.60 (sheet W/L = diode_w/diode_gap; the
generator rescales Rs/Xs if either is ever varied). Boundary name `BAR90_02LRH`.

Geometry (designer decision, corrected): island keeps island_w = 1.30 mm and is SHIFTED to the 0.30 mm diode gap
(x = -1.15 to +0.15, centre x = -0.50). The island via sits at the island CENTROID (x and y) in EVERY row:
x = -slot_w/2 + 0.30 + island_w/2, y = island_offset (nominal -0.50, 0). It is derived by generate_sweeps.py, so the
island-via OAT group was removed (SCREEN_3 = disc via + stub width only). SCREEN island_w = 1.0-1.6 (>= 0.9 keeps the
0.6 mm pad >= 0.15 mm inside). Pads lie wholly inside the island / patch metal (checked); they only constrain geometry.
SCREEN states: OFF (C 0.15 pF), ON, and OFF_C018 (C 0.18 pF, datasheet typical) -> 120/84/42 rows; SENS_nominal.csv (8 rows)
holds C 0.13/0.18 and Rp 2k/5k at nominal geometry. Report OFF results at both 0.15 and 0.18 pF.

SCREEN gates, all enforced BEFORE AEDT launches: PASS `--verify-gate` (BAR50 verify.json), `ARTEFACT_TEST.json` present,
and the project's build_summary.json must show the 377-ohm absorbing floor and the same `--diode`:
```
python cell_param_pyaedt.py --build --diode bar90 --floor impedance377 --air-below 30 --project <new>.aedt --output screen_bar90
python cell_param_pyaedt.py --sweep SCREEN_1.csv --diode bar90 --project <new>.aedt --output screen_bar90 `
       --verify-gate verify_cal_pec\verify.json
```
Sweep files regenerate with `python generate_sweeps.py bar90` (120/84/42 rows + SENS 8); the BAR50 profile regenerates
the old files bit-for-bit in data columns (legacy copies in `bar50_legacy/`).
