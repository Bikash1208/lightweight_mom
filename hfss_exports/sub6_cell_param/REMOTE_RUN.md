# QUEUE_B — run on the lab machine (session S3, LOCAL only)

This file is everything the lab session needs. There is **no connection between the lab machine and the main server**: you get the code
from git, run locally, and push results to git. Do not connect to any other machine, AEDT instance or share.

## What you run (and what the main server runs)

The BAR90 cell study (HFSS Eigenmode, periodic cell, 11.6 mm pitch) is split in two queues so nothing runs twice.

| Queue | Where | Rows |
|---|---|---|
| QUEUE_A | main server | SCREEN_4 (no-stub 90 deg + 80/100 supplement), layout benchmark, SCREEN_5 (periodic bias trace) |
| **QUEUE_B** | **this lab machine** | **`SCREEN_2.csv`** (144 rows, full cell: stub_l 0.3 / 1 / 2 mm, fr4_t, disc_r, via_drill, pad_dia, nominal) + **`SCREEN_6_short_nostub.csv`** (96) + **`SCREEN_6_short_trace.csv`** (96): island_l 3.5 / 3.0 / 2.5 / 2.0 mm x island_w 1.3 / 1.0 mm, island via at the island centroid, no-stub and trace projects |
| held | nobody yet | `SCREEN_1.csv`, `SCREEN_3.csv` (the decision session picks a subset) — do NOT run them |

QUEUE_B = 336 rows, every geometry at phases **80 / 90 / 100 deg**, states **OFF and ON**, both directions (across = x phase, along = y phase).
Each row is one HFSS Eigenmode solve (8 modes from 4.5 GHz). Expect roughly 3.5-5 min per row per process (small model, ~27 k tets).

## Fixed setup (already encoded in the scripts — do not change)

* AEDT **2025 R1 (v251)** HFSS + HPC licences. If your install is another version: `set SUB6_AEDT_VERSION=<x.y>` and **state it prominently in your notes**; the VERIFY comparison below is the gate.
* Python 3.12, `pyaedt==1.9.0`, `numpy`. Non-graphical AEDT, one AEDT process per solver process (launched by the scripts).
* Closed PEC box, **AIR bottom = 377-ohm/square impedance sheet at z = -34.01 mm (30 mm below the FR4 bottom)**, Perfect-E lid at +150 mm, lattice pairs on the four side faces. (Radiation boundaries are illegal in Eigenmode; PML is not generated.)
* **No manual mesh operations** (they inflate the model to > 1 M tets and ~10x the run time). Variables: `design_variables_bar90.csv`.
* BAR90-02LRH diode as an impedance sheet 0.30 x 0.60 mm: OFF Rs/Xs = 22.23 / -335.4 ohm/sq, ON = 1.60 / +29.2 ohm/sq.
* Every solved row stores the values read back from AEDT; a row whose read-back differs from the intended variables is refused. A failing row is **retried once** (full AEDT reconnect),
  then logged FAILED and **skipped** — it never stops the queue.

## 0. Safety rules

1. Local only. Never connect to another machine; never create or mount a share to the main server.
2. Other people may be simulating on this machine. **Do not stop, attach to or modify any process you did not start.** Only stop PIDs derived from the parent chain of a process YOU launched.
3. Size cores to the CPU actually free (see step 5). Check `Get-Process ansysedt, hf3d` and Task Manager before you start; re-check between files.
4. Never commit AEDT projects, logs, `*.aedt*`, `lab_B_work/`, other users' project names or paths. Only `results/sub6_cell_param/lab_B/**` is pushed.
5. If anything unexpected happens (licence denied, repeated AEDT crashes, VERIFY fails): STOP, write `results/sub6_cell_param/lab_B/NOTES.md` with the facts, and push that. Do not improvise model changes.

## 1. Environment (PowerShell)

```powershell
git fetch --all
git checkout -b compute/sub6-cell-param-lab-b origin/compute/sub6-cell-param   # or hub/ / github/ : the tip that contains this file
py -3.12 -m venv .venv_lab
.\.venv_lab\Scripts\python -m pip install pyaedt==1.9.0 numpy
cd hfss_exports\sub6_cell_param
..\..\.venv_lab\Scripts\python run_queue_b.py --step env
```
`--step env` writes `lab_B_work/lab_B_machine.json` (python, pyaedt, logical CPUs, whether v251 was found; no host name). Fix `ANSYSEM_ROOT251` if v251 is not found.

## 2. Build + 4-row VERIFY (hard gate, ~15-20 min)

```powershell
..\..\.venv_lab\Scripts\python run_queue_b.py --step verify
```
Builds the user's BAR50 calibration cell (closed PEC box, floor z = -33.21 mm, no mesh ops) and solves four rows (OFF/ON x along/across, phase 90 deg).
The lab VERIFY must (a) pass the user targets (f <= 1 %, Q <= 30 %) **and** (b) match the main server's VERIFY within **0.5 % in frequency** (`verify_reference/`):

| state / dir | f (GHz) | Q |
|---|---|---|
| OFF along  | 5.7746 | 154.8 |
| OFF across | 6.2030 | 193.7 |
| ON along   | 4.8662 | 130.4 |
| ON across  | 4.9314 | 87.2 |

`check_verify_vs_reference.py` decides (exit 0 = PASS). **If it fails: STOP. Do not run the queue.** Report the table in `NOTES.md` and push it (a version difference is the first suspect).

## 3. Build the three BAR90 projects (one at a time, a few minutes each)

```powershell
..\..\.venv_lab\Scripts\python run_queue_b.py --step build
```
Projects: full (stubs), `nostub` (stub 0 mm), `trace` (periodic 0.25 mm bias trace on the FR4 bottom replacing the stubs). Each must report `status BUILT`; the build summary lists `mesh_operations: []`
and boundary `ABSORBING_FLOOR_377`.

## 4. Optional layout benchmark (~15 min)

```powershell
..\..\.venv_lab\Scripts\python run_queue_b.py --step benchmark
```
Solves the same 8 real rows as 8 processes x 6 cores and as 4 x 13 cores and prints the faster one (`lab_B_work/layout_benchmark.json`). On the main server the winner was chosen this way; pick your own for this machine.

## 5. Run QUEUE_B (several hours; resumable)

```powershell
..\..\.venv_lab\Scripts\python run_queue_b.py --step queue --cores-each 8 --max-proc 16
```
`--cores-each` = cores per solver process. The launcher measures the CPU free **now**, runs `floor(free / cores_each)` processes, and re-measures before each file.
Order: SCREEN_2 (gives the stub 0.3 / 1 / 2 mm points) -> short-island no-stub -> short-island trace. It prints the licence checkout counts it finds in the Ansys client log.
**If the queue is interrupted, rerun the same command**: solved rows are skipped, rows already logged FAILED are not repeated.

Row rules: per-row hard limit 30 min; read-back checked before the solve; retry once then FAILED; every FAILED row is reported by geometry/state, and a rate > 5 % or clustering on one
state / direction / phase / geometry family is flagged as a possible physics or setup issue.

## 6. Collect and push results (no code changes!)

```powershell
..\..\.venv_lab\Scripts\python run_queue_b.py --step collect
cd ..\..
git add results/sub6_cell_param/lab_B
git commit -m "SUB6 QUEUE_B lab_B results: VERIFY + SCREEN_2 + short islands (<n> rows, <k> FAILED)"
git push hub compute/sub6-cell-param-lab-b
git push github compute/sub6-cell-param-lab-b
```
Push only to the new branch `compute/sub6-cell-param-lab-b`, never force-push, never push to `compute/sub6-cell-param`. `results/sub6_cell_param/lab_B/` contains:
`rows/*.json` (one per solved row, with the AEDT read-back), `verify/` (4-row VERIFY + comparison), `build_summary_*.json`, `parallel_report_*.json`, `FAILED_ROWS_*.json`,
`failed_report/FAILED_ROWS_REPORT.md|json`, `lab_B_machine.json`, `layout_benchmark.json`, `run_queue_b.log`, and `NOTES.md` if you wrote one.

The main-server session fetches that branch and merges the rows with `merge_lab_b.py`, which accepts a row only if the lab VERIFY is within 0.5 % of the reference, the row is in QUEUE_B, its variables equal
the CSV row, the project variant and floor depth are right, and the AEDT read-back matches — the same gates as local rows.

## Known issues (all handled by the scripts)

* PyAEDT sometimes raises `KeyError: 'HFSS Terminal Network'` after an unhealthy solve -> the row is retried once after a reconnect, then logged FAILED.
* AEDT start-up can stall for minutes under load; the launcher staggers process starts by 40 s. Kill **only** a stuck process you started.
* Do not run two `--build` steps at the same time (PyAEDT validation uses a shared temp file).
* `Radiation`/PML boundaries are not available in Eigenmode; the 377-ohm floor is the accepted absorbing approximation (floor-dependence of cell modes is about +-1 %).
