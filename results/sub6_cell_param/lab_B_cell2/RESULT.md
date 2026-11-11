# Cell 2 result — 2026-11-11

## STATUS

COMPLETE. Six geometries × 28 rows = 168/168 eigenmode solves. Zero current FAILED rows. Local AEDT returned version `2025.2`. No scientific performance pass/fail calls were made.

## GATES

- All six analytic geometry and periodic-wrap DC gates passed; all six native AEDT builds passed validation.
- SD: `LEFT(n)` and `RIGHT(n)` are separate without the diode. Optional y lines connect `LEFT(n)↔LEFT(n+y)` and `RIGHT(n)↔RIGHT(n+y)` only.
- ID: without the diode, `LANDING(n)↔DISC(n+x)` across the x face, and the optional y line connects `DISC(n)↔DISC(n+y)` only. `DISC(n)` reaches `DISC(n+x)` through the diode.
- The shifted ID diode ends 0.20 mm inside the +x face. `ID_noline` OFF/along/70° returned eight eigenmodes on its one requested probe; that row is included in the 168.
- Every row passed exact intended-variable and AEDT readback comparison and exported eight finite, positive-Q modes. Retry-once/FAILED handling remained active.

## DATA

| Geometry | Rows |
|---|---:|
| SD_w1p0_noline | 28 |
| SD_w1p0_ylines | 28 |
| SD_w2p9_noline | 28 |
| SD_w2p9_ylines | 28 |
| ID_noline | 28 |
| ID_yline | 28 |

The phase samples are 70, 75, 80, 85, 90, 95 and 100 degrees for OFF/ON × along/across. The 168 row records contain 1,344 AEDT eigenmode entries. Results were pushed after each 28-row geometry batch.

## FLAGS

- Clearances: SD disc rim to cell face 0.365 mm, diode gap 0.30 mm, cut to yline edge 1.45 mm, minimum line-to-rim overlap 0.15 mm. Shifted ID disc left/right rim to face 0.165/0.565 mm; diode to +x face 0.20 mm; landing to local disc rim 0.365 mm; minimum yline rim overlap 0.15 mm. The 0.1267 mm tab-to-disc overlap on each ID side is an intended copper join, not an isolation clearance.
- Choices not specified in the screen: SD and ID line starts were set so the outer edge of the full 0.10 mm strip overlaps its disc rim by 0.15 mm; RT5880 uses relative permittivity 2.2 and loss tangent 0.0009 from the previous framework; copper and ground were represented as PEC at 35 µm thickness; eigenmode setup used minimum 4.5 GHz, eight modes, 3–15 adaptive passes and `MaxDeltaFreq=0.1`; free logical CPUs were redistributed before each batch across at most four workers, with a 40-second startup stagger; the successful ID probe requested 11 cores; ID's unused `split_w` bookkeeping variable remained 1.0 mm.
- The original ID geometry probe failure is retained in `evidence/`; it is not a current FAILED row. The corrected ID probe passed. Option 1b remains paused, and no G-DISC row records remain on this branch.

## FILES

- `rows/CELL2_*.json`: all 168 result/readback records.
- `geometry_dc_gates.json`: six wrapped graphs and clearance checks.
- `builds/*.json`: six native AEDT build summaries.
- `audit.json`: full row-set, readback, and mode-count audit.
- `parallel_reports/*.json`: free-CPU samples, per-worker core allocations, row and FAILED counts.
- `status.json`, `FAILED_ROWS.json`, `row_manifest.json`: current state and expected cases.
