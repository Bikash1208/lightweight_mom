# Phase 27 — tensor-CP design flow in the HMSA Designer app (RIS-ready plumbing)

**Status:** active brief, written 2026-10-16. Read `docs/START_HERE.md`, `docs/STATUS.md`, `CLAUDE.md`, then
`docs/briefs/PHASE26_LEAKY_TENSOR_CP_HEXAGON.md` and `results/phase26/README.md` for what the flow is. Branch
`phase27-app-tensor-cp` from **`phase26-leaky-tensor-cp`** (the tensor-CP code is not on `main` yet). Phase 26 stays open
on its own branch only for the board A HFSS comparison and the board B export; do that work on the Phase 26 branch when
the user's files arrive, not here. All repo rules apply (one commit per item, `make test-fast` + snapshot gates green,
physics tests, honest labels, no ML in the solver core, frozen modules reused unmodified).

## Why
The user's goal is to move to RIS (Phase 20, re-activated in `docs/briefs/PHASE20_RIS_ROADMAP.md`). The Phase 26 tensor-CP
flow (slotted-circle tensor cells → supercell β/α table → pure-RHCP hologram synthesis → leaky forward model → far field →
AEDT/DXF export → HFSS comparison) lives only in `scripts/research/` and `src/dispersion/` + `src/hologram/`. The app
(`src/app/`, HMSA Designer) still reflects the 1-bit scalar flow. Putting the tensor-CP flow in the app gives the user a
validated CP design path now, and — the RIS reason — forces a **generic cell-family / state-library registry** that the
RIS design mode (roadmap R5) will plug into instead of being bolted on later.

## What already exists (check before writing code)
- Tensor extraction and families: `src/dispersion/tensor_sheet.py`, `slotted_cells.py`, `tensor_floquet.py`; supercell
  table JSONs in `results/phase26/supercell_table_*.json`; `CrossSlottedCircularPatch` geometry.
- Forward model: `src/hologram/leaky_forward.py` (`march_rays(..., density_fn=)`), `src/dispersion/leakage_model.py`,
  `src/dispersion/pin_launch.py` (HFSS-calibrated feed: R × 1.30, effective radius 1.5 mm, `HFSS_CALIBRATED`).
- Design + export scripts for `hex150_tensor_cp_single`, `hex4_tensor_cp_static` and board B (combined beam);
  comparison scripts `scripts/research/phase26_compare_single_hex.py`, `phase26_compare_fourhex.py`.
- App: FastAPI routers in `src/app/routers/` (`hmsa.py` is the existing HMSA flow, `export.py`, `pattern.py`,
  `sidelobes.py`, `hfss_compare.py`), job runner `src/app/jobs.py`, library registry `src/app/library_registry.py`,
  static UI `src/app/static/`.

## Work items (one commit each; stop point ■)

**1. Library module from the scripts (no physics change).** Move the tensor-CP design steps the scripts perform into an
importable module (e.g. `src/hologram/tensor_cp_design.py`): inputs → layout + predictions. The scripts become thin
wrappers. **Gate:** regenerating `hex150_tensor_cp_single`, `hex4_tensor_cp_static` and the board B design through the
module reproduces the committed layouts (cell list: positions, D, slot L×W, angle) **bit-identically** and the
`predictions.json` numbers exactly.

**2. Cell-family / state-library registry (the RIS-ready part).** A registry describing each cell family the app can
design with, by its *domain of validity*, not by code path: name; geometry type; substrate (εr, tanδ, h); period;
frequencies with data; what the data is (scalar Γ, 2×2 Jones/tensor, supercell β/α table, discrete switchable states);
validation label (`HFSS_VALIDATED_ONE_RUN`, `SOLVER_SUPERCELL`, `MODEL_ONLY`, …) and known biases; fabrication rules
(e.g. bridge ≥ 0.35 mm). Register: the existing 1-bit scalar pair(s), the slotted-circle tensor family. Design it so a
**discrete-state RIS library** (N states per cell, each a 2×2 reflection matrix vs f, θ, φ — roadmap R4) registers
without schema changes; include a schema test using a synthetic 4-state entry. Requests outside a family's domain
(other substrate, frequency without a table) return a named `NOT_AVAILABLE_FROM_CURRENT_DATA` state, never extrapolate.

**3. App endpoint + job.** A "Tensor CP (monopole-fed)" design mode: board = single hexagon / hex4 static (4 beams) /
hex4 combined-beam; per-port beam (θ, φ); handedness RHCP/LHCP; pin height; frequency restricted to what the registry
holds (10.54 GHz today). Runs as a background job. Outputs: cell map, predicted RHCP/LHCP patterns, AR cut, power budget
(beam / sidelobes / grazing / back / LHCP), gain/AR/XPD/SLL table — each number labelled with its source and the
registry's known biases (today: SLL optimistic by ≈ 4 dB, AR/XPD pessimistic, vs one HFSS run). The existing 1-bit HMSA
mode is unchanged (its app tests stay green untouched).

**4. Exports from the app.** DXF, layout CSV/JSON, AEDT build script (with the `DRIVE_PORT` option), `predictions.json`
and README, as the Phase 26 export folders produce today. **Gate:** app export of the single hexagon matches
`hfss_exports/hex150_tensor_cp_single/` files byte-for-byte except timestamps.

**5. HFSS comparison in the app.** Upload the user's HFSS exports (gain CSV 1°, antenna-params txt, S-parameter CSV or
.s1p/.s4p, AR cut) and get the Phase 26 comparison report (named states, power budget). Reuse the comparison scripts'
logic via the module, not a copy. Test with the committed `hfss_results/hex150_tensor_cp_single/` files: reproduces
23.03 dBi, XPD 23.0 dB, S11 −15.0 dB at 10.54 GHz.

**6. ■ Close-out.** End-to-end app test (design → export → compare) on the single hexagon. Update `docs/STATUS.md`
(capability row "tensor-CP design in the app", validation label as in Phase 26), `CHANGELOG.md`, the app section of
`docs/START_HERE.md`. Stop and report; the user then starts Phase 20 R1.

## Do not
- Change any Phase 26 physics, table or calibration; this phase moves and wires code, it does not re-model.
- Remove or alter the 1-bit flow, its presets or its tests.
- Build any RIS UI or diode model (that is R2/R5); only the registry schema must accept it.
- Run long solves in the request thread; use the job runner.
