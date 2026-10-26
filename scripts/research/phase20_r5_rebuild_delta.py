"""Phase 20 R5: rebuild the run-2 HFSS delta table (HFSS - library, structure frame) against the full-grid library v2. Same method, same HFSS data, same plane correction as
`phase20_r3b_rescore.py` (run-1 plane from the bare-slab fit minus the data-measured run-to-run shift, `PLANE_CORRECTION.md`); only the library the delta is taken against changes.
The user approved the method for this cell on 2026-10-23 (HFSS_CALIBRATED); this is that method applied to the new library, recorded in the table's source string.

    python scripts/research/phase20_r5_rebuild_delta.py [out.npz]
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.hfss_calibration import DeltaTable, fit_reference_plane_mm, read_hfss_reim, run_shift_mm  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402

HF = ROOT / "hfss_results" / "phase20_r3_cp_cell"


def main(out):
    lib = StateLibrary.load(ROOT / "results/phase20_r4/cp2_library_v2_fullgrid.npz")
    d1, _ = fit_reference_plane_mm(HF / "hfss_bare_slab_theta0_30.csv", ROOT / "hfss_exports/phase20_r3_cp_cell/bare_grounded_slab_reference.csv")
    r1 = read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30.csv", 0.0)
    shift, lo, hi = run_shift_mm(r1, read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30_run2_fine.csv", 0.0))
    d_res = d1 - shift
    r2 = read_hfss_reim(HF / "hfss_cp2_states1-4_theta0_30_run2_fine.csv", d_res)
    src = (
        f"hfss_cp2_states1-4_theta0_30_run2_fine.csv (user HFSS run 2); plane {d_res:.3f} mm = run-1 plane {d1:.3f} - measured shift {shift:.3f} mm (16-84 % {lo:.3f}..{hi:.3f}); "
        "delta against library v2 (full grid, results/phase20_r4/cp2_library_v2_fullgrid.npz); cross sign +1; method approved by the user 2026-10-23, rebuilt 2026-10-24 for the new library"
    )
    DeltaTable.build(lib, r2, src, approved_by_user=True).save(out)
    print("wrote", out, src)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "results/phase20_r5/hfss_delta_table_run2_lib_v2.npz")
