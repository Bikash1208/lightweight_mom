"""Phase 20 R5-6: the 3 x 3 SMOKE-TEST package for AEDT (no predictions, no comparison): it only checks that `hfss_build_board.py` builds the model in AEDT and that the new calls
(`AssignLumpedRLC`, the plane-wave assignment) work before the real 10 x 10 board is generated.

    python scripts/research/phase20_r5_smoke_package.py [out_dir]     # default hfss_exports/phase20_r5_ris_board_smoke3x3
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_opt_cp2 import make_cell  # noqa: E402

from src.export import ris_board_export as rb  # noqa: E402
from src.ris import coding as cd  # noqa: E402
from src.ris.calibrated_library import CalibratedLibrary  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402

SMOKE_NOTE = """**SMOKE TEST ONLY -- NO PREDICTIONS, NOTHING TO COMPARE WITH.** This 3 x 3 board (36 diodes) exists to check that `hfss_build_board.py` builds the model in AEDT and that the two calls that have never
been run there work: `AssignLumpedRLC` (36 boundaries; ON sheets R only, OFF sheets C only) and the incident-wave assignment. The coding is the plain 2-bit steering coding for 30 deg in the phi = 0 plane
(any valid coding would do); it carries no design meaning at 3 x 3 cells. What to look at after the run: (1) the build finishes without the message "Lumped RLC assignment failed" or "Incident wave not assigned";
(2) HFSS > Boundaries lists `Load_1_RLC` ... `Load_36_RLC`, nine with a resistance of 5.2 ohm and 27 with a 0.025 pF capacitance, each with its current line along the stub;
(3) the excitation list shows the incident wave with the intended direction (the specular direction is theta_i = 0 here, i.e. normal incidence). Report what failed (exact message) so the script can be fixed."""


def main(out):
    cal = CalibratedLibrary.load()
    cell = make_cell([4.601469081921106, 2.475026504986856, 1.3749256942062944], 15e-3)  # the retained baseline CP-2
    coding = cd.design(cal, "plain", 10.5, 30.0, 0.0, rows=3, cols=3)
    board = rb.make_board(cell, madp_kit(), coding.states - 1, name="ris_smoke3x3", f_ghz=10.5)
    names = rb.export_board(out, board, with_step=True, predictions=SMOKE_NOTE)
    print(out, names)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "hfss_exports" / "phase20_r5_ris_board_smoke3x3")
