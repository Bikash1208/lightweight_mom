"""Phase 28: the HFSS run of the loop cell R4 (hfss_results/phase28_loop_R4) against the MoM. The HFSS project's state labels are swapped relative to the physics (OFF.csv carries the 2.4 ohm boundary: lossy,
isotropic; ON.csv the 28.2 fF one: open, anisotropic); with the physical pairing every sheet-admittance component agrees within 5 % (loss within 25 %), with the project labels it does not."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
needs = pytest.mark.skipif(not (ROOT / "hfss_results/phase28_loop_R4/ON.csv").exists() or not (ROOT / "results/phase28/f10/loop_R4/measured_rows.jsonl").exists(), reason="HFSS csv / MoM rows not in the tree")


@needs
def test_physical_pairing_agrees_with_the_mom_and_the_project_labels_do_not():
    import phase28_hfss_cell_compare as c

    fo, so = c.read_hfss("ON.csv")
    ff, sf = c.read_hfss("OFF.csv")
    yc, yr = c.sheet_y(fo, so), c.sheet_y(ff, sf)
    mom = c.mom_states()
    for f in (9.5, 10.0, 10.5):
        i = int(np.argmin(np.abs(fo - f)))
        rx, ry, _ = mom[("ON_10mA_Ronly", f, 60)]
        cx, cy, _ = mom[("OFF_0V_Conly", f, 60)]
        rel = lambda a, b: abs(a - b) / abs(b)  # noqa: E731
        # mode 2 = x (along the diode), mode 1 = y
        assert rel(yr[i, 1].imag, rx.imag) < 0.05 and rel(yr[i, 0].imag, ry.imag) < 0.05
        assert rel(yc[i, 1].imag, cx.imag) < 0.06 and rel(yc[i, 0].imag, cy.imag) < 0.05
        assert rel(yr[i, 1].real, rx.real) < 0.25  # the R loss
        # project labels (ON.csv taken as the R state) fail by a large margin
        assert rel(yc[i, 1].imag, rx.imag) > 0.5


def test_calibration_hook():
    from src.ris import loop_cell_calibration as lc

    yr, yc = lc.hfss_tensors(10.0)
    assert yr[0, 0].imag * 1e3 == pytest.approx(6.54, abs=0.05) and yc[0, 0].imag * 1e3 == pytest.approx(2.33, abs=0.05)
    assert (yr[0, 0].imag - yc[0, 0].imag) * 1e3 == pytest.approx(4.22, abs=0.05)  # Delta B along the diode: the advisor's 4.22 mS
    assert yr[0, 0].real * 1e3 == pytest.approx(0.105, rel=0.05) and yc[0, 0].real * 1e3 < 0.01  # the loss sits in the R state, along x
    with pytest.raises(ValueError):
        lc.hfss_tensors(12.0)
    d = lc.delta_to_mom(10.0)
    assert np.allclose(lc.calibrate(np.zeros((2, 2)), "R", 10.0), d["R"])
    assert 2.6 < lc.alpha_exact(10.0) < 3.1
