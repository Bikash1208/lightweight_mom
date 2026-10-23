"""Phase 20 R3 item 2: the LP-1, LP-2 and CP-2 starting cells -- structure, state symmetry, energy conservation, STEP export."""

import math
import warnings
from dataclasses import replace

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.materials import Substrate
from src.ris.cells import OFF, ON, Cp2Cell, Lp1Cell, Lp2Cell, ideal_kit, madp_kit, smp1345_kit
from src.ris.polarization import R, circular_reflection
from src.sparameters import compute_s_parameters_oblique

PLANES = ReferencePlanes(0.0, None, 1.575e-3)
MESH = {"lp1": 48, "lp2": 48, "cp2": 54}


def _solve(cell, state, kit, theta=0.0, phi=0.0, lossless=False, f=10e9):
    uc = cell.unit_cell(state, kit)
    if lossless:
        uc = replace(uc, substrate=Substrate("lossless", 2.2, 0.0, 1.575e-3))
    n = MESH[cell.kind]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique(uc, f, MeshSettings(n, n, n // 2 + 1), PLANES, math.radians(theta), phi)


def test_cells_fit_in_the_period_and_have_the_stated_number_of_distinct_states():
    for cell, n_states in ((Lp1Cell(), 4), (Lp2Cell(), 4), (Cp2Cell(), 4)):
        sts = cell.states()
        assert len(sts) == len(set(sts)) == n_states
        uc = cell.unit_cell(sts[0], smp1345_kit())
        x0, x1, y0, y1 = uc.geometry.bbox()
        p = uc.lattice.period_x
        assert -p / 2 < x0 and x1 < p / 2 and -p / 2 < y0 and y1 < p / 2
        assert len(uc.loads) == cell.n_diodes
        for ld in uc.loads:  # every diode region is conductor in the pattern
            gx, gy = np.meshgrid(
                *(
                    np.linspace(a, b, 7)[1:-1]
                    for a, b in ((ld.region()[0], ld.region()[1]), (ld.region()[2], ld.region()[3]))
                ),
                indexing="ij",
            )
            keep = ld.contains(gx, gy)
            assert np.all(uc.is_metal(gx[keep], gy[keep]))
    assert Cp2Cell().states()[2] == (OFF, OFF, ON, OFF)  # one-hot
    with pytest.raises(ValueError):
        Lp1Cell().unit_cell((ON,), smp1345_kit())
    with pytest.raises(ValueError):
        Cp2Cell().unit_cell(("ON", "OFF", "OFF", "maybe"), smp1345_kit())
    with pytest.raises(ValueError):  # an axis-aligned stub reaching the cell boundary is refused
        Cp2Cell(line=3e-3).unit_cell(Cp2Cell().states()[0], smp1345_kit())


@pytest.mark.parametrize("cell", [Lp1Cell(), Lp2Cell()], ids=["lp1", "lp2"])
def test_lp_cells_give_four_distinct_states_and_leave_the_cross_polarisation_alone(cell):
    kit = smp1345_kit()
    tm = []
    te = []
    for st in cell.states():
        s = _solve(cell, st, kit).s11
        tm.append(s[1, 1])
        te.append(s[0, 0])
    for i in range(4):
        for j in range(i + 1, 4):
            assert abs(tm[i] - tm[j]) > 0.1, (i, j)  # all four combinations are distinct (asymmetric cell)
    assert (
        max(abs(te[i] - te[j]) for i in range(4) for j in range(4)) < 0.2
    )  # polarisation across the gaps is the controlled one


def test_cp2_states_are_rotations_pancharatnam_berry_ladder():
    """Ideal switches (ON short, OFF open): the four one-hot states are the cell rotated by 45 deg steps, so the handedness-preserving
    term r_RR advances by 2 x 45 = 90 deg per state. The 45-degree states are staircased on the Cartesian mesh: tolerance 6 deg."""
    cell = Cp2Cell()
    r = [circular_reflection(_solve(cell, st, ideal_kit()).s11)[R, R] for st in cell.states()]
    mags = np.abs(r)
    assert (
        mags.min() > 0.1 and mags.max() / mags.min() < 1.3
    )  # rotation-equal magnitudes (weak anisotropy of the starting geometry)
    steps = [math.degrees(np.angle(r[(k + 1) % 4] / r[k])) % 360 for k in range(4)]
    assert all(abs(s - 90.0) < 6.0 for s in steps), steps
    assert (
        abs((math.degrees(np.angle(r[2] / r[0])) % 360) - 180.0) < 3.0
    )  # exact 90-degree rotation of the axis-aligned pair


@pytest.mark.parametrize("cell", [Lp1Cell(), Lp2Cell(), Cp2Cell()], ids=["lp1", "lp2", "cp2"])
def test_energy_conservation(cell):
    # lossless substrate, reactive diodes only (all OFF, MADP with a lossless OFF): total reflected power is the incident power
    off_state = tuple(OFF for _ in range(cell.n_diodes))
    s = _solve(cell, off_state, madp_kit(0.1e-9), theta=30.0, phi=0.3, lossless=True)
    assert np.allclose(s.power_out_per_incident(), 1.0, atol=1e-8)
    # with resistive ON diodes and the real substrate: never more than the incident power
    for st in cell.states()[:2]:
        s2 = _solve(cell, st, smp1345_kit(), theta=30.0, phi=0.3)
        assert np.all(s2.power_out_per_incident() <= 1.0 + 1e-9)


def test_smp1345_and_madp_give_different_off_states():
    cell = Cp2Cell()
    st = cell.states()[0]
    a = circular_reflection(_solve(cell, st, smp1345_kit()).s11)[R, R]
    b = circular_reflection(_solve(cell, st, madp_kit()).s11)[R, R]
    assert abs(a - b) > 0.02  # the diode model matters (C_T 0.19 pF + 0.7 nH vs 0.025 pF, L_s = 0 assumed)


@pytest.mark.parametrize(
    "cell,solids,loads", [(Lp1Cell(), 3, 2), (Lp2Cell(), 3, 2), (Cp2Cell(), 5, 4)], ids=["lp1", "lp2", "cp2"]
)
def test_step_export_has_load_sheets_and_the_expected_metal_pieces(cell, solids, loads, tmp_path):
    pytest.importorskip("cadquery")
    from src.cad_export import export_unit_cell_step, import_named_parts

    uc = cell.unit_cell(cell.states()[0], smp1345_kit())
    res = export_unit_cell_step(uc, tmp_path / "cell.step")
    parts = import_named_parts(tmp_path / "cell.step")
    assert parts["Metal"]["n_solids"] == solids  # pieces separated by the diode regions cut out of the metal
    assert len(res.load_bboxes_mm) == loads and all(f"Load_{k + 1}" in parts for k in range(loads))
    if cell.kind == "cp2":  # the sheets of the diagonal stubs are rotated rectangles
        bb = res.load_bboxes_mm[1]
        assert (bb[1] - bb[0]) == pytest.approx(bb[3] - bb[2], rel=1e-6)  # square bbox of a 45-degree rectangle


def test_cp2_separate_diagonal_stub_parameters_are_additive():
    """Phase 20 R3b: `line_diag` / `stub_width_diag` default to None (identical geometry); setting them changes only the 45 / 135 deg stubs."""
    from dataclasses import replace

    base = Cp2Cell(period=15e-3, circumradius=4.6e-3, stub_width=1.37e-3, line=1.75e-3)
    same = replace(base, line_diag=base.line, stub_width_diag=base.stub_width)
    assert base.diode_sites() == same.diode_sites()
    xs = np.linspace(-7e-3, 7e-3, 141)
    gx, gy = np.meshgrid(xs, xs)
    assert np.array_equal(base.geometry().in_feature(gx, gy), same.geometry().in_feature(gx, gy))
    longer = replace(base, line_diag=3.0e-3, stub_width_diag=2.0e-3)
    sb, sl = base.diode_sites(), longer.diode_sites()
    assert sb[0] == sl[0] and sb[2] == sl[2]  # axis stubs untouched
    assert (
        sl[1][3] == 2.0e-3 and sl[3][3] == 2.0e-3 and sb[1][3] == 1.37e-3
    )  # diagonal diode width follows the diagonal stub width
    with pytest.raises(ValueError, match="diagonal stub"):
        replace(base, line_diag=9e-3).unit_cell(base.states()[0], madp_kit())
