"""Phase 12A: one-port (grounded) SParameters -- compute_s_parameters, scattering_from_fields,
stack_s_parameters, at_planes, and the surface-wave-proximity diagnostic.
"""

import numpy as np
import pytest

from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.green_function import LayeredGreenFunction
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import (
    compute_s_parameters,
    scattering_from_fields,
    stack_s_parameters,
    surface_wave_interval,
    surface_wave_proximity,
)

LATTICE = RectangularLattice(4.8e-3, 4.8e-3)
SUBSTRATE_LOSSLESS = Substrate("x", 2.2, 0.0, 1.575e-3)
SUBSTRATE_LOSSY = Substrate("x", 2.2, 0.0009, 1.575e-3)
MESH = MeshSettings(8, 8, 16)
PLANES = ReferencePlanes(5.0e-3, None, 1.575e-3)
F = 12.5e9


def grounded_cell(substrate=SUBSTRATE_LOSSLESS, length=2.0e-3, width=2.0e-3):
    return UnitCell(LATTICE, Geometry("rectangle", length, width), substrate, Metal(), backing="pec")


def test_compute_s_parameters_is_grounded_and_one_port():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    assert s.grounded is True
    assert s.s12 is None and s.s21 is None and s.s22 is None
    assert s.s11.shape == (2, 2)
    assert s.incident_amplitudes.shape == (1, 2)
    assert s.mean_currents.shape == (1, 2, 2)


def test_undefined_families_raise_not_silently_zero():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    for name in ("S21_xx", "S12_xy", "S22_yy"):
        with pytest.raises(ValueError, match="undefined"):
            s.component(name)


def test_matrix_raises_for_grounded_result():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    with pytest.raises(ValueError, match="no port 2"):
        s.matrix()


def test_lossless_patterned_grounded_cell_is_fully_reflective():
    """|S11_xx|^2 + |S11_yx|^2 = 1 for a lossless grounded stack: S21 == 0 by construction (no
    port 2 exists), so all power that goes anywhere at all comes back out of port 1."""
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSLESS), F, MESH, PLANES)
    total = abs(s.component("S11_xx")) ** 2 + abs(s.component("S11_yx")) ** 2
    assert total == pytest.approx(1.0, abs=1e-8)


def test_power_out_per_incident_shape_and_bound():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    power = s.power_out_per_incident()
    assert power.shape == (1, 2)
    assert np.all(power <= 1.0 + 1e-9)


def test_reciprocity_s11_symmetric_at_normal_incidence():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    assert np.max(np.abs(s.s11 - s.s11.T)) < 1e-9


def test_bare_grounded_slab_is_lossless_unit_reflection():
    green = LayeredGreenFunction(SUBSTRATE_LOSSLESS, F, backing="pec")
    s = stack_s_parameters(green, PLANES)
    assert s.grounded is True
    assert abs(s.component("S11_xx")) == pytest.approx(1.0, abs=1e-9)


def test_scattering_from_fields_rejects_wrong_mean_currents_shape():
    green = LayeredGreenFunction(SUBSTRATE_LOSSLESS, F, backing="pec")
    with pytest.raises(ValueError, match="shape"):
        scattering_from_fields(green, PLANES, np.zeros((2, 2, 2)))  # two-port shape, wrong here


def test_check_planes_rejects_backing_planes_mismatch():
    green_pec = LayeredGreenFunction(SUBSTRATE_LOSSLESS, F, backing="pec")
    two_port_planes = ReferencePlanes(5.0e-3, 5.0e-3, 1.575e-3)
    with pytest.raises(ValueError):
        stack_s_parameters(green_pec, two_port_planes)
    green_air = LayeredGreenFunction(SUBSTRATE_LOSSLESS, F, backing="air")
    with pytest.raises(ValueError):
        stack_s_parameters(green_air, PLANES)


def test_at_planes_grounded_rereferences_phase_but_not_magnitude():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSLESS), F, MESH, PLANES)
    new_planes = ReferencePlanes(8.0e-3, None, 1.575e-3)
    s2 = s.at_planes(new_planes)
    assert s2.magnitude("S11_xx") == pytest.approx(s.magnitude("S11_xx"), abs=1e-9)
    assert s2.phase_deg("S11_xx") != pytest.approx(s.phase_deg("S11_xx"))


def test_at_planes_rejects_mixed_topology():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSLESS), F, MESH, PLANES)
    with pytest.raises(ValueError, match="grounded"):
        s.at_planes(ReferencePlanes(8.0e-3, 5.0e-3, 1.575e-3))


def test_to_record_only_contains_s11_for_grounded():
    s = compute_s_parameters(grounded_cell(SUBSTRATE_LOSSY), F, MESH, PLANES)
    rec = s.to_record()
    assert all(k.startswith("S11_") or k == "frequency_hz" for k in rec)
    assert len(rec) == 1 + 4 * 2  # frequency_hz + 4 components * (mag, phase)


def test_full_pec_patch_makes_the_ground_plane_irrelevant():
    """A cell whose metal covers the whole period: the air side sees a PEC mirror at z=0
    regardless of what is below it, so |S11_xx| = 1 exactly as for the ungrounded full-PEC sheet
    (S_PARAMETERS.md Sec. 7: S11 = -exp(-2j k0 z_top), a pure phase factor -- only the magnitude
    is backing-independent, which is what this checks)."""
    full_pec = UnitCell(
        LATTICE,
        Geometry("rectangle", 4.8e-3 - 1e-9, 4.8e-3 - 1e-9),
        SUBSTRATE_LOSSLESS,
        Metal(),
        backing="pec",
    )
    s = compute_s_parameters(full_pec, F, MESH, PLANES)
    assert abs(s.component("S11_xx")) == pytest.approx(1.0, abs=5e-2)  # staircase-limited mesh


# --------------------------------------------------------------------- surface-wave diagnostic


def test_surface_wave_interval_has_no_cutoff_structure():
    green = LayeredGreenFunction(SUBSTRATE_LOSSLESS, F, backing="pec")
    lo, hi = surface_wave_interval(green)
    assert lo == pytest.approx(green.k0)
    assert hi == pytest.approx(green.k_substrate.real)
    assert hi > lo


def test_compute_s_parameters_warns_near_a_surface_wave_pole():
    """A deliberately large period/high-eps substrate so a Floquet harmonic falls inside the TM0
    interval; compute_s_parameters must warn (like the existing grating-order warning)."""
    big_lattice = RectangularLattice(30.0e-3, 30.0e-3)
    sub = Substrate("x", 10.2, 0.0023, 3.0e-3)
    uc = UnitCell(big_lattice, Geometry("rectangle", 5.0e-3, 5.0e-3), sub, Metal(), backing="pec")
    mesh = MeshSettings(8, 8, 16)
    planes = ReferencePlanes(5.0e-3, None, 3.0e-3)
    with pytest.warns(UserWarning, match="surface-wave"):
        compute_s_parameters(uc, F, mesh, planes, warn_grating=False)


def test_compute_s_parameters_warn_surface_wave_false_suppresses_it():
    big_lattice = RectangularLattice(30.0e-3, 30.0e-3)
    sub = Substrate("x", 10.2, 0.0023, 3.0e-3)
    uc = UnitCell(big_lattice, Geometry("rectangle", 5.0e-3, 5.0e-3), sub, Metal(), backing="pec")
    mesh = MeshSettings(8, 8, 16)
    planes = ReferencePlanes(5.0e-3, None, 3.0e-3)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_s_parameters(uc, F, mesh, planes, warn_grating=False, warn_surface_wave=False)


def test_surface_wave_proximity_reports_a_margin_for_the_production_library_stack():
    """4.8 mm period, 10-15 GHz, normal incidence: report the margin, do not assume it (brief's
    own instruction) -- this just checks the function runs and returns a self-consistent dict."""
    uc = grounded_cell(SUBSTRATE_LOSSLESS)
    for f in (10.0e9, 12.5e9, 15.0e9):
        green = LayeredGreenFunction(SUBSTRATE_LOSSLESS, f, backing="pec")
        result = surface_wave_proximity(uc, green, floquet_order=16)
        assert result["k_substrate_re"] > result["k0"]
        assert isinstance(result["any_in_band"], bool)
        assert isinstance(result["any_near_band"], bool)
        assert result["nearest_margin_over_k0"] >= 0.0
