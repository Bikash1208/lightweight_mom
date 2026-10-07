"""Phase 13 Part 2: solve diagnostics, on a small/cheap real solve (this module wraps
`sparameters.compute_s_parameters` unmodified -- the physics is already covered by that module's
own tests; these check the diagnostics themselves)."""

import pytest

from src import diagnostics as diag
from src import errors
from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice


def _small_ungrounded_cell():
    lattice = RectangularLattice(4.8e-3, 4.8e-3)
    substrate = Substrate("sub", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    metal = Metal("pec")
    geometry = Geometry(kind="rectangle", length=2.0e-3, width=2.0e-3)
    return UnitCell(lattice, geometry, substrate, metal)


def test_instrumented_solve_returns_result_and_diagnostics():
    cell = _small_ungrounded_cell()
    mesh = MeshSettings(cells_x=8, cells_y=8, floquet_order=8)
    planes = ReferencePlanes(
        top_air_thickness=10e-3, bottom_air_thickness=10e-3, substrate_thickness=cell.substrate.thickness
    )
    result, d = diag.instrumented_compute_s_parameters(cell, 10e9, mesh, planes, run_id="test-run")

    assert result.magnitude("S11_xx") <= 1.0 + 1e-9
    assert d.run_id == "test-run"
    assert d.cells_x == 8 and d.cells_y == 8
    assert d.n_unknowns > 0
    assert d.total_time_s > 0.0
    assert d.condition_estimate is not None  # small mesh, below CONDITION_ESTIMATE_MAX_UNKNOWNS
    assert d.energy_conservation_residual is not None and d.energy_conservation_residual < 0.2
    assert d.reciprocity_residual is not None and d.reciprocity_residual < 1e-9  # exact reciprocity
    assert d.surface_wave_margin is None  # ungrounded stack
    # PHASE13B A9: platform-normalized to kB everywhere -- a tiny test-process solve is nowhere
    # near the ~1.4e9 "kB" the live review saw from an un-normalized macOS ru_maxrss (bytes).
    assert d.peak_memory_kb is not None and 0.0 < d.peak_memory_kb < 1.0e7


def test_ru_maxrss_kb_normalizes_macos_bytes_but_not_linux_kilobytes(monkeypatch):
    """PHASE13B A9: POSIX leaves ru_maxrss's unit unspecified -- Linux reports kB, macOS/BSD report
    bytes. `_ru_maxrss_kb` must divide by 1024 only on the platforms that need it."""

    class _FakeUsage:
        ru_maxrss = 1_400_000_000  # the live review's own reported (buggy, un-normalized) value

    monkeypatch.setattr(diag.resource, "getrusage", lambda who: _FakeUsage())

    monkeypatch.setattr(diag.sys, "platform", "darwin")
    assert diag._ru_maxrss_kb() == pytest.approx(1_400_000_000 / 1024.0)

    monkeypatch.setattr(diag.sys, "platform", "linux")
    assert diag._ru_maxrss_kb() == pytest.approx(1_400_000_000.0)


def test_condition_estimate_skipped_above_threshold(monkeypatch):
    monkeypatch.setattr(diag, "CONDITION_ESTIMATE_MAX_UNKNOWNS", 1)
    cell = _small_ungrounded_cell()
    mesh = MeshSettings(cells_x=8, cells_y=8, floquet_order=8)
    planes = ReferencePlanes(
        top_air_thickness=10e-3, bottom_air_thickness=10e-3, substrate_thickness=cell.substrate.thickness
    )
    _, d = diag.instrumented_compute_s_parameters(cell, 10e9, mesh, planes)
    assert d.condition_estimate is None
    assert "too expensive" in d.condition_estimate_skipped_reason


def test_solve_failure_raises_typed_lwmom_error(monkeypatch):
    """A low-level solver failure (any exception from `compute_s_parameters`) must come out as a
    typed `LwmomError` with the original text preserved, not propagate raw or vanish."""

    def _boom(*args, **kwargs):
        raise RuntimeError("singular impedance matrix")

    monkeypatch.setattr(diag, "compute_s_parameters", _boom)
    cell = _small_ungrounded_cell()
    mesh = MeshSettings(cells_x=8, cells_y=8, floquet_order=8)
    planes = ReferencePlanes(
        top_air_thickness=10e-3, bottom_air_thickness=10e-3, substrate_thickness=cell.substrate.thickness
    )
    with pytest.raises(errors.LwmomError) as exc_info:
        diag.instrumented_compute_s_parameters(cell, 10e9, mesh, planes)
    assert exc_info.value.code == "SOLVE_FAILED"
    assert exc_info.value.debug_detail == "singular impedance matrix"
