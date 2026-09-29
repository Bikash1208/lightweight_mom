"""Phase 7 sweep construction and solve orchestration: `runner.run_sweep` (Phase 5, unmodified)
reused exactly as `library_builder.py` wires it, at a small/tractable scale -- these tests never
build the full 324x18x18x201 production library (see `examples/build_unit_cell_library.py` for
that, run separately and explicitly)."""

import warnings

import pytest

from src.basis_functions import MeshSettings
from src.discretization import DEFAULT_MAX_AREA_ERROR, check_cell
from src.geometry import Geometry, UnitCell
from src.library_builder import build_library_database, build_sweep_definition
from src.library_schema import DEFAULT_MESH, LibraryConfig
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def small_config(l_values_mm=(1.0, 2.4, 4.2), w_values_mm=(1.0, 2.4, 4.2), **overrides):
    kwargs = dict(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=12.0e9,
        f_max_hz=13.0e9,
        f_step_hz=0.5e9,
        f_design_hz=12.5e9,
        l_values_mm=l_values_mm,
        w_values_mm=w_values_mm,
    )
    kwargs.update(overrides)
    return LibraryConfig(**kwargs)


# ------------------------------------------------------------------------- sweep construction


def test_build_sweep_definition_is_length_width_cartesian_product():
    cfg = small_config()
    definition = build_sweep_definition(cfg)
    names = [a.name for a in definition.axes]
    assert names == ["length_m", "width_m"]
    assert definition.axes[0].values == tuple(v * 1e-3 for v in cfg.l_values_mm)
    assert definition.axes[1].values == tuple(v * 1e-3 for v in cfg.w_values_mm)
    assert definition.on_invalid == "skip"  # req. 4: record, don't abort
    assert definition.validation.strict is False  # high-area-error cells still solve


def test_build_sweep_definition_frequency_grid():
    cfg = small_config()
    definition = build_sweep_definition(cfg)
    freqs = definition.frequency.frequencies()
    assert freqs[0] == pytest.approx(cfg.f_min_hz)
    assert freqs[-1] == pytest.approx(cfg.f_max_hz)


# ------------------------------------------------------------------------------- solve


def test_build_library_database_solves_every_geometry():
    cfg = small_config()
    result = build_library_database(cfg, out_dir=None)
    assert result.n_ok_cells == 9
    assert result.n_invalid_cells == 0
    assert len(result.database.cells) == 9
    assert len(result.database.records) == 9 * 3  # 3 frequencies


def test_every_ok_cell_solves_at_every_frequency():
    cfg = small_config()
    result = build_library_database(cfg, out_dir=None)
    n_freq = len(result.database.frequencies)
    for cell in result.database.cells:
        assert cell.status == "ok"
        assert len(result.database.records_of(cell.index)) == n_freq


def test_build_raises_if_geometry_count_does_not_match_config():
    """A belt-and-braces invariant (req. 4): the builder itself checks its own output count."""
    cfg = small_config()
    result = build_library_database(cfg, out_dir=None)
    assert len(result.database.cells) == cfg.n_geometries


# ---------------------------------------------------------------------------- reproducibility


def test_two_builds_of_the_same_config_are_numerically_identical(tmp_path):
    """Req. 27: no random ordering, no hidden global state -- repeated builds of the same small
    config give the same cells in the same order with the same S-matrices."""
    cfg = small_config()
    r1 = build_library_database(cfg, out_dir=tmp_path / "run1", resume=False)
    r2 = build_library_database(cfg, out_dir=tmp_path / "run2", resume=False)
    assert [c.cell_id for c in r1.database.cells] == [c.cell_id for c in r2.database.cells]
    assert r1.database == r2.database


def test_resume_reuses_checkpoint_without_resolving(tmp_path):
    cfg = small_config()
    out = tmp_path / "resumable"
    r1 = build_library_database(cfg, out_dir=out, resume=True)
    assert r1.computed == 9 * 3 and r1.reused == 0

    r2 = build_library_database(cfg, out_dir=out, resume=True)
    assert r2.computed == 0 and r2.reused == 9 * 3
    assert r1.database == r2.database


# ------------------------------------------------------------ mesh choice: documented pitfall


def _rect_cell(length_mm, width_mm, cells, floquet_order):
    sub = Substrate(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    lat = RectangularLattice(period_x=4.8e-3, period_y=4.8e-3)
    geo = Geometry(kind="rectangle", length=length_mm * 1e-3, width=width_mm * 1e-3)
    uc = UnitCell(lattice=lat, geometry=geo, substrate=sub, metal=Metal())
    mesh = MeshSettings(cells_x=cells, cells_y=cells, floquet_order=floquet_order)
    return uc, mesh


def test_24x24_mesh_has_large_area_error_for_some_sweep_geometries():
    """Documents, with a runnable check (not just prose), why `library_schema.DEFAULT_MESH` is
    32x32 and not the seemingly-ideal 24x24 (whose 0.2 mm cell pitch exactly matches the L/W
    sweep step): a 1.0x1.0 mm patch at 24x24 misclassifies an entire extra mesh row/column due to
    floating-point rounding in the (unmodified) mesh boundary comparison, giving a metal area
    error far above this project's own `DEFAULT_MAX_AREA_ERROR` (5%) policy -- not an ordinary,
    small staircase error."""
    uc, mesh = _rect_cell(1.0, 1.0, cells=24, floquet_order=32)
    chk = check_cell(uc, mesh, DEFAULT_MAX_AREA_ERROR)
    assert chk["metal_area_error"] > 0.30  # measured ~44%, nowhere near an ordinary staircase


def test_32x32_mesh_has_no_such_pathology_for_the_same_geometry():
    """The chosen mesh (`DEFAULT_MESH`) does not exhibit the 24x24 pathology above: its area
    error for the same worst-case (smallest) patch is large only in the ordinary,
    proportional-to-dx/L staircase sense, not a multiple-row misclassification."""
    uc, mesh = _rect_cell(1.0, 1.0, cells=DEFAULT_MESH["cells_x"], floquet_order=DEFAULT_MESH["floquet_order"])
    chk = check_cell(uc, mesh, DEFAULT_MAX_AREA_ERROR)
    # dx = 4.8/32 = 0.15 mm on a 1.0 mm patch: an ordinary staircase error is a few dx/L, not the
    # ~9x-worse jump measured for 24x24 (0.44) on the same geometry.
    assert chk["metal_area_error"] < 0.25
