"""STEP export of a solved library (req. 16-17): deterministic filenames, independent bounds
verification, and every produced solid re-validated as a genuine topological solid on re-import
-- reusing `cad_export.py` (Phase 6, unmodified), not a second geometry pipeline.

Needs the optional cadquery dependency; skipped entirely otherwise.
"""

import warnings

import pytest

cadquery = pytest.importorskip("cadquery")

import math  # noqa: E402

from src.cad_export import solid_validity  # noqa: E402
from src.geometry import Geometry, UnitCell  # noqa: E402
from src.geometry_export import check_geometry_bounds, export_library_geometries, step_filename  # noqa: E402
from src.library_builder import build_library_database  # noqa: E402
from src.library_schema import LibraryConfig  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def _cell(length_mm, width_mm, period_mm=4.8, rotation_deg=0.0, center=(0.0, 0.0)):
    sub = Substrate(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    lat = RectangularLattice(period_x=period_mm * 1e-3, period_y=period_mm * 1e-3)
    geo = Geometry(
        kind="rectangle",
        length=length_mm * 1e-3,
        width=width_mm * 1e-3,
        rotation=math.radians(rotation_deg),
        center_x=center[0] * 1e-3,
        center_y=center[1] * 1e-3,
    )
    return UnitCell(lattice=lat, geometry=geo, substrate=sub, metal=Metal())


# ------------------------------------------------------------------------------- filenames


def test_step_filename_format():
    assert step_filename(1.0, 1.0) == "UC_L1.0_W1.0.step"
    assert step_filename(4.4, 4.4) == "UC_L4.4_W4.4.step"
    assert step_filename(2.6, 3.4) == "UC_L2.6_W3.4.step"


# --------------------------------------------------------------------------- bounds checking


def test_check_geometry_bounds_ok_for_a_centered_patch():
    uc = _cell(2.0, 2.0)
    bounds = check_geometry_bounds(uc)
    assert bounds.ok
    assert bounds.l_in_range and bounds.w_in_range and bounds.within_cell
    assert bounds.patch_x_min_mm == pytest.approx(-1.0)
    assert bounds.patch_x_max_mm == pytest.approx(1.0)
    assert bounds.cell_x_min_mm == pytest.approx(-2.4)
    assert bounds.cell_x_max_mm == pytest.approx(2.4)


def test_check_geometry_bounds_rotated_patch_uses_rotated_bbox():
    """The bounds check's patch bbox must be the *rotated* one (L|cos|+W|sin|)/2, not the raw
    (unrotated) half-length -- checked against the closed form directly, not a vague inequality
    (for L=3, W=1 at 45 deg the rotated x half-extent, 1.414 mm, happens to be *smaller* than
    the unrotated 1.5 mm, which is itself informative: it shows this is really using the rotated
    formula, not merely "always bigger")."""
    uc = _cell(3.0, 1.0, rotation_deg=45.0)
    bounds = check_geometry_bounds(uc)
    assert bounds.ok
    c = s = math.sqrt(2.0) / 2.0
    expected_hx = (3.0 * c + 1.0 * s) / 2.0
    assert bounds.patch_x_max_mm == pytest.approx(expected_hx)
    assert bounds.patch_x_max_mm != pytest.approx(1.5)  # not the unrotated half-length


def test_check_geometry_bounds_flags_a_too_large_patch():
    """0 < L < period_x must fail for a patch as wide as the period itself (constructed
    independently of UnitCell's own guard, at the boundary)."""
    sub = Substrate(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    lat = RectangularLattice(period_x=4.8e-3, period_y=4.8e-3)
    geo = Geometry(kind="rectangle", length=4.8e-3 - 1e-9, width=2.0e-3)  # just inside, still valid
    uc = UnitCell(lattice=lat, geometry=geo, substrate=sub, metal=Metal())
    bounds = check_geometry_bounds(uc)
    assert bounds.l_in_range  # strictly < period_x holds numerically


# -------------------------------------------------------------------------- full export/validation


@pytest.fixture(scope="module")
def small_export(tmp_path_factory):
    cfg = LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=12.0e9,
        f_max_hz=13.0e9,
        f_step_hz=1.0e9,
        f_design_hz=12.5e9,
        l_values_mm=(1.0, 2.4, 4.2),
        w_values_mm=(1.0, 2.4, 4.2),
    )
    database = build_library_database(cfg, out_dir=None).database
    out = tmp_path_factory.mktemp("geom")
    results = export_library_geometries(database, cfg, out)
    return cfg, database, out, results


def test_export_writes_one_step_file_per_ok_cell(small_export):
    cfg, database, out, results = small_export
    n_ok = sum(1 for c in database.cells if c.status == "ok")
    assert len(results) == n_ok == 9
    for r in results.values():
        assert r.step_path.exists()
        assert r.step_path.name == step_filename(r.l_mm, r.w_mm)


def test_every_exported_solid_is_topologically_valid(small_export):
    _, _, _, results = small_export
    for r in results.values():
        imported = cadquery.importers.importStep(str(r.step_path))
        assert solid_validity(imported)
        assert len(imported.val().Solids()) == 2  # Substrate + Metal


def test_every_export_passes_its_own_bounds_check(small_export):
    _, _, _, results = small_export
    for r in results.values():
        assert r.bounds.ok


def test_step_geometry_matches_solved_unit_cell(small_export):
    """The exported STEP's metal footprint corresponds to the same (L, W) the solver used for
    that cell -- cross-checked via the independent closed-form bbox (cad_export, Phase 6)."""
    from src.cad_export import expected_feature_bbox_mm, import_and_summarize

    _, database, _, results = small_export
    for cell_id, r in results.items():
        uc = database.unit_cell([c for c in database.cells if c.cell_id == cell_id][0].index)
        expected = expected_feature_bbox_mm(uc.geometry)
        recovered = import_and_summarize(r.step_path)
        xmin, xmax, ymin, ymax, _, _ = recovered.metal.bbox_mm
        assert (xmin, xmax, ymin, ymax) == pytest.approx(expected, abs=1e-6)
