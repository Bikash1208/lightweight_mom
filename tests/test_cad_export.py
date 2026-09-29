"""Phase 6 geometry-consistency tests: STEP parameters == MoM geometry parameters.

Builds `UnitCell`s directly (the same objects `solver.py`/`mom_matrix.py` consume, Phase 1,
unmodified), exports them with `cad_export.export_unit_cell_step`, then independently re-imports
the file from disk (`cad_export.import_and_summarize`) and checks the recovered solids' bounding
boxes/volumes/centroids against closed-form expectations computed from the *same* `UnitCell`
object by a route that does not go through cadquery (`expected_feature_bbox_mm`, `unit_cell.
metal_area`, `units.to_mm`). Every case also proves "opens as a valid solid model"
(`cad_export.solid_validity`, OpenCASCADE's own `BRepCheck_Analyzer`) on both the pre-export and
the re-imported solids.

Needs the optional `cadquery` dependency (`requirements-cad.txt`); skipped entirely otherwise so
the base suite (`pip install -r requirements.txt` only) never depends on it.
"""

import math

import pytest

cadquery = pytest.importorskip("cadquery")

from src import cad_export, units
from src.geometry import Geometry, UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.validation_cases import CASES

_SUB = Substrate(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
_METAL_MAT = Metal()
_LAT = RectangularLattice(period_x=8.0e-3, period_y=8.0e-3)


def _cell(geometry, lattice=_LAT, substrate=_SUB):
    return UnitCell(lattice=lattice, geometry=geometry, substrate=substrate, metal=_METAL_MAT)


# A handful of ad hoc cells covering both polarities, rotation and an off-centre feature, plus
# every one of the 10 Phase 6 validation cases (built from the actual project config, not
# hand-rolled) -- so the consistency check exercises both "any UnitCell" and "the shipped cases".
_AD_HOC = {
    "metal_unrotated": _cell(Geometry(kind="rectangle", length=5.0e-3, width=2.0e-3)),
    "metal_rotated_30": _cell(
        Geometry(kind="dipole", length=6.0e-3, width=1.2e-3, rotation=math.radians(30.0))
    ),
    "metal_rotated_offset": _cell(
        Geometry(
            kind="rectangle", length=4.0e-3, width=3.0e-3,
            center_x=1.0e-3, center_y=-0.7e-3, rotation=math.radians(50.0),
        )
    ),
    "metal_rotation_90": _cell(Geometry(kind="rectangle", length=5.0e-3, width=2.0e-3, rotation=math.radians(90.0))),
    "slot_unrotated": _cell(Geometry(kind="slot", length=5.0e-3, width=2.0e-3)),
    "slot_rotated": _cell(Geometry(kind="slot", length=4.0e-3, width=1.5e-3, rotation=math.radians(20.0))),
}

_CASE_CELLS = {c.case_id: c.config.unit_cell for c in CASES}

ALL_CELLS = {**_AD_HOC, **{f"validation:{k}": v for k, v in _CASE_CELLS.items()}}


@pytest.fixture(params=list(ALL_CELLS), ids=list(ALL_CELLS))
def named_cell(request):
    return request.param, ALL_CELLS[request.param]


# ---------------------------------------------------------------------- pre-export solid checks

def test_build_solids_valid_and_countable(named_cell):
    """Before any export: exactly one substrate solid, one metal solid, both topologically valid."""
    _, uc = named_cell
    substrate, metal = cad_export.build_unit_cell_solids(uc)
    assert len(substrate.val().Solids()) == 1
    assert len(metal.val().Solids()) == 1
    assert cad_export.solid_validity(substrate)
    assert cad_export.solid_validity(metal)


# -------------------------------------------------------------------- export/import round trip

def test_step_round_trip_matches_ground_truth(named_cell, tmp_path):
    """export_unit_cell_step's own returned properties equal what import_and_summarize recovers
    independently from the STEP bytes on disk -- the export did not silently lose or distort
    anything between the in-memory solid and the written file."""
    name, uc = named_cell
    path = tmp_path / f"{name.replace(':', '_')}.step"
    exported = cad_export.export_unit_cell_step(uc, path)
    assert path.exists() and path.stat().st_size > 0

    recovered = cad_export.import_and_summarize(path)

    assert recovered.substrate.n_solids == exported.substrate.n_solids == 1
    assert recovered.metal.n_solids == exported.metal.n_solids == 1
    for a, b in zip(exported.substrate.bbox_mm, recovered.substrate.bbox_mm):
        assert a == pytest.approx(b, abs=1e-9)
    for a, b in zip(exported.metal.bbox_mm, recovered.metal.bbox_mm):
        assert a == pytest.approx(b, abs=1e-9)
    assert exported.substrate.volume_mm3 == pytest.approx(recovered.substrate.volume_mm3, rel=1e-9)
    assert exported.metal.volume_mm3 == pytest.approx(recovered.metal.volume_mm3, rel=1e-9)


def test_step_solids_are_valid_after_reimport(named_cell, tmp_path):
    """The 'every exported STEP opens as a valid solid model' acceptance criterion, checked with
    OpenCASCADE's own topology validator on the solids as read back from disk."""
    name, uc = named_cell
    path = tmp_path / f"{name.replace(':', '_')}.step"
    cad_export.export_unit_cell_step(uc, path)
    imported = cadquery.importers.importStep(str(path))
    assert cad_export.solid_validity(imported)


# --------------------------------------------------------------- STEP geometry == MoM geometry

def test_substrate_bbox_matches_lattice_and_thickness(named_cell, tmp_path):
    """Substrate solid recovered from the STEP file is exactly the Px x Py x h box the solver's
    own stack description implies (GREEN_FUNCTION.md): z in [-h, 0], centred at the cell origin."""
    name, uc = named_cell
    path = tmp_path / f"{name.replace(':', '_')}.step"
    cad_export.export_unit_cell_step(uc, path)
    recovered = cad_export.import_and_summarize(path)

    px_mm, py_mm = units.to_mm(uc.lattice.period_x), units.to_mm(uc.lattice.period_y)
    h_mm = units.to_mm(uc.substrate.thickness)
    xmin, xmax, ymin, ymax, zmin, zmax = recovered.substrate.bbox_mm
    assert (xmin, xmax) == pytest.approx((-px_mm / 2, px_mm / 2), abs=1e-9)
    assert (ymin, ymax) == pytest.approx((-py_mm / 2, py_mm / 2), abs=1e-9)
    assert (zmin, zmax) == pytest.approx((-h_mm, 0.0), abs=1e-9)


def test_metal_bbox_matches_closed_form_rotation(named_cell, tmp_path):
    """Metal bounding box recovered from the STEP file matches the closed-form rotated-rectangle
    bbox computed independently (not via cadquery) from the same Geometry object -- for metal
    (non-aperture) features, whose bbox is exactly the feature's own bbox."""
    name, uc = named_cell
    if uc.geometry.polarity != "metal":
        pytest.skip("aperture (slot) metal bbox is the full period plate, checked separately")
    path = tmp_path / f"{name.replace(':', '_')}.step"
    cad_export.export_unit_cell_step(uc, path)
    recovered = cad_export.import_and_summarize(path)

    expected = cad_export.expected_feature_bbox_mm(uc.geometry)
    xmin, xmax, ymin, ymax, _, _ = recovered.metal.bbox_mm
    assert (xmin, xmax, ymin, ymax) == pytest.approx(expected, abs=1e-9)


def test_metal_volume_matches_unit_cell_metal_area(named_cell, tmp_path):
    """metal solid volume / export thickness == UnitCell.metal_area (mm^2), for both metal and
    aperture (slot) polarity -- the same complement rule `UnitCell.metal_area`/`is_metal` use."""
    name, uc = named_cell
    path = tmp_path / f"{name.replace(':', '_')}.step"
    exported = cad_export.export_unit_cell_step(uc, path)
    recovered = cad_export.import_and_summarize(path)

    t_mm = units.to_mm(cad_export.DEFAULT_METAL_EXPORT_THICKNESS_M)
    # metal_area is an area in m^2; to_mm only scales lengths, so convert to mm^2 directly (1 m^2 = 1e6 mm^2).
    expected_area_mm2 = uc.metal_area * 1.0e6
    assert exported.metal.volume_mm3 / t_mm == pytest.approx(expected_area_mm2, rel=1e-6)
    assert recovered.metal.volume_mm3 / t_mm == pytest.approx(expected_area_mm2, rel=1e-6)


def test_metal_centroid_matches_feature_center(named_cell, tmp_path):
    """Recovered metal bounding-box centre equals (center_x, center_y) for metal-polarity
    features (a rectangle rotated about its own centroid has a bbox centred on that centroid,
    independent of the rotation angle); for slot/aperture polarity the metal is the full-period
    plate minus the feature, so its bbox centre is the cell centre instead."""
    name, uc = named_cell
    path = tmp_path / f"{name.replace(':', '_')}.step"
    cad_export.export_unit_cell_step(uc, path)
    recovered = cad_export.import_and_summarize(path)

    if uc.geometry.polarity == "metal":
        expected = (units.to_mm(uc.geometry.center_x), units.to_mm(uc.geometry.center_y))
    else:
        expected = (0.0, 0.0)
    assert recovered.metal.centroid_xy_mm == pytest.approx(expected, abs=1e-9)


# ------------------------------------------------------------------------------ sensitivity

def test_rotation_changes_recovered_bbox():
    """Non-vacuousness check: two otherwise-identical cells that differ only in rotation must
    give measurably different recovered metal bounding boxes -- the consistency checks above are
    actually sensitive to the geometry, not simply passing regardless of input."""
    base = _cell(Geometry(kind="rectangle", length=6.0e-3, width=2.0e-3, rotation=0.0))
    rotated = _cell(Geometry(kind="rectangle", length=6.0e-3, width=2.0e-3, rotation=math.radians(35.0)))

    bbox_base = cad_export.expected_feature_bbox_mm(base.geometry)
    bbox_rotated = cad_export.expected_feature_bbox_mm(rotated.geometry)
    assert bbox_base != pytest.approx(bbox_rotated, abs=1e-6)

    # And the STEP-recovered bboxes actually differ too, not just the closed-form expectation.
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p0, p1 = Path(d) / "base.step", Path(d) / "rotated.step"
        cad_export.export_unit_cell_step(base, p0)
        cad_export.export_unit_cell_step(rotated, p1)
        r0 = cad_export.import_and_summarize(p0)
        r1 = cad_export.import_and_summarize(p1)
    assert r0.metal.bbox_mm != pytest.approx(r1.metal.bbox_mm, abs=1e-6)


def test_offset_changes_recovered_centroid():
    """Same idea for center_x/center_y: an offset feature must recover a shifted centroid."""
    centered = _cell(Geometry(kind="rectangle", length=3.0e-3, width=3.0e-3))
    offset = _cell(Geometry(kind="rectangle", length=3.0e-3, width=3.0e-3, center_x=1.5e-3, center_y=-1.0e-3))

    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p0, p1 = Path(d) / "centered.step", Path(d) / "offset.step"
        cad_export.export_unit_cell_step(centered, p0)
        cad_export.export_unit_cell_step(offset, p1)
        r0 = cad_export.import_and_summarize(p0)
        r1 = cad_export.import_and_summarize(p1)
    assert r0.metal.centroid_xy_mm == pytest.approx((0.0, 0.0), abs=1e-9)
    assert r1.metal.centroid_xy_mm == pytest.approx((1.5, -1.0), abs=1e-9)
    assert r0.metal.centroid_xy_mm != pytest.approx(r1.metal.centroid_xy_mm, abs=1e-6)


def test_thickness_changes_metal_volume_but_not_bbox_footprint():
    """metal_export_thickness_m is a CAD-only convenience (module docstring): changing it must
    change the metal solid's volume/z-extent but leave its x/y footprint untouched."""
    uc = _cell(Geometry(kind="rectangle", length=4.0e-3, width=2.0e-3))
    thin, thick = 10e-6, 100e-6

    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p_thin, p_thick = Path(d) / "thin.step", Path(d) / "thick.step"
        r_thin = cad_export.export_unit_cell_step(uc, p_thin, metal_export_thickness_m=thin)
        r_thick = cad_export.export_unit_cell_step(uc, p_thick, metal_export_thickness_m=thick)

    xmin_t, xmax_t, ymin_t, ymax_t, _, _ = r_thin.metal.bbox_mm
    xmin_k, xmax_k, ymin_k, ymax_k, _, _ = r_thick.metal.bbox_mm
    assert (xmin_t, xmax_t, ymin_t, ymax_t) == pytest.approx((xmin_k, xmax_k, ymin_k, ymax_k), abs=1e-9)
    assert r_thick.metal.volume_mm3 == pytest.approx(r_thin.metal.volume_mm3 * (thick / thin), rel=1e-6)


# -------------------------------------------------------------------------------- error handling

def test_export_rejects_non_positive_thickness(tmp_path):
    uc = _cell(Geometry(kind="rectangle", length=4.0e-3, width=2.0e-3))
    with pytest.raises(ValueError):
        cad_export.export_unit_cell_step(uc, tmp_path / "bad.step", metal_export_thickness_m=-1.0)
    with pytest.raises(ValueError):
        cad_export.export_unit_cell_step(uc, tmp_path / "bad.step", metal_export_thickness_m=0.0)


def test_import_and_summarize_rejects_missing_file(tmp_path):
    with pytest.raises(Exception):
        cad_export.import_and_summarize(tmp_path / "does_not_exist.step")


def test_require_cadquery_error_names_requirements_file(monkeypatch):
    """The ImportError guard (exercised here without actually uninstalling cadquery) names the
    optional requirements file, so a user without it gets an actionable message."""
    monkeypatch.setattr(cad_export, "cq", None)
    monkeypatch.setattr(cad_export, "_IMPORT_ERROR", ImportError("no module named cadquery"))
    with pytest.raises(ImportError, match="requirements-cad.txt"):
        cad_export._require_cadquery()
