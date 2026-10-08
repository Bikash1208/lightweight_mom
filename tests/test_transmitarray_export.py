"""Phase 19 Part D -- lens STEP export on a small lens: one valid solid, bbox, volume vs an independent raster."""

import pytest

pytest.importorskip("cadquery")

from src.transmitarray import export as ex  # noqa: E402
from src.transmitarray.layout import LensLayout  # noqa: E402

D, P, PLATE, TOTAL = 30.0, 5.0, 0.6, 25.0


def small_cells():
    lay = LensLayout(D * 1e-3, P * 1e-3, 36e-3)
    _ix, _iy, x, y = lay.cells()
    cells = []
    for i, (xx, yy) in enumerate(zip(x * 1e3, y * 1e3, strict=True)):
        if i % 3 == 0:
            cells.append(dict(x_mm=float(xx), y_mm=float(yy), kind="hole", radius_mm=1.6 + 0.05 * (i % 5)))
        elif i % 3 == 1:
            cells.append(dict(x_mm=float(xx), y_mm=float(yy), kind="pillar", radius_mm=1.2 + 0.1 * (i % 4)))
        else:
            cells.append(dict(x_mm=float(xx), y_mm=float(yy), kind="pillar", radius_mm=0.0))
    return cells


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    cells = small_cells()
    body = ex.build_lens(cells, D, P, PLATE, TOTAL)
    path = tmp_path_factory.mktemp("lens") / "lens.step"
    ex.save_step_named(body, path)
    return cells, path


def test_step_is_one_valid_solid_with_the_designed_bbox_and_name(exported):
    cells, path = exported
    v = ex.verify_step(path, cells, D, P, PLATE, TOTAL)
    assert v["n_solids"] == 1 and v["valid"]
    for got, want in zip(v["bbox_mm"], v["bbox_expected_mm"], strict=True):
        assert got == pytest.approx(want, abs=1e-6)
    assert v["name_in_step"]


def test_step_volume_matches_the_independent_raster(exported):
    cells, path = exported
    v = ex.verify_step(path, cells, D, P, PLATE, TOTAL)
    assert v["volume_rel_diff"] < 0.004


def test_pillar_top_faces_match_the_layout(exported):
    cells, path = exported
    v = ex.verify_step(path, cells, D, P, PLATE, TOTAL)
    assert v["pillar_top_faces_found"] == v["pillar_cells"] > 0


def test_designed_vertices_are_rounded_and_the_step_has_no_float_noise(exported):
    _cells, path = exported
    for px, py in ex.octagon_pts(3.0, -7.5, 1.77777777):
        assert abs(px * 1e4 - round(px * 1e4)) < 1e-6 and abs(py * 1e4 - round(py * 1e4)) < 1e-6
    import re

    text = path.read_text(errors="ignore")
    noisy = 0
    for m in re.finditer(r"CARTESIAN_POINT\('[^']*',\(([^)]*)\)\)", text):
        for tok in m.group(1).split(","):
            d = abs(float(tok) * 1e4 - round(float(tok) * 1e4)) * 1e-4
            noisy += 1e-13 < d < 1e-9  # 4.9999999999-style noise: almost, but not exactly, a rounded value
    assert noisy == 0


def test_stl_is_written(exported, tmp_path):
    cells, _path = exported
    body = ex.build_lens(cells, D, P, PLATE, TOTAL)
    ex.save_stl(body, tmp_path / "lens.stl")
    assert (tmp_path / "lens.stl").stat().st_size > 1000
