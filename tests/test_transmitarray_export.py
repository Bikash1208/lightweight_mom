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


# ---------------------------------------------------------------------------------------------------------------------
# Part D-bis: exact coordinates in the written STEP
# ---------------------------------------------------------------------------------------------------------------------
import math  # noqa: E402
import re  # noqa: E402

import numpy as np  # noqa: E402


def cells_for(d, p):
    lay = LensLayout(d * 1e-3, p * 1e-3, 1.2 * d * 1e-3)
    _ix, _iy, x, y = lay.cells()
    out = []
    for i, (xx, yy) in enumerate(zip(x * 1e3, y * 1e3, strict=True)):
        k = i % 3
        out.append(
            dict(
                x_mm=float(xx),
                y_mm=float(yy),
                kind="hole" if k == 0 else "pillar",
                radius_mm=(1.5 + 0.07 * (i % 5)) if k < 2 else 0.0,
            )
        )
    return out


@pytest.fixture(scope="module", params=[(30.0, 5.0), (24.0, 4.0)], ids=["pitch5", "pitch4"])
def exact_export(request, tmp_path_factory):
    d, p = request.param
    cells = cells_for(d, p)
    body = ex.build_lens(cells, d, p, PLATE, TOTAL)
    path = tmp_path_factory.mktemp("exact") / "lens.step"
    report, intended = ex.export_exact_step(body, path, cells, d, p, PLATE, TOTAL)
    return dict(d=d, p=p, cells=cells, path=path, report=report, intended=intended, text=path.read_text())


def _coords(text, kinds=("CARTESIAN_POINT",)):
    rx = ex._regexes()
    out = []
    for name in ("point", "direction"):
        if name == "point" and "CARTESIAN_POINT" not in kinds:
            continue
        if name == "direction" and "DIRECTION" not in kinds:
            continue
        for m in rx[name].finditer(text):
            out += [t.strip() for t in m.group(2).split(",")]
    return out


def test_every_exact_dimension_parses_to_exactly_that_double(exact_export):
    e = exact_export
    toks = _coords(e["text"])
    doubles = {float(t) for t in toks}
    exact_dims = {
        "0": 0.0,
        "plate": PLATE,
        "total": TOTAL,
        "half aperture": e["d"] / 2,
        "pitch": e["p"],
        "half pitch": e["p"] / 2,
    }
    for name, want in exact_dims.items():
        assert want in doubles, f"{name} = {want!r} never appears as a coordinate"  # compares doubles
        assert ex.step_real(want) in toks, f"{name}: literal text {ex.step_real(want)!r} not written"
    # and nothing near a design value is the wrong double (the contract, for every intended coordinate)
    assert ex.inexact_literals(e["text"], e["intended"]) == []
    assert e["report"]["inexact_literals_after"] == 0


def test_the_kernel_text_really_was_inexact_so_this_test_has_power(exact_export):
    assert exact_export["report"]["inexact_literals_before"] > 0


def test_no_noise_literals_for_the_4mm_case(exact_export):
    if exact_export["p"] != 4.0:
        pytest.skip("4 mm pitch only")
    toks = _coords(exact_export["text"])
    assert "4." in toks and "2." in toks
    assert not any(re.fullmatch(r"3\.9999\d*|4\.0000\d*[1-9]\d*", t) for t in toks)


@pytest.mark.parametrize("target", ["0.6", "25.", "0."])
def test_one_ulp_error_is_detected(exact_export, target):
    e = exact_export
    pat = re.compile(r"(CARTESIAN_POINT\('',\([^)]*?)(?<![\d.E-])" + re.escape(target) + r"(?=[,)])")
    assert pat.search(e["text"]), f"{target} not found in a point"
    v = float(target)
    bumped = float(np.nextafter(v, np.inf))
    assert bumped != v
    broken = pat.sub(lambda m: m.group(1) + repr(bumped), e["text"], count=1)
    bad = ex.inexact_literals(broken, e["intended"])
    assert len(bad) >= 1 and bad[0][2] == v  # the predicate names the intended double the literal missed
    assert ex.inexact_literals(e["text"], e["intended"]) == []  # unperturbed file still passes


def test_a_real_feature_just_outside_the_snap_tolerance_is_not_moved():
    intended = dict(xy=np.array([0.0, 2.5]), z=np.array([0.0, 25.0]), dir=np.array([0.0, 1.0]))
    txt = "#1 = CARTESIAN_POINT('',(2.5000002,1.E-8,24.99999995));"
    new, stats = ex.snap_step_text(txt, intended, uncertainty_mm=1e-6)
    assert "2.5000002" in new  # 2e-7 off: a real feature, left alone
    assert "(2.5000002,0.,25.)" in new  # 1e-8 and 1e-7 are noise-sized: snapped to exact 0 and 25
    assert stats["n_snapped"] == 2


def test_step_real_is_shortest_round_trip():
    rng = np.random.default_rng(1)
    for v in [4.0, 0.085, 25.0, 1e-16, -2.5, 0.0, 123456.789, 0.1 + 0.2, 1 / 3, math.sqrt(0.5)]:
        s = ex.step_real(v)
        assert float(s) == v
        assert re.fullmatch(r"-?\d+\.\d*(E-?\d+)?", s), s  # valid STEP REAL (has a point, upper-case E)
    assert ex.step_real(4.0) == "4." and ex.step_real(-0.0) == "0." and ex.step_real(0.085) == "0.085"
    for v in rng.normal(size=500) * 10.0 ** rng.integers(-8, 3, size=500):
        assert float(ex.step_real(v)) == v
        assert len(ex.step_real(v)) <= len(repr(float(v))) + 2


def test_header_uncertainty_is_written_explicitly(exact_export):
    text = exact_export["text"]
    vals = re.findall(r"UNCERTAINTY_MEASURE_WITH_UNIT\(LENGTH_MEASURE\(([^)]*)\)", text)
    assert vals and all(float(v) == ex.UNCERTAINTY_MM for v in vals)
    assert ex.UNCERTAINTY_MM >= 10 * ex.SNAP_TOL_MM  # an order of magnitude above the largest permitted snap
    assert exact_export["report"]["vertex_check"]["min_distinct_vertex_separation_mm"] / 100 > ex.UNCERTAINTY_MM


def test_snap_keeps_the_solid_valid_and_the_volume_to_1e9(exact_export):
    r = exact_export["report"]
    assert r["no_distinct_vertices_merged"] and r["vertex_check"]["distinct_merge_violations"] == []
    assert r["vertex_check"]["collapsed_edges"] == 0
    assert r["n_solids"] == 1 and r["valid"] and r["self_intersection_free"]
    assert r["volume_rel_diff"] < 1e-9 and r["volume_within_1e9"]
    assert r["snap"]["max_snap_mm"] < ex.SNAP_TOL_MM


def test_self_intersection_analyzer_flags_an_overlapping_compound():
    import cadquery as cq

    one = cq.Workplane("XY").box(2, 2, 2).val()
    two = cq.Workplane("XY").box(2, 2, 2).translate((1, 0, 0)).val()
    assert ex.self_intersection_free(one)
    assert not ex.self_intersection_free(cq.Compound.makeCompound([one, two]))


def test_vertex_merge_check_catches_a_snap_that_collapses_distinct_vertices():
    raw = (
        "#1 = CARTESIAN_POINT('',(0.,0.,0.));\n#2 = CARTESIAN_POINT('',(0.00005,0.,0.));\n"
        "#3 = VERTEX_POINT('',#1);\n#4 = VERTEX_POINT('',#2);\n#5 = EDGE_CURVE('',#3,#4,#9,.T.);\n"
    )
    snapped = raw.replace("0.00005", "0.")  # a (bad) snap that moves a point 5e-5 mm
    rep = ex.verify_snap(raw, snapped)
    assert rep["distinct_merge_violations"] and rep["collapsed_edges"] == 1
