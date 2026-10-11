"""Phase 21 D: the exact multi-port board export -- STEP from exact vertex coordinates (verified after
writing), layout CSV, the AEDT script and the README -- checked against the hand-built reference board."""

import ast
import csv
import re
from pathlib import Path

import pytest

from src.export.aedt_script import generate_aedt_script
from src.export.board_readme import build_board_readme
from src.export.exact_board import exact_board_from_ports_result, verify_step_exact, write_exact_step, write_layout_csv
from tests.test_app_phase21_multiport import _Job, _grid_request
from src.app.routers.hmsa import _run_synthesize_ports_job

REFERENCE = Path("hfss_exports/four_port_1bit_square_300mm_gap28")


@pytest.fixture(scope="module")
def small_result(tmp_path_factory):
    import os

    os.environ["LWMOM_APP_DATA_DIR"] = str(tmp_path_factory.mktemp("appdata"))
    return _run_synthesize_ports_job(_Job("p21export_small"), _grid_request(32.0, 28.0))


@pytest.fixture(scope="module")
def small_board(small_result):
    return exact_board_from_ports_result(small_result, 6.4, 1.575)


def test_cells_sit_at_each_ports_own_exact_positions(small_board):
    from decimal import Decimal as D

    # 2x2 grid of 32 mm ports (5 x 5 cells), 28 mm gap: port 1 spans x 14..46, so its first cell starts at 14.275
    port1 = [p for p in small_board.patches if p.port == 1]
    assert all(p.x1 - p.x0 == D("5.85") for p in port1)
    assert min(p.x0 for p in port1) == D("14.275")
    assert small_board.x_min == -small_board.x_max and small_board.width_mm == D(32 * 2 + 28)
    assert small_board.max_rounding_mm < 1e-9


def test_step_is_exact_and_the_verifier_catches_residue(small_board, tmp_path):
    step = write_exact_step(small_board, tmp_path / "board.step")
    assert verify_step_exact(step) > 100
    bad = tmp_path / "bad.step"
    bad.write_text("#1 = CARTESIAN_POINT('',(1.9999999999,0.0,2.0000001));\n")
    with pytest.raises(ValueError, match="not exact"):
        verify_step_exact(bad)


def test_aedt_script_is_valid_python_with_every_patch_and_feed(small_board, small_result):
    text = generate_aedt_script(small_board, "HMSA test", 2.2, 0.0009, 10.54e9)
    ast.parse(text)
    assert "@@" not in text
    assert text.count("(") > len(small_board.patches)
    patches_block = text.split("PATCHES = [")[1].split("\n]")[0]
    assert len([ln for ln in patches_block.splitlines() if ln.startswith("(")]) == len(small_board.patches)
    assert "HMSA_test" in text and '"Frequency:=", "10.54GHz"' in text and "mono_h" in text and "AssignLumpedPort" in text


def _zero_conductivity_solve_inside_violations(source):
    """attrs(name, material, solve_inside, ...) calls pairing vacuum/air with SolveInside False."""
    bad = []
    source = re.sub(r"@@\w+@@", "0", source)  # template placeholders
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "attrs" and len(node.args) >= 3:
            material, inside = node.args[1], node.args[2]
            if (isinstance(material, ast.Constant) and str(material.value).lower() in ("vacuum", "air")
                    and not (isinstance(inside, ast.Constant) and inside.value is True)):
                bad.append((node.lineno, material.value))
    return bad


def test_aedt_script_never_pairs_vacuum_or_air_with_solve_inside_false(small_board):
    text = generate_aedt_script(small_board, "HMSA test", 2.2, 0.0009, 10.54e9)
    assert text.count("attrs(") >= 5  # the check below really sees the calls
    assert _zero_conductivity_solve_inside_violations(text) == []
    template = Path(__file__).resolve().parents[1] / "scripts" / "adhoc" / "square4" / "aedt_template.py"
    if template.exists():
        assert _zero_conductivity_solve_inside_violations(template.read_text()) == []
    # the helper itself catches the 2025 R1 failure mode
    assert _zero_conductivity_solve_inside_violations('attrs("s", "vacuum", False)') == [(1, "vacuum")]


def test_readme_names_the_targeted_aedt_versions(small_board, small_result):
    text = build_board_readme(small_board, small_result, "x", 2.2, 0.0009, 10.54e9, ["hfss_build_board.py"])
    assert "2019 R1" in text and "2025 R1" in text and "SolveInside" in text


def test_readme_has_the_port_beam_and_cell_tables(small_board, small_result):
    text = build_board_readme(small_board, small_result, "x", 2.2, 0.0009, 10.54e9, ["board.step", "layout.csv"])
    for needle in ("## Ports and feeds", "## Beams", "## Cell states", "Edit Sources", "empty (bare substrate)"):
        assert needle in text


def test_annular_cells_are_refused_by_the_exact_path(small_result):
    import copy

    bad = copy.deepcopy(small_result)
    bad["ports"][0]["patches_mm"]["cells"][0]["geometry"] = {"kind": "annular_slot", "a_mm": 4.0}
    with pytest.raises(Exception, match="rectangular patches only"):
        exact_board_from_ports_result(bad, 6.4, 1.575)


@pytest.mark.skipif(not (REFERENCE / "layout.csv").is_file(), reason="reference export folder not present")
def test_layout_matches_the_hand_built_reference_board(tmp_path):
    """The app's own grid layout + coded 1-bit design reproduces the reference board's patch layout."""
    import os

    os.environ["LWMOM_APP_DATA_DIR"] = str(tmp_path / "appdata")
    result = _run_synthesize_ports_job(_Job("p21export_ref"), _grid_request(300.0, 28.0))
    board = exact_board_from_ports_result(result, 6.4, 1.575)
    path = write_layout_csv(board, tmp_path / "layout.csv")
    ours = {tuple(r) for r in list(csv.reader(open(path)))[1:]}
    ref = {tuple(r) for r in list(csv.reader(open(REFERENCE / "layout.csv")))[1:]}
    assert len(ours) == len(ref) == 4416
    assert ours == ref


# ======================================================================= Phase 22: feeds in the exact STEP and the package


@pytest.fixture(scope="module", params=[0.0, 28.0], ids=["gap0", "gap28"])
def full_board(request, tmp_path_factory):
    import os

    os.environ["LWMOM_APP_DATA_DIR"] = str(tmp_path_factory.mktemp("appdata_full"))
    gap = request.param
    result = _run_synthesize_ports_job(_Job(f"p22full{int(gap)}"), _grid_request(300.0, gap))
    return gap, result, exact_board_from_ports_result(result, 6.4, 1.575)


def _import_feed_faces(step_path):
    import cadquery as cq

    shape = cq.importers.importStep(str(step_path)).val()
    faces = [f for f in shape.Faces() if abs(f.Center().z) < 1e-9 and abs(f.Area() - 0.36) < 1e-6]
    return sorted((round(f.Center().x, 6), round(f.Center().y, 6)) for f in faces)


def test_step_has_feed_markers_at_the_exact_feed_coordinates(full_board, tmp_path):
    gap, _result, board = full_board
    step = write_exact_step(board, tmp_path / "board.step")
    text = step.read_text()
    for k in range(1, 5):
        assert f"PRODUCT('Feed{k}'" in text
    assert "Monopole" not in text
    c = 150.4 + gap / 2.0
    assert [(float(x), float(y)) for x, y in board.feeds_mm] == [(c, c), (-c, c), (-c, -c), (c, -c)]
    expected = sorted([(c, c), (-c, c), (-c, -c), (c, -c)])
    assert _import_feed_faces(step) == [(round(x, 6), round(y, 6)) for x, y in expected]


def test_step_verifies_with_and_without_pins_and_pins_are_exact(full_board, tmp_path):
    from src.export.exact_board import DEFAULT_MONO_H_MM, DEFAULT_R_PIN_MM

    gap, _result, board = full_board
    plain = write_exact_step(board, tmp_path / "plain.step")
    assert verify_step_exact(plain) > 100
    pinned = write_exact_step(board, tmp_path / "pins.step", include_pins=True)
    text = pinned.read_text()
    assert all(f"PRODUCT('Monopole{k}'" in text and f"PRODUCT('Feed{k}'" in text for k in range(1, 5))
    c = 150.4 + gap / 2.0
    pins = [(sx * c, sy * c, float(DEFAULT_R_PIN_MM), float(DEFAULT_MONO_H_MM)) for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1))]
    assert verify_step_exact(pinned, strict_points=False, pins=pins) > 100
    with pytest.raises(ValueError, match="monopole pin"):
        verify_step_exact(pinned, strict_points=False, pins=[(c, c, 0.635, 7.2)])
    # the pin geometry is on the grid: pin vertices at z = 0 and exactly z = 7.1, at feed +- 0.635 in x
    from src.export.exact_board import _step_points

    points, vertices = _step_points(text)
    tops = [points[i] for i in vertices if i in points and len(points[i]) == 3 and points[i][2] == 7.1]
    assert sorted(round(v[0], 9) for v in tops if abs(v[1] - c) < 1e-9) == sorted([round(-c + 0.635, 9), round(c + 0.635, 9)])


def test_layout_files_list_the_feeds(full_board, tmp_path):
    from src.export.exact_board import write_feeds_csv, write_layout_dxf, write_layout_json

    gap, _result, board = full_board
    c = 150.4 + gap / 2.0
    feeds = write_feeds_csv(board, tmp_path / "feeds.csv").read_text().splitlines()
    assert feeds[0] == "port,x_mm,y_mm" and feeds[1:] == [f"1,{c:g},{c:g}", f"2,{-c:g},{c:g}", f"3,{-c:g},{-c:g}", f"4,{c:g},{-c:g}"]
    import json

    data = json.loads(write_layout_json(board, tmp_path / "layout.json").read_text())
    assert [f["port"] for f in data["feeds"]] == [1, 2, 3, 4] and data["feeds"][0]["x_mm"] == f"{c:g}"
    assert len(data["patches"]) == 4416
    dxf = write_layout_dxf(board, tmp_path / "layout.dxf").read_text().splitlines()
    circles = [i for i, t in enumerate(dxf) if t == "CIRCLE"]
    assert len(circles) == 4 and all(dxf[i + 2] == "FEEDS" for i in circles)
    assert float(dxf[circles[0] + 4]) == pytest.approx(c) and float(dxf[circles[0] + 6]) == pytest.approx(c)
    assert {"P1", "P2", "P3", "P4"} <= set(dxf)


def test_standard_multi_port_export_lists_the_feeds_too(tmp_path):
    from src.export.array_export import export_multi_port_array
    from src.materials import Metal, Substrate

    grid = [[(0.0, 0.0), (5.85e-3, 5.85e-3)], [(5.85e-3, 5.85e-3), (0.0, 0.0)]]
    result = export_multi_port_array(grid, 12.8e-3, 12.8e-3, Substrate("s", 2.2, 0.0009, 1.575e-3), Metal(), tmp_path,
                                     include_step=False, feeds_mm=[(-3.2, -3.2), (3.2, 3.2)])
    assert Path(result.feeds_csv_path).read_text().splitlines() == ["port,x_mm,y_mm", "1,-3.200,-3.200", "2,3.200,3.200"]
    import json

    assert json.loads(Path(result.layout_json_path).read_text())["feeds"][1] == {"port": 2, "x_mm": 3.2, "y_mm": 3.2}
    dxf = Path(result.layout_dxf_path).read_text().splitlines()
    assert dxf.count("FEEDS") == 4 and "P2" in dxf  # two CIRCLE + two TEXT entities, layer FEEDS
