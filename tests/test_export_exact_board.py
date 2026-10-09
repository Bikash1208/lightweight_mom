"""Phase 21 D: the exact multi-port board export -- STEP from exact vertex coordinates (verified after
writing), layout CSV, the AEDT script and the README -- checked against the hand-built reference board."""

import ast
import csv
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
