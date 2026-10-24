"""Phase 20 R5-6 infrastructure: the RIS board exporter on a small 3 x 3 case (no real package, no predictions)."""

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
from phase20_r3_opt_cp2 import make_cell  # noqa: E402

from src.export import ris_board_export as rb  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402
from src.ris.supercell import build_supercell  # noqa: E402

GRID = ((0, 1, 2), (3, 0, 1), (2, 3, 0))


@pytest.fixture(scope="module")
def board():
    cell = make_cell(
        [4.601469081921106, 2.475026504986856, 1.3749256942062944], 15e-3
    )  # the retained baseline CP-2 parameters
    return rb.make_board(cell, madp_kit(), GRID, name="unit_3x3")


def _inside(poly, x, y):
    inside = False
    n = len(poly)
    for k in range(n):
        (x0, y0), (x1, y1) = poly[k], poly[(k + 1) % n]
        if (y0 > y) != (y1 > y) and x < (x1 - x0) * (y - y0) / (y1 - y0) + x0:
            inside = not inside
    return inside


def test_pieces_plus_diode_gaps_are_exactly_the_cell_metal(board):
    pieces, sites = rb.cell_pieces_mm(board.cell)
    geo = board.cell.geometry()
    rng = np.random.default_rng(0)
    gap_polys = [rb._rect(x, y, span, w, np.radians(a)) for x, y, span, w, a in sites]
    for x, y in rng.uniform(-7.5, 7.5, (4000, 2)):
        in_pieces = any(_inside(p, x, y) for _, p in pieces)
        in_gap = any(_inside(g, x, y) for g in gap_polys)
        assert bool(geo.in_feature(x * 1e-3, y * 1e-3)) == (in_pieces or in_gap), (x, y)
        assert not (in_pieces and in_gap)  # the gap is a gap: no metal piece covers it


def test_diode_table_matches_the_supercell_loads_in_order_and_state(board):
    grid_states = [[board.cell.states()[board.grid[j][i]] for i in range(board.nx)] for j in range(board.ny)]
    loads = build_supercell(board.cell, grid_states, board.kit).loads
    d = rb.diodes(board)
    assert len(d) == len(loads) == 36
    for ld, dd in zip(loads, d, strict=True):
        assert dd.x_mm == pytest.approx(ld.x * 1e3, abs=1e-5) and dd.y_mm == pytest.approx(ld.y * 1e3, abs=1e-5)
        assert dd.span_mm == pytest.approx(ld.length * 1e3) and dd.width_mm == pytest.approx(ld.width * 1e3)
        assert dd.on == (
            ld.C is None and pytest.approx(5.2) == ld.R
        )  # ON: R only; OFF: C only (single-type boundaries)
        assert dd.on != (pytest.approx(board.kit.off.C_T) == ld.C)
    for j in range(3):
        for i in range(3):
            cell_d = [x for x in d if x.i == i and x.j == j]
            assert [x.stub for x in cell_d if x.on] == [GRID[j][i]]  # exactly one ON diode, the state's own stub
            assert [x.k for x in cell_d] == list(range(4 * (j * 3 + i) + 1, 4 * (j * 3 + i) + 5))


def test_layout_and_diodes_csv(board):
    rows = rb.layout_csv_text(board).splitlines()
    assert rows[0].startswith("row_j,col_i") and len(rows) == 10
    assert rows[1].split(",")[2:5] == [
        "-15",
        "-15",
        "1",
    ]  # cell (0, 0) of the 3 x 3 board is at (-15, -15) mm in state 1
    dl = rb.diodes_csv_text(board).splitlines()
    assert len(dl) == 37 and dl[1].startswith("Load_1,") and sum(",ON," in x for x in dl) == 9


def test_aedt_script_is_valid_python_with_every_piece_and_load(board):
    text = rb.aedt_script_text(board)
    assert "@@" not in text
    ast.parse(text)
    assert (
        text.count("'Load_") == 36
        and text.count("('C0_0_Oct'") == 1
        and "AssignLumpedRLC" in text
        and "NOT yet run in AEDT" in text
    )
    ns = {}
    code = text.split("oProject.GetDefinitionManager()")[0].split("oEditor.SetModelUnits")[1].split("\n", 1)[1]
    exec(compile(code, "aedt", "exec"), ns)  # the data tables evaluate without AEDT
    assert len(ns["PIECES"]) == 9 * 9 and len(ns["LOADS"]) == 36 and sum(1 for ld in ns["LOADS"] if ld[1]) == 9


def test_readme_has_no_invented_numbers_and_no_prediction_files(board, tmp_path):
    names = rb.export_board(tmp_path, board, with_step=False)
    assert set(names) == {"layout.csv", "diodes.csv", "hfss_build_board.py", "README.md"}
    assert not (tmp_path / "mom_reference.csv").exists() and not (tmp_path / "predictions.json").exists()
    text = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert "PREDICTIONS_NOT_YET_GENERATED" in text and "NOT_HFSS_VALIDATED" in text and "9 cells" in text
    assert (
        b"\r\n" not in (tmp_path / "hfss_build_board.py").read_bytes()
        and b"\r\n" not in (tmp_path / "README.md").read_bytes()
    )


def test_step_has_load_sheets_in_diode_order(board, tmp_path):
    pytest.importorskip("cadquery")
    res = rb.export_step(board, tmp_path / "board.step")
    assert len(res.load_bboxes_mm) == 36
    text = (tmp_path / "board.step").read_text(errors="ignore")
    for k in (1, 18, 36):
        assert f"Load_{k}'" in text or f"Load_{k}" in text
    assert "Load_37" not in text


def test_refusals():
    cell = make_cell([4.6, 2.4, 1.3], 15e-3)
    with pytest.raises(ValueError, match="refused"):
        rb.make_board(cell, madp_kit(), [[0] * 13] * 2)
    with pytest.raises(ValueError, match="0..3"):
        rb.make_board(cell, madp_kit(), [[0, 4]])
    with pytest.raises(ValueError):
        rb.make_board(cell, madp_kit(), [[0, 1], [0]])
