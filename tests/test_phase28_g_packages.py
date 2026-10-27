"""Phase 28 HFSS-first path: the G cell polygons reproduce the MoM geometry; the generated AEDT scripts run against a recording stub of the AEDT API and make the expected
number of objects and boundaries (ON = R only, OFF = C only, two ON diodes per cell and state); the codings are consistent with the layout."""

import csv
import re
import sys
from pathlib import Path
from unittest import mock

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_g_geometry as geo  # noqa: E402

BOARD = ROOT / "hfss_exports" / "phase28_coded_small"
CELL = ROOT / "hfss_exports" / "phase28_cell_G"
needs_pkg = pytest.mark.skipif(not (BOARD / "build_board.py").exists(), reason="HFSS packages not generated")


def test_polygons_and_necks_reproduce_the_mom_mask():
    n_bad, n_total, _ = geo.check_against_mask(401)
    assert n_bad == 0 and n_total > 100000


def test_states_turn_the_on_pair_by_45_degrees():
    on = geo.state_on_sites()
    assert [len(s) for s in on] == [2, 2, 2, 2]
    assert sorted(sum(on, [])) == list(range(8))  # every site is ON in exactly one state
    dv = geo.diode_values()
    assert dv == {"r_on_ohm": 5.2, "c_off_pf": 0.025, "l_s_nh": 0.0}


def _run(path, extra=None, cwd=None):
    """Execute an AEDT script against a recording stub; returns the recorded Assign/Edit calls."""
    stub = mock.MagicMock()
    glob = {"oDesktop": stub, "__name__": "__aedt__"}
    glob.update(extra or {})
    with mock.patch.dict(sys.modules, {"ScriptEnv": mock.MagicMock()}):
        code = compile(Path(path).read_text(encoding="utf-8"), str(path), "exec")
        exec(code, glob)  # noqa: S102 - generated, committed script, stub environment
    return glob


def _rlc_calls(glob):
    bnd = glob["oBnd"]
    return [c.args[0] for c in bnd.AssignLumpedRLC.call_args_list]


@needs_pkg
def test_no_f_strings_or_py3_only_syntax_in_the_aedt_scripts():
    for p in list(BOARD.glob("*.py")) + list(CELL.glob("*.py")):
        text = p.read_text(encoding="utf-8")
        assert not re.search(r"""\bf["']""", text), p.name
        assert ":=" in text and "print(" not in text


@needs_pkg
def test_board_script_builds_the_expected_objects_and_rlc_boundaries():
    g = _run(BOARD / "build_board.py")
    cells = g["CELLS"]
    calls = _rlc_calls(g)
    assert len(calls) == 8 * len(cells) <= 1200
    n_on = sum(1 for a in calls if a[a.index("UseResist:=") + 1])
    assert n_on == 2 * len(cells)
    for a in calls:
        on = a[a.index("UseResist:=") + 1]
        cap = a[a.index("UseCap:=") + 1]
        assert on != cap  # one element type each
        assert a[a.index("Resistance:=") + 1] == ("5.2ohm" if on else "0ohm")
        assert a[a.index("Capacitance:=") + 1] == ("0pF" if on else "0.025pF")
        assert a[a.index("UseInduct:=") + 1] is False
    names = [a[0] for a in calls]
    assert len(set(names)) == len(names)
    # 9 metal pieces + 8 diode sheets per cell
    assert g["oEditor"].CreatePolyline.call_count == 17 * len(cells)


@needs_pkg
@pytest.mark.parametrize("coding", sorted((BOARD / "codings").glob("coding_*.csv")))
def test_set_coding_applies_two_on_diodes_per_cell(coding):
    text = (BOARD / "set_coding.py").read_text(encoding="utf-8")
    rel = coding.as_posix()
    text = text.replace('CODING_CSV = r"codings/coding_beam_030_225.csv"', f'CODING_CSV = r"{rel}"')
    tmp = ROOT / "tests" / "_tmp_set_coding_run.py"
    tmp.write_text(text, encoding="utf-8")
    try:
        g = _run(tmp, {})
    finally:
        tmp.unlink()
    bnd = g["oBnd"]
    edits = [c.args for c in bnd.EditLumpedRLC.call_args_list]
    n_cells = len(g["CELLS"])
    assert len(edits) == 8 * n_cells
    assert sum(1 for _, a in edits if a[a.index("UseResist:=") + 1]) == 2 * n_cells
    rows = list(csv.DictReader(coding.open()))
    assert [int(r["cell_id"]) for r in rows] == [c[0] for c in g["CELLS"]]
    assert {int(r["state_0to3"]) for r in rows} <= {0, 1, 2, 3}


@needs_pkg
def test_codings_share_the_layout_and_the_three_beams_differ():
    lay = list(csv.DictReader((BOARD / "layout.csv").open()))
    codes = {}
    for p in sorted((BOARD / "codings").glob("coding_*.csv")):
        rows = list(csv.DictReader(p.open()))
        assert [(r["cell_id"], r["x_mm"], r["y_mm"]) for r in rows] == [(r["cell_id"], r["x_mm"], r["y_mm"]) for r in lay]
        assert all(int(r["axis_deg"]) == 45 * int(r["state_0to3"]) for r in rows)
        codes[p.name] = [r["state_0to3"] for r in rows]
    assert len(codes) == 3 and len({tuple(v) for v in codes.values()}) == 3


@needs_pkg
def test_cell_script_builds_eight_boundaries_with_two_on():
    g = _run(CELL / "build_cell.py")
    calls = _rlc_calls(g)
    assert len(calls) == 8 and sum(1 for a in calls if a[a.index("UseResist:=") + 1]) == 2
    assert g["oEditor"].CreatePolyline.call_count == 9 + 8
