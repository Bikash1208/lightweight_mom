"""Phase 28 loop cell R4 HFSS package: the AEDT scripts run against a recording stub of the AEDT API and create the expected objects and single-type lumped boundaries (R only for ON, C only for OFF, an
inductor in series when L_s > 0); the rotated MoM geometry is the unrotated one turned by the rotation angle."""

import math
import re
import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

PKG = ROOT / "hfss_exports" / "phase28_loop_R4"
needs_pkg = pytest.mark.skipif(not (PKG / "build_cell.py").exists(), reason="loop package not generated")


def _run(text, **overrides):
    for k, v in overrides.items():
        text = re.sub(rf"^{k} = .*?(?=\s+#|\n)", f"{k} = {v!r}", text, count=1, flags=re.M)
    glob = {"oDesktop": mock.MagicMock(), "__name__": "__aedt__"}
    with mock.patch.dict(sys.modules, {"ScriptEnv": mock.MagicMock()}):
        exec(compile(text, "build_cell.py", "exec"), glob)  # noqa: S102 - generated, committed script, stub environment
    return glob


def _boundaries(glob):
    return [c.args[0] for c in glob["oBnd"].AssignLumpedRLC.call_args_list]


def _flag(args, key):
    return args[args.index(key) + 1]


@needs_pkg
def test_no_py3_only_syntax_in_the_aedt_scripts():
    for p in PKG.glob("*.py"):
        text = p.read_text(encoding="utf-8")
        assert not re.search(r"""\bf["']""", text), p.name
        assert "print(" not in text


@needs_pkg
@pytest.mark.parametrize("ls,state,rot", [(0.0, "ON", 0.0), (0.0, "OFF", 0.0), (0.2, "ON", 0.0), (0.2, "OFF", 20.0)])
def test_build_cell_boundaries(ls, state, rot):
    g = _run((PKG / "build_cell.py").read_text(encoding="utf-8"), LS_NH=ls, STATE=state, ROT_DEG=rot)
    bnd = _boundaries(g)
    assert len(bnd) == (2 if ls > 0 else 1)
    first = bnd[0]
    assert _flag(first, "UseResist:=") == (state == "ON") and _flag(first, "UseCap:=") == (state == "OFF") and _flag(first, "UseInduct:=") is False
    assert _flag(first, "Resistance:=") == ("2.4ohm" if state == "ON" else "0ohm") and _flag(first, "Capacitance:=") == ("0pF" if state == "ON" else "0.0282pF")
    if ls > 0:
        assert _flag(bnd[1], "UseInduct:=") is True and _flag(bnd[1], "Inductance:=") == "0.2nH" and not _flag(bnd[1], "UseResist:=") and not _flag(bnd[1], "UseCap:=")
    # five ring pieces (+ island when L_s > 0) and the diode sheets
    assert g["oEditor"].CreatePolyline.call_count == (5 if ls == 0 else 6) + len(bnd)


@needs_pkg
def test_current_line_follows_the_rotation():
    g = _run((PKG / "build_cell.py").read_text(encoding="utf-8"), LS_NH=0.0, STATE="ON", ROT_DEG=20.0)
    a = _boundaries(g)[0]
    line = a[a.index("NAME:CurrentLine") - 0 + 0] if False else [x for x in a if isinstance(x, list) and x and x[0] == "NAME:CurrentLine"][0]
    s, e = line[2], line[4]
    num = lambda t: float(t[:-2])  # noqa: E731
    ang = math.degrees(math.atan2(num(e[1]) - num(s[1]), num(e[0]) - num(s[0])))
    assert ang == pytest.approx(20.0, abs=1e-3)


@needs_pkg
def test_set_state_edits_the_first_element():
    text = (PKG / "set_state.py").read_text(encoding="utf-8")
    for state, res in (("ON", True), ("OFF", False)):
        g = _run(text, STATE=state)
        call = g["oBnd"].EditLumpedRLC.call_args
        assert call.args[0] == "D0_RLC" and _flag(call.args[1], "UseResist:=") is res


def test_rotated_loop_mask_is_the_unrotated_mask_turned():
    import phase28_loop_rotation as lr

    a = math.radians(20.0)
    x, y = np.meshgrid(np.linspace(-4.4e-3, 4.4e-3, 181), np.linspace(-4.4e-3, 4.4e-3, 181))
    rotated = lr.RotLoop(lr.SO, lr.WR, 20.0).in_feature(x, y)
    base = lr.RotLoop(lr.SO, lr.WR, 0.0).in_feature(x * math.cos(a) + y * math.sin(a), -x * math.sin(a) + y * math.cos(a))
    assert np.array_equal(rotated, base)
    assert lr.RotLoop(lr.SO, lr.WR, 20.0).corners().max() < lr.PITCH / 2
