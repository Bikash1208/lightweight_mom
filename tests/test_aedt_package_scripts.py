"""Every generated AEDT script of the Phase 28 packages, run against a recording stub of the AEDT API (the 2026-10-27 AEDT run of phase28_loop_R4 failed on two things that this test pins down):
(1) AssignMaster / AssignSlave / AssignFloquetPort take exactly ONE face per boundary; the periodic walls are one master/slave pair per face, each pair with an origin in the plane of its face;
(2) no fallback may use a hard-coded design name or raise: every message goes through note() (oDesign.GetName()) and a failing boundary assignment or a failing AddMessage must not stop the script."""

import re
import sys
from pathlib import Path
from unittest import mock

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(p for d in ("phase28_cell_G", "phase28_coded_small", "phase28_loop_R4", "phase28_launcher") for p in (ROOT / "hfss_exports" / d).glob("*.py"))
pytestmark = pytest.mark.skipif(not SCRIPTS, reason="HFSS packages not generated")


def _desktop(fail_assign=False, fail_message=False):
    desk = mock.MagicMock()
    bnd = desk.NewProject.return_value.SetActiveDesign.return_value.GetModule.return_value  # every GetModule() returns this one mock
    if fail_assign:
        for name in ("AssignMaster", "AssignSlave", "AssignFloquetPort", "AssignLumpedRLC", "EditLumpedRLC", "DeleteBoundaries"):
            getattr(bnd, name).side_effect = Exception("Assignment must contain exactly one selection.")
    if fail_message:
        desk.AddMessage.side_effect = Exception("message failed")
    return desk, bnd


def _run(path, **kw):
    desk, bnd = _desktop(**kw)
    # set_state / set_coding work on an open model
    active = desk.GetActiveProject.return_value.GetActiveDesign.return_value
    active.GetModule.return_value = bnd
    glob = {"oDesktop": desk, "__name__": "__aedt__"}
    text = path.read_text(encoding="utf-8")
    if path.name == "set_coding.py":
        text = text.replace('CODING_CSV = r"codings/coding_beam_030_225.csv"', f'CODING_CSV = r"{(ROOT / "hfss_exports/phase28_coded_small/codings/coding_beam_030_225.csv").as_posix()}"')
    with mock.patch.dict(sys.modules, {"ScriptEnv": mock.MagicMock()}):
        exec(compile(text, str(path), "exec"), glob)  # noqa: S102 - generated, committed script, stub environment
    return desk, bnd, glob


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_hard_coded_design_name_in_messages_and_no_f_strings(path):
    text = path.read_text(encoding="utf-8")
    assert not re.search(r'oDesktop\.AddMessage\(oProject\.GetName\(\), "', text), "hard-coded design name in a message"
    assert not re.search(r"""\bf["']""", text) and "print(" not in text
    desk, _, _ = _run(path)
    for c in desk.AddMessage.call_args_list:
        assert not isinstance(c.args[1], str)  # oDesign.GetName() (a stub object here), never a literal


@pytest.mark.parametrize("path", [p for p in SCRIPTS if p.name.startswith("build")], ids=lambda p: f"{p.parent.name}/{p.name}")
def test_every_boundary_gets_exactly_one_face_and_periodic_pairs_match(path):
    _, bnd, _ = _run(path)
    masters = [c.args[0] for c in bnd.AssignMaster.call_args_list]
    slaves = [c.args[0] for c in bnd.AssignSlave.call_args_list]
    for a in masters + slaves + [c.args[0] for c in bnd.AssignFloquetPort.call_args_list]:
        faces = a[a.index("Faces:=") + 1]
        assert isinstance(faces, list) and len(faces) == 1, a[1]
    if "cell" in path.name:  # single-cell Floquet packages: 4 master/slave pairs (substrate X, air X, substrate Y, air Y), one Floquet port
        assert len(masters) == len(slaves) == 4 and bnd.AssignFloquetPort.call_count == 1
        assert sorted(a[0] for a in masters) == ["NAME:MasterXa", "NAME:MasterXs", "NAME:MasterYa", "NAME:MasterYs"]
        for s in slaves:
            assert s[s.index("Master:=") + 1] == "Master" + s[0].split("Slave")[1]
        for a in masters + slaves:  # origin and u point of a pair lie in the plane of its face (same z), one P long
            v = a[a.index("CoordSysVector:=") + 1]
            o, u = v[2], v[4]
            assert o[2] == u[2]
    else:
        assert not masters and not slaves


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_failing_assignments_and_failing_messages_never_stop_the_script(path):
    _run(path, fail_assign=True, fail_message=True)
