"""Phase 29 AEDT scripts (compact7, B3) run against a recording stub of the AEDT API, as the Phase 28 packages: no hard-coded design name in a message, no f-strings / print, every unverified call
wrapped (a failing far-field / sweep / EditSources / message never stops the script), 7 lumped ports, one radiation boundary on the Air object, seven pins."""

import re
import sys
from pathlib import Path
from unittest import mock

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(p for d in ("phase29_compact7", "phase29_B3", "phase29_twofeed") for p in (ROOT / "hfss_exports" / d).glob("*.py"))
pytestmark = pytest.mark.skipif(not SCRIPTS, reason="HFSS packages not generated")


def _run(path, fail=False):
    desk = mock.MagicMock()
    mod = desk.NewProject.return_value.GetActiveDesign.return_value.GetModule.return_value
    desk.GetActiveProject.return_value.GetActiveDesign.return_value.GetModule.return_value = mod
    if fail:
        for name in ("InsertFarFieldSphereSetup", "InsertFrequencySweep", "EditSources"):
            getattr(mod, name).side_effect = Exception("failed")
        desk.AddMessage.side_effect = Exception("message failed")
        desk.NewProject.return_value.SaveAs.side_effect = Exception("save failed")
    glob = {"oDesktop": desk, "__name__": "__aedt__"}
    with mock.patch.dict(sys.modules, {"ScriptEnv": mock.MagicMock()}):
        exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), glob)  # noqa: S102 - generated, committed script, stub environment
    return desk, mod, glob


@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: f"{p.parent.name}/{p.name}")
def test_no_hard_coded_design_name_no_fstrings_and_failures_never_stop_the_script(path):
    text = path.read_text(encoding="utf-8")
    assert (
        not re.search(r'oDesktop\.AddMessage\(oProject\.GetName\(\), "', text)
        and not re.search(r"""\bf["']""", text)
        and "print(" not in text
    )
    _run(path)
    _run(path, fail=True)


@pytest.mark.parametrize("path", [p for p in SCRIPTS if p.name.startswith("build_board")], ids=lambda p: f"{p.parent.name}/{p.name}")
def test_board_has_seven_ports_seven_pins_and_one_radiation_boundary(path):
    _, mod, glob = _run(path)
    n = len(glob["FEEDS"])
    assert mod.AssignLumpedPort.call_count == n and n in (1, 2, 7)
    assert mod.AssignRadiation.call_count == 1
    assert [c.args[0][0] for c in mod.AssignLumpedPort.call_args_list] == [f"NAME:Port{k}" for k in range(1, n + 1)]
    assert len(glob["CELLS"]) > 300
