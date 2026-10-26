"""Phase 20 R5-4: the supercell job script -- task lists and a one-cell consistency check (a 1 x 1 'supercell' must reproduce the isolated-cell prediction)."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import phase20_r5_supercell as sc  # noqa: E402

from src.ris.calibrated_library import CalibratedLibrary  # noqa: E402


def test_realisable_gradients_are_the_documented_set():
    got = {(q, n) for q, n, _ in sc.realisable_gradients()}
    assert (1, 4) in got and (2, 5) in got and (3, 8) in got and (1, 2) not in got  # 71.8 deg is beyond the 50 deg design range
    ang = {(q, n): a for q, n, a in sc.realisable_gradients()}
    assert ang[(1, 4)] == pytest.approx(28.42, abs=0.01) and ang[(2, 5)] == pytest.approx(49.59, abs=0.01)


def test_task_lists_cover_the_planned_cases():
    cal = CalibratedLibrary.load()
    g, r = sc.gradient_tasks(cal), sc.random_tasks(cal)
    assert len(g) == 11 and sum(t["ny"] > 1 for t in g) == 2
    assert sum(len(t["codings"]) for t in r) == 5 + 7 * 4
    assert all(len(c["grid"]) == 4 and len(c["grid"][0]) == 4 for t in r for c in t["codings"])
    assert all(0 <= v <= 3 for t in r for c in t["codings"] for row in c["grid"] for v in row)


def test_one_cell_supercell_reproduces_the_isolated_prediction(tmp_path):
    task = {"name": "unit", "nx": 1, "ny": 1, "codings": [{"label": "state 2", "grid": [[2]]}]}
    name, msg = sc.run_task((task, str(tmp_path), 8))
    rec = json.loads((tmp_path / "unit.json").read_text())
    o = rec["results"][0]["orders"][0]
    assert o["order"] == [0, 0] and abs(o["err_db"]) < 1e-6 and abs(o["err_deg"]) < 1e-4
    assert rec["floquet_factor"] == 2.0 and rec["provenance"]["git_commit"]
    assert sc.run_task((task, str(tmp_path), 8))[1].startswith("skipped")
