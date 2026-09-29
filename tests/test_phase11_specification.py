"""Sec. 27, 28: the required JSON output writers, exercised against small synthetic objects."""

import json

from src.phase11 import specification as spec_mod
from src.phase11.convergence import ConvergenceResult
from src.phase11.classification import classify_overall
from src.phase11.coupling_radius import CouplingRadiusPoint, ShellIncrement
from src.phase11.interaction import PairwiseInteractionResult


def test_write_coupling_radius_results_json_round_trips(tmp_path):
    points = {"C3": [CouplingRadiusPoint(1, 1 + 2j, 24), CouplingRadiusPoint(3, 1.1 + 2.1j, 136)]}
    path = spec_mod.write_coupling_radius_results_json(points, tmp_path / "cr.json")
    doc = json.loads(path.read_text())
    assert doc["C3"]["points"][0]["response"] == [1.0, 2.0]


def test_write_incremental_shell_results_json(tmp_path):
    incs = {"C3": [ShellIncrement(1, 3, 0.5, 0.1)]}
    path = spec_mod.write_incremental_shell_results_json(incs, tmp_path / "shells.json")
    doc = json.loads(path.read_text())
    assert doc["C3"][0]["delta_relative"] == 0.1


def test_write_convergence_json(tmp_path):
    cr = [ConvergenceResult("C3", "CONVERGED", 1e-9, 1e-9, 3)]
    classification = classify_overall(cr, pattern_dependence_detected=False)
    path = spec_mod.write_convergence_json(cr, classification, tmp_path / "conv.json")
    doc = json.loads(path.read_text())
    assert doc["classification"]["overall_case"] == "CASE_A_LOCAL"
    assert doc["noise_floor_relative_declared_before_evaluation"] is True


def test_write_pairwise_interaction_json(tmp_path):
    results = [PairwiseInteractionResult("N+S", 1 + 0j, 1 + 0j, 2 + 0j, 2 + 0j, 0.0, 0.0)]
    path = spec_mod.write_pairwise_interaction_json(0j, results, tmp_path / "pw.json")
    doc = json.loads(path.read_text())
    assert doc["results"][0]["pair_label"] == "N+S"


def test_write_phase11_summary_and_gate_status_json_are_consistent(tmp_path):
    from src.phase11 import gate as gate_mod
    from src.phase11.classification import CASE_A_LOCAL
    from src.phase11.reduction import REDUCTION_VALIDATED

    geometry_conditions = {k: True for k in gate_mod.GEOMETRY_INVERSION_CONDITION_KEYS}
    step_conditions = {k: True for k in gate_mod.STEP_EXPORT_CONDITION_KEYS}
    decision = gate_mod.decide(CASE_A_LOCAL, True, True, REDUCTION_VALIDATED, geometry_conditions, step_conditions)
    summary_path = spec_mod.write_phase11_summary_json(decision, tmp_path / "summary.json")
    status_path = spec_mod.write_phase11_gate_status_json(decision, tmp_path / "status.json")
    summary = json.loads(summary_path.read_text())
    status = json.loads(status_path.read_text())
    assert summary["geometry_inversion_status"] == status["geometry_inversion_status"]
    assert summary["geometry_inversion_status"] == "GEOMETRY_INVERSION_UNLOCKED"
