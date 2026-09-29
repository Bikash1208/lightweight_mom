"""Sec. 27-28: Phase 9 output writers round-trip correctly."""

import csv
import json

from src.phase9.control_variable import SolverReferenceRow
from src.phase9.gate import decide, evaluate_gate_a, evaluate_gate_b, GATE_A_CHECKLIST_KEYS, GATE_B_CHECKLIST_KEYS
from src.phase9.hfss_reference import HfssRecord
from src.phase9.local_periodicity import LpaResult
from src.phase9.specification import (
    write_hfss_reference_cells_csv,
    write_lpa_error_metrics_csv,
    write_phase9_gate_status_json,
    write_phase9_summary_json,
    write_s_parameter_comparison_csv,
    write_supercell_definitions_csv,
    write_supercell_results_csv,
)


def _row():
    return SolverReferenceRow("CELL-01", 1.0, 1.0, 12.5, 0.1 + 0.2j, 0.9 + 0.05j)


def test_write_hfss_reference_cells_csv_reports_not_available_without_hfss(tmp_path):
    path = write_hfss_reference_cells_csv([_row()], None, tmp_path / "hfss.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert "NOT_AVAILABLE" in rows[1]


def test_write_hfss_reference_cells_csv_includes_real_hfss_values_when_present(tmp_path):
    hfss = (HfssRecord("CELL-01", 1.0, 1.0, 12.5, 0.11 + 0.19j, 0.88 + 0.04j),)
    path = write_hfss_reference_cells_csv([_row()], hfss, tmp_path / "hfss.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert "NOT_AVAILABLE" not in rows[1]


def test_write_s_parameter_comparison_csv_not_available_without_hfss(tmp_path):
    path = write_s_parameter_comparison_csv([_row()], None, "s11", tmp_path / "s11.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert rows[1][-1] == "NOT_AVAILABLE"


def test_write_supercell_definitions_csv_row_count(tmp_path):
    grid = [[(2.8e-3, 2.8e-3)] * 3 for _ in range(3)]
    path = write_supercell_definitions_csv([("case1", 3, grid)], tmp_path / "defs.csv")
    with open(path) as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 1 + 9


def test_write_supercell_results_and_lpa_metrics_csv(tmp_path):
    result = LpaResult("3x3_weak", 3, 1 + 1j, 1.1 + 1.1j, 0.14, 0.1, 5.0, 0.1, 40, 300)
    p1 = write_supercell_results_csv([result], tmp_path / "results.csv")
    p2 = write_lpa_error_metrics_csv([result], tmp_path / "metrics.csv")
    with open(p1) as fh:
        assert len(list(csv.reader(fh))) == 2
    with open(p2) as fh:
        assert len(list(csv.reader(fh))) == 2


def test_write_phase9_summary_and_gate_status_json_round_trip(tmp_path):
    gate_a = evaluate_gate_a({k: True for k in GATE_A_CHECKLIST_KEYS})
    gate_b = evaluate_gate_b({k: True for k in GATE_B_CHECKLIST_KEYS})
    decision = decide(gate_a, gate_b, "UNVALIDATED", "UNVALIDATED", "INVALIDATED")
    p1 = write_phase9_summary_json(decision, tmp_path / "summary.json")
    p2 = write_phase9_gate_status_json(decision, tmp_path / "gate.json")
    summary = json.loads(p1.read_text())
    gate_status = json.loads(p2.read_text())
    assert summary["final_status"] == "LOCAL_PERIODICITY_INVALIDATED"
    assert gate_status["geometry_inversion"]["status"] == "LOCKED"
