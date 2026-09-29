"""Phase 12A: the required JSON/CSV writers round-trip the data they are given."""

import csv
import json

from src.phase12a import gate, specification
from src.phase12a.analytic_checks import bare_grounded_slab_checks, vanishing_thickness_limit
from src.phase12a.sweep_study import coverage_by_frequency, run_grounded_patch_sweep

POINTS = run_grounded_patch_sweep(l_values_mm=(1.0, 2.0), frequencies_hz=(12.5e9,))
COVERAGE = coverage_by_frequency(POINTS)


def test_sweep_csv_has_one_row_per_point(tmp_path):
    path = specification.write_sweep_csv(POINTS, tmp_path / "sweep.csv")
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    assert len(rows) == 1 + len(POINTS)
    assert rows[0] == ["length_width_mm", "frequency_hz", "s11_xx_mag", "s11_xx_phase_deg", "s11_yx_mag"]


def test_sweep_json_round_trips_points_and_coverage(tmp_path):
    path = specification.write_sweep_json(POINTS, COVERAGE, tmp_path / "sweep.json")
    doc = json.loads(path.read_text())
    assert len(doc["points"]) == len(POINTS)
    assert str(12.5e9) in doc["coverage_by_frequency"]


def test_analytic_checks_json_includes_every_point(tmp_path):
    points = bare_grounded_slab_checks(frequencies_hz=(12.5e9,))
    thin = vanishing_thickness_limit()
    path = specification.write_analytic_checks_json(points, thin, tmp_path / "a.json")
    doc = json.loads(path.read_text())
    assert len(doc["bare_grounded_slab"]) == len(points)
    assert "vanishing_thickness_limit" in doc


def test_gate_status_json_round_trips():
    decision = gate.decide(1e-15, 0.001, 0.005, False, 200.0)
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        path = specification.write_gate_status_json(decision, Path(d) / "g.json")
        doc = json.loads(path.read_text())
    assert doc["status"] == decision.status
    assert doc["reason"] == decision.reason
