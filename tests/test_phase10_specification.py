"""Sec. 45, 47: the required CSV/JSON writers, exercised against small synthetic objects (no
solves needed -- these functions only format already-computed results)."""

import json

from src.phase10 import specification as spec_mod
from src.phase10.dataset import ALL_POSITIONS, EnvironmentSample
from src.phase10.error_metrics import compute_error_distribution
from src.phase10.inverse import InverseResult


def _sample(sample_id, center, split_source="structured"):
    neighbors = {p: center for p in ALL_POSITIONS}
    return EnvironmentSample(sample_id, "C3", center, neighbors, "test", 12.5e9, 1.0 + 1.0j, split_source)


def test_write_dataset_split_csv_only_writes_the_requested_split(tmp_path):
    samples = [_sample("a", (2.8, 2.8)), _sample("b", (2.8, 2.8))]
    splits = {"a": "train", "b": "validation"}
    path = spec_mod.write_dataset_split_csv(samples, splits, "train", tmp_path / "train.csv")
    rows = path.read_text().strip().split("\n")
    assert len(rows) == 2  # header + sample "a" only
    assert rows[1].startswith("a,")


def test_write_dataset_metadata_json_reports_split_counts(tmp_path):
    samples = [_sample("a", (2.8, 2.8)), _sample("b", (2.8, 2.8)), _sample("c", (4.2, 4.2))]
    splits = {"a": "train", "b": "validation", "c": "blind_test"}
    path = spec_mod.write_dataset_metadata_json(samples, splits, tmp_path / "meta.json")
    doc = json.loads(path.read_text())
    assert doc["n_train"] == 1 and doc["n_validation"] == 1 and doc["n_blind_test"] == 1
    assert doc["n_total"] == 3


def test_write_error_metrics_csv_covers_every_model(tmp_path):
    report = compute_error_distribution([1.0 + 0j], [1.1 + 0j])
    reports = {"baseline1": report, "baseline2": report, "surrogate": report}
    path = spec_mod.write_error_metrics_csv(reports, tmp_path / "err.csv")
    text = path.read_text()
    for model in reports:
        assert model in text


def test_write_target_responses_and_mom_verification_csvs_do_not_shadow_the_writer(tmp_path):
    """Regression test for a real bug caught during development: an inner ``l, w = ...`` unpack
    once shadowed the csv ``writer`` variable named ``w`` and crashed on the second row."""
    results = [
        InverseResult(
            1.0 + 1.0j, {p: (2.8, 2.8) for p in ALL_POSITIONS}, (2.8, 2.8), 1.0 + 1.0j, 1.0 + 1.0j, 0.0, "REACHED"
        ),
        InverseResult(
            2.0 + 2.0j, {p: (1.0, 4.4) for p in ALL_POSITIONS}, None, None, None, float("inf"), "TARGET_UNREACHABLE"
        ),
    ]
    p1 = spec_mod.write_target_responses_csv(results, tmp_path / "targets.csv")
    p2 = spec_mod.write_inverse_results_csv(results, tmp_path / "inverse.csv")
    p3 = spec_mod.write_mom_verification_csv(results, tmp_path / "mom.csv")
    for p in (p1, p2, p3):
        rows = p.read_text().strip().split("\n")
        assert len(rows) == 3  # header + 2 result rows


def test_write_phase10_summary_and_gate_status_json_are_consistent(tmp_path):
    from src.phase10 import gate as gate_mod

    keys_and_evals = [
        (gate_mod.GATE_10A_KEYS, gate_mod.evaluate_gate_10a),
        (gate_mod.GATE_10B_KEYS, gate_mod.evaluate_gate_10b),
        (gate_mod.GATE_10C_KEYS, gate_mod.evaluate_gate_10c),
        (gate_mod.GATE_10D_KEYS, gate_mod.evaluate_gate_10d),
        (gate_mod.GATE_10E_KEYS, gate_mod.evaluate_gate_10e),
    ]
    gates = [ev({k: True for k in keys}) for keys, ev in keys_and_evals]
    decision = gate_mod.decide(
        *gates, surrogate_beats_baseline1=True, inverse_mostly_reachable=True, beam_level_acceptable=True
    )
    summary_path = spec_mod.write_phase10_summary_json(decision, tmp_path / "summary.json")
    status_path = spec_mod.write_phase10_gate_status_json(decision, tmp_path / "status.json")
    summary = json.loads(summary_path.read_text())
    status = json.loads(status_path.read_text())
    assert summary["final_status"] == status["final_status"]
    assert summary["geometry_inversion"]["status"] == "UNLOCKED"
