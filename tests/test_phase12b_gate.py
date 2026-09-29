"""Phase 12B gate/status-code logic."""

from src.phase12b import gate


def test_mesh_alignment_validated_within_tolerance():
    g = gate.decide_mesh_alignment(0.005, 0.001)
    assert g.status == gate.ALIGNED_MESH_VALIDATED


def test_mesh_alignment_invalidated_over_tolerance():
    g = gate.decide_mesh_alignment(0.05, 0.001)
    assert g.status == gate.ALIGNED_MESH_INVALIDATED


def test_audit_no_conclusion_changed_when_nothing_found():
    g = gate.decide_audit({}, {"leaks": ()}, ())
    assert g.status == gate.NO_CONCLUSION_CHANGED
    assert g.at_risk_findings == ()


def test_audit_conclusions_at_risk_when_duplicates_found():
    g = gate.decide_audit({}, {"leaks": ()}, (("C2", "C3"),))
    assert g.status == gate.CONCLUSIONS_AT_RISK
    assert len(g.at_risk_findings) == 1


def test_period_study_sufficient_above_threshold():
    g = gate.decide_period_study({4.8: 320.0}, threshold_deg=300.0)
    assert g.status == gate.COVERAGE_SUFFICIENT


def test_period_study_insufficient_below_threshold():
    g = gate.decide_period_study({4.8: 245.0, 8.0: 200.0}, threshold_deg=300.0)
    assert g.status == gate.COVERAGE_INSUFFICIENT
    assert "4.8" in g.reason


def test_hfss_status_reexported_and_not_available_without_a_file(tmp_path):
    assert gate.hfss_status(tmp_path / "none.csv") == gate.HFSS_NOT_AVAILABLE
