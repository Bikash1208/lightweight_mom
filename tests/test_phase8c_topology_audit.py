"""Phase 8C Sec. 6-10: topology audit, ground-plane test, architecture classification."""

import pytest

from src.phase8c import power_balance as pb
from src.phase8c import topology_audit as ta


def test_ground_plane_is_absent_for_a_config_with_no_ground_key():
    context = {"substrate": {"eps_r": 2.2}, "metal": {"model": "pec"}, "reference_planes": {}}
    label, evidence = ta._ground_plane_test(context)
    assert label == "ABSENT"
    assert "PEC" in evidence or "pec" in evidence


def test_ground_plane_is_present_for_a_pec_backed_context():
    """Phase 12A: the audit is now a real per-database read of context['backing'], not a
    project-wide constant."""
    context = {
        "substrate": {"eps_r": 2.2}, "metal": {"model": "pec"}, "reference_planes": {},
        "backing": "pec",
    }
    label, evidence = ta._ground_plane_test(context)
    assert label == "PRESENT"
    assert "pec" in evidence


def test_ground_plane_unknown_for_an_unrecognized_backing_value():
    context = {"substrate": {}, "metal": {}, "reference_planes": {}, "backing": "copper"}
    label, evidence = ta._ground_plane_test(context)
    assert label == "UNKNOWN"


def test_audit_topology_raises_without_context():
    class _Fake:
        metadata = {}
        cells = ()

    with pytest.raises(ValueError, match="context"):
        ta.audit_topology(_Fake(), {"mean_R_co": 0, "mean_R_cross": 0, "mean_T_co": 0, "mean_T_cross": 0})


def test_architecture_classification_is_never_pure_reflective_or_transmissive():
    label, note = ta._architecture_classification(0.9, 0.0, 0.05, 0.0)
    assert label == "PARTIALLY_REFLECTIVE_TRANSMISSIVE"
    assert "more reflective" in note


def test_architecture_classification_transmissive_case():
    label, note = ta._architecture_classification(0.05, 0.9, 0.0, 0.0)
    assert "more transmissive" in note


def test_render_topology_audit_md_includes_ground_plane_status():
    audit = ta.TopologyAudit("ABSENT", "evidence text", "metal", "substrate", "air", "PARTIALLY_REFLECTIVE_TRANSMISSIVE", "note")
    md = ta.render_topology_audit_md(audit)
    assert "ABSENT" in md
    assert "PARTIALLY_REFLECTIVE_TRANSMISSIVE" in md


def test_real_library_audit_finds_no_ground_plane(phase8c_full_database):
    db = phase8c_full_database
    summary = pb.library_power_summary(db, 100)   # index 100 = 12.5 GHz for the locked 201-pt grid
    audit = ta.audit_topology(db, summary)
    assert audit.ground_plane == "ABSENT"


def test_s_parameter_definition_md_cites_project_conventions():
    md = ta.render_s_parameter_definition_md()
    assert "exp(+j" in md or "Ports" in md
    assert "reference plane" in md.lower()
