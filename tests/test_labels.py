"""Phase 13 Part 2: the fixed validation-label vocabulary."""

from src import labels


def test_parameterized_labels_bake_in_the_value():
    assert labels.hfss_residual_at_resonance(3.7) == "HFSS_RESIDUAL_AT_RESONANCE(~3.7deg)"
    assert labels.phase_coverage_shortfall(186.0) == "PHASE_COVERAGE_SHORTFALL(174deg)"
    assert labels.steep_gradient(5.8) == "STEEP_GRADIENT(>5.8deg/cell)"


def test_labeled_result_accumulates_labels_and_warnings_immutably():
    r = labels.LabeledResult(value=42)
    r2 = r.with_label(labels.MESH_CONVERGED)
    r3 = r2.with_warning("mesh is coarse")
    assert r.labels == ()
    assert r2.labels == (labels.MESH_CONVERGED,)
    assert r3.labels == (labels.MESH_CONVERGED,)
    assert r3.warnings == ("mesh is coarse",)
    assert r3.value == 42


def test_has_matches_plain_and_parameterized_labels():
    r = labels.LabeledResult(value=None, labels=(labels.hfss_residual_at_resonance(3.7),))
    assert r.has("HFSS_RESIDUAL_AT_RESONANCE")
    assert not r.has(labels.MESH_CONVERGED)
