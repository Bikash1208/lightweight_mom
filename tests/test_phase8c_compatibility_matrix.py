"""Phase 8C Sec. 37-38, 53, 62: the factual compatibility matrix (never a ranking) and the final
status decision tree."""

from src.phase8c.compatibility_matrix import FINAL_STATUS_TAXONOMY, final_status, library_adequacy


def test_final_status_is_never_a_bare_pass():
    assert "PASS" not in FINAL_STATUS_TAXONOMY


def test_library_adequacy_unvalidated_when_control_variable_not_validated():
    assert library_adequacy("INSUFFICIENT_EVIDENCE", 24.84) == "UNVALIDATED"


def test_library_adequacy_incompatible_when_control_variable_invalid():
    assert library_adequacy("INVALID", None) == "INCOMPATIBLE"


def test_final_status_architecture_unresolved_short_circuits_everything():
    status = final_status(
        architecture_resolved=False, s11_investigation=None, s21_investigation=None,
        formulation_a_status="LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE",
        reciprocity_ok=True, propagating_modes_ok=True,
    )
    assert status == "ARCHITECTURE_UNRESOLVED"


def test_final_status_scattering_data_insufficient_when_modes_or_reciprocity_fail():
    class _Inv:
        validated = "INSUFFICIENT_EVIDENCE"

    for reciprocity_ok, modes_ok in [(False, True), (True, False)]:
        status = final_status(
            architecture_resolved=True, s11_investigation=_Inv(), s21_investigation=_Inv(),
            formulation_a_status="LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE",
            reciprocity_ok=reciprocity_ok, propagating_modes_ok=modes_ok,
        )
        assert status == "SCATTERING_DATA_INSUFFICIENT"


def test_final_status_control_variable_unvalidated_for_this_projects_real_finding():
    class _Inv:
        validated = "INSUFFICIENT_EVIDENCE"

    status = final_status(
        architecture_resolved=True, s11_investigation=_Inv(), s21_investigation=_Inv(),
        formulation_a_status="LIBRARY_CONTROL_VARIABLE_INCOMPATIBLE",
        reciprocity_ok=True, propagating_modes_ok=True,
    )
    assert status == "CONTROL_VARIABLE_UNVALIDATED"


def test_final_status_does_not_reward_large_phase_coverage_alone():
    # Sec. 52: a "validated" outcome requires physical validity, not phase coverage; this test
    # exists to document that `final_status`'s signature never even takes a coverage number as an
    # argument capable of flipping the decision -- it cannot pass "because coverage is large".
    import inspect

    sig = inspect.signature(final_status)
    assert "coverage" not in " ".join(sig.parameters).lower()
