"""Sec. 11-12: architecture compatibility, against the real Phase 7 checkpoint's own audit."""

import pytest

from src.phase8c import data_source, power_balance, topology_audit
from src.phase8e.architecture_compatibility import assess_all
from src.phase_extraction import design_frequency_index

pytestmark = pytest.mark.skipif(
    not data_source.DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present"
)


def test_formulation_a_is_architecture_incompatible_no_ground_plane(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    mean_power = power_balance.library_power_summary(phase8c_full_database, idx)
    audit = topology_audit.audit_topology(phase8c_full_database, mean_power)
    results = assess_all(audit, mean_power)
    assert results["A"].status == "ARCHITECTURE_INCOMPATIBLE"


def test_formulation_b_and_c_are_architecture_compatible_unvalidated_given_real_power_split(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    mean_power = power_balance.library_power_summary(phase8c_full_database, idx)
    audit = topology_audit.audit_topology(phase8c_full_database, mean_power)
    results = assess_all(audit, mean_power)
    assert results["B"].status == "ARCHITECTURE_COMPATIBLE_UNVALIDATED"
    assert results["C"].status == "ARCHITECTURE_COMPATIBLE_UNVALIDATED"


def test_formulation_d_is_architecture_incompatible_given_measured_negligible_cross_pol(phase8c_full_database):
    idx = design_frequency_index(phase8c_full_database, 12.5e9)
    mean_power = power_balance.library_power_summary(phase8c_full_database, idx)
    audit = topology_audit.audit_topology(phase8c_full_database, mean_power)
    results = assess_all(audit, mean_power)
    assert results["D"].status == "ARCHITECTURE_INCOMPATIBLE"


def test_formulation_d_becomes_compatible_when_cross_pol_is_not_negligible():
    from src.phase8e.architecture_compatibility import assess_formulation_d

    result = assess_formulation_d({"mean_R_cross": 0.05, "mean_T_cross": 0.01})
    assert result.status == "ARCHITECTURE_COMPATIBLE_UNVALIDATED"


def test_formulation_b_becomes_incompatible_with_zero_reflected_power():
    from src.phase8e.architecture_compatibility import assess_formulation_b

    result = assess_formulation_b({"mean_R_co": 0.0, "mean_R_cross": 0.0})
    assert result.status == "ARCHITECTURE_INCOMPATIBLE"
