"""Phase 8C Sec. 22: the phase-only justification test must never invent a threshold."""

from src.phase8c.phase_only_test import phase_only_justification_test
from src.phase8c.scattering_analysis import CandidateCoverage


def test_phase_only_test_reports_unresolved_not_a_forced_yes_or_no():
    coverage = CandidateCoverage("S11xx", 324, 24.84, 335.16, 0.23, 0.83, 0.5, 0.1, 0.2)
    result = phase_only_justification_test(coverage)
    assert result.status == "PHASE_ONLY_APPROXIMATION_UNRESOLVED"
    assert result.status not in ("PHASE_ONLY_APPROXIMATION_JUSTIFIED", "PHASE_ONLY_APPROXIMATION_NOT_JUSTIFIED")
    assert result.threshold_available is False


def test_phase_only_test_reports_the_full_magnitude_data_even_though_unresolved():
    coverage = CandidateCoverage("S21xx", 324, 43.85, 316.15, 0.556, 0.973, 0.8, 0.05, None)
    result = phase_only_justification_test(coverage)
    assert result.magnitude_range == (0.556, 0.973)
    assert "neither Phase 8B candidate formulation" in result.note
