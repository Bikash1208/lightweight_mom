"""Phase 12D gate logic: explicit named states, never a bare pass/fail."""

from dataclasses import dataclass

from src.phase12d import gate


@dataclass
class _FakeStepCase:
    name: str
    achieved_step_deg: float
    reflection_phase_error_deg: float
    reflection_magnitude_error_db: float


def test_decide_step_limit_measured_when_some_case_within_tolerance():
    cases = [
        _FakeStepCase("a", 2.0, 0.5, 0.05),
        _FakeStepCase("b", 50.0, 20.0, 3.0),
    ]
    result = gate.decide_step_limit(cases)
    assert result.status == gate.LOCAL_PERIODICITY_STEP_LIMIT_MEASURED
    assert result.max_step_within_design_tolerance_deg == 2.0


def test_decide_step_limit_not_resolved_when_nothing_within_tolerance():
    cases = [_FakeStepCase("a", 80.0, 30.0, 5.0)]
    result = gate.decide_step_limit(cases)
    assert result.status == gate.LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED
    assert result.max_step_within_design_tolerance_deg is None


def test_decide_step_limit_empty_input():
    result = gate.decide_step_limit([])
    assert result.status == gate.LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED


def test_decide_step_limit_reports_strict_tolerance_separately():
    cases = [_FakeStepCase("a", 2.0, 7.0, 0.7)]  # within design (10/1) but not strict (5/0.5)
    result = gate.decide_step_limit(cases)
    assert result.status == gate.LOCAL_PERIODICITY_STEP_LIMIT_MEASURED
    assert result.max_step_within_design_tolerance_deg == 2.0
    assert result.max_step_within_strict_tolerance_deg is None


@dataclass
class _FakeApertureResult:
    n_cells: int
    fraction_intended: float
    spurious_orders: tuple


def test_decide_aperture_adequacy_adequate():
    results = [_FakeApertureResult(6, 0.95, (((1, 1), 0.01),))]
    g = gate.decide_aperture_adequacy(results, requested_sizes=(6,))
    assert g.status == gate.LOOKUP_TABLE_ADEQUATE


def test_decide_aperture_adequacy_inadequate():
    results = [_FakeApertureResult(6, 0.3, (((1, 1), 0.2),))]
    g = gate.decide_aperture_adequacy(results, requested_sizes=(6,))
    assert g.status == gate.LOOKUP_TABLE_INADEQUATE


def test_decide_aperture_adequacy_records_infeasible_size_without_dropping_it():
    results = [_FakeApertureResult(6, 0.95, ())]
    g = gate.decide_aperture_adequacy(results, requested_sizes=(6, 8, 12))
    attempted = {n: a for n, _, _, _, a in g.per_size}
    assert attempted == {6: True, 8: False, 12: False}


def test_decide_aperture_adequacy_untested_when_none_feasible():
    g = gate.decide_aperture_adequacy([], requested_sizes=(6, 8, 12))
    assert g.status == gate.LOOKUP_TABLE_UNTESTED


@dataclass
class _FakeComparison:
    case_id: str
    meets_tolerance: bool


def test_decide_hfss_status_partial_when_one_case_validated():
    g = gate.decide_hfss_status([_FakeComparison("G1", True)], ("G1", "G2", "G3", "G4", "G5", "G6"), ("G4",))
    assert g.status == gate.HFSS_SINGLE_CELL_VALIDATED_PARTIAL
    assert g.cases_with_data == ("G1",)
    assert "G2" in g.cases_pending


def test_decide_hfss_status_not_validated_when_no_data():
    g = gate.decide_hfss_status([], ("G1",), ())
    assert g.status == "HFSS_SINGLE_CELL_NOT_VALIDATED"
