"""Phase 12C gates: explicit named states, and the two negative-test convention this project's
gates all follow (a clearly-failing input must not be classified as passing)."""

from src.phase12c import gate
from src.phase12c.locality_study import GradientCaseResult


def _case(name, achieved_step_deg, phase_error_deg, isolated_px, embedded_px, relative_complex_error=0.0):
    return GradientCaseResult(
        name=name,
        n=3,
        center_lw_mm=(2.0, 2.0),
        grid_mm=(),
        isolated_px=isolated_px,
        embedded_px=embedded_px,
        complex_error=0.0,
        magnitude_error=0.0,
        phase_error_deg=phase_error_deg,
        relative_complex_error=relative_complex_error,
        achieved_step_deg=achieved_step_deg,
        n_unknowns_isolated=10,
        n_unknowns_embedded=90,
    )


def test_magnitude_error_db_is_zero_for_identical_magnitudes():
    assert gate.magnitude_error_db(1.0 + 0j, 1.0 + 0j) == 0.0


def test_every_case_within_tolerance_gives_valid_status():
    cases = [_case("a", 5.0, 1.0, 1.0 + 0j, 1.0 * (0.99 + 0.001j))]
    result = gate.decide_local_periodicity(cases)
    assert result.status == gate.GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS
    assert result.max_step_within_tolerance_deg == 5.0


def test_no_case_within_tolerance_gives_invalidated_status():
    """Mandatory negative test: a case with a large phase error must not be classified VALID."""
    cases = [_case("bad", 100.0, 90.0, 1.0 + 0j, 0.1 - 0.1j, relative_complex_error=0.9)]
    result = gate.decide_local_periodicity(cases)
    assert result.status == gate.GROUNDED_LOCAL_PERIODICITY_INVALIDATED
    assert result.max_step_within_tolerance_deg is None


def test_mixed_cases_report_the_largest_step_that_still_holds():
    good = _case("good", 5.0, 1.0, 1.0 + 0j, 0.99 + 0.001j)
    bad = _case("bad", 80.0, 90.0, 1.0 + 0j, 0.1 - 0.1j, relative_complex_error=0.9)
    result = gate.decide_local_periodicity([good, bad])
    assert result.status == gate.GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS
    assert result.max_step_within_tolerance_deg == 5.0


def test_spot_check_gate_reports_changed_when_c2_c3_mask_is_no_longer_identical():
    from src.phase12c.spot_check import C2C3Recheck, CaseDRecheck, CouplingSpotCheckPoint, Phase9SpotCheck

    c2_c3 = C2C3Recheck(24, False, 1.0 + 0j, 1.5 + 0j, 0.5)
    p9 = Phase9SpotCheck(48, 3, 1.0 + 0j, 0.01 + 0j, 0.99, 10, 90, 0.997)
    case_d = CaseDRecheck(48, 3, 1.0 + 0j, 2.0 + 0j, 1.0, True)
    result = gate.decide_spot_check(
        c2_c3, p9, case_d, [CouplingSpotCheckPoint(1, 1.0 + 0j, 10), CouplingSpotCheckPoint(3, 1.1 + 0j, 90)]
    )
    assert result.status == gate.EARLIER_CONCLUSIONS_CHANGED
    names = [f[0] for f in result.findings]
    assert "phase10_c2_c3" in names


def test_spot_check_gate_reports_hold_when_nothing_changed():
    from src.phase12c.spot_check import C2C3Recheck, CaseDRecheck, CouplingSpotCheckPoint, Phase9SpotCheck

    c2_c3 = C2C3Recheck(24, True, 1.0 + 0j, 1.0 + 0j, 0.0)
    p9 = Phase9SpotCheck(48, 3, 1.0 + 0j, 0.001 + 0j, 0.999, 10, 90, 0.997)
    case_d = CaseDRecheck(48, 3, 1.0 + 0j, 2.0 + 0j, 1.0, True)
    result = gate.decide_spot_check(
        c2_c3, p9, case_d, [CouplingSpotCheckPoint(1, 1.0 + 0j, 10), CouplingSpotCheckPoint(3, 2.0 + 0j, 90)]
    )
    assert result.status == gate.EARLIER_CONCLUSIONS_HOLD


def test_hfss_export_status_requires_every_file_to_exist(tmp_path):
    from src.phase12c.exports import SingleCellExport

    missing = SingleCellExport(1.0, str(tmp_path / "nonexistent.step"), 0.9, -10.0, 8, 8)
    assert gate.hfss_export_status([missing], [], None) == "HFSS_EXPORT_NOT_READY"
