"""Phase 8B Sec. 22-25: the control-variable investigation, including a re-derivation of the
case01 dB figures it cites (from the actual solved-case file, not from the literal constants in
`control_variable.py`)."""

import json

import numpy as np
import pytest

from src.phase8b import control_variable as cv

CASE01_JSON = "results/validation_cases/case01_non_resonant_patch/mom_results.json"


def test_investigate_s11xx_is_insufficient_evidence_not_a_forced_binary_answer():
    result = cv.investigate_s11xx()
    assert result.validated == "INSUFFICIENT_EVIDENCE"
    assert result.validated not in ("VALIDATED", "INVALID")


def test_investigate_s21xx_is_insufficient_evidence():
    result = cv.investigate_s21xx()
    assert result.validated == "INSUFFICIENT_EVIDENCE"
    # Sec. 23's explicit prohibition must be honoured: larger magnitude alone is not a justification.
    assert "not select S21 merely because" in result.justification or "not sufficient justification" in result.justification


def test_investigate_surface_impedance_is_invalid_with_two_independent_reasons():
    result = cv.investigate_surface_impedance()
    assert result.validated == "INVALID"
    assert "ground plane" in result.justification
    assert "dispersion" in result.justification


def test_investigate_reflection_coefficient_is_validated():
    # Gamma_xx = S11_xx is a definitional identity, unlike the surface-impedance interpretation.
    result = cv.investigate_reflection_coefficient()
    assert result.validated == "VALIDATED"


def test_all_candidates_are_independently_evaluable():
    results = [f() for f in cv.ALL_CANDIDATES]
    assert len(results) == 4
    assert {r.name for r in results} == {"S11_xx", "S21_xx", "surface_impedance_Zs", "reflection_coefficient_Gamma"}


@pytest.mark.skipif(not __import__("pathlib").Path(CASE01_JSON).exists(), reason="case01 validation export not built")
def test_case01_db_ranges_are_re_derived_from_the_actual_file():
    doc = json.loads(open(CASE01_JSON).read().strip())
    cell = doc["cells"][0]
    mags11 = np.array([r["S"]["S11"]["xx"]["mag"] for r in cell["records"]])
    mags21 = np.array([r["S"]["S21"]["xx"]["mag"] for r in cell["records"]])
    db11 = 20 * np.log10(mags11)
    db21 = 20 * np.log10(mags21)
    # Cross-check against the literal constants control_variable.py cites (documented tolerance:
    # those constants were rounded to 2 decimals from exactly this computation).
    assert db11.min() == pytest.approx(cv.CASE01_S11XX_DB_RANGE[0], abs=0.01)
    assert db11.max() == pytest.approx(cv.CASE01_S11XX_DB_RANGE[1], abs=0.01)
    assert db21.min() == pytest.approx(cv.CASE01_S21XX_DB_RANGE[0], abs=0.01)
    assert db21.max() == pytest.approx(cv.CASE01_S21XX_DB_RANGE[1], abs=0.01)
    # The actual claim being tested: this structure is mostly transmissive, not mostly reflective.
    assert db21.min() > db11.max()
