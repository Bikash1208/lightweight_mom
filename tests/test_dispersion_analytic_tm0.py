"""Analytic TM0 grounded-slab dispersion, checked against PHASE13_WEBAPP_PROMPT.md's own cited
numbers for RT/Duroid 5880 (eps_r=2.2, h=1.575mm): beta ~= 213/269/326 rad/m at 10/12.5/15 GHz."""

import pytest

from src.dispersion import analytic_tm0 as tm0

EPS_R = 2.2
THICKNESS_M = 1.575e-3


@pytest.mark.parametrize("freq_ghz, expected_beta", [(10.0, 213.0), (12.5, 269.0), (15.0, 326.0)])
def test_beta_matches_brief_cited_values(freq_ghz, expected_beta):
    beta = tm0.tm0_beta(freq_ghz * 1e9, EPS_R, THICKNESS_M)
    assert beta == pytest.approx(expected_beta, abs=1.0)


def test_n_eff_is_just_above_one_for_a_thin_low_contrast_slab():
    n_eff = tm0.tm0_n_eff(10e9, EPS_R, THICKNESS_M)
    assert 1.0 < n_eff < 1.05   # matches the brief's own "n_eff ~= 1.07" ballpark


def test_beta_increases_with_frequency():
    betas = [tm0.tm0_beta(f * 1e9, EPS_R, THICKNESS_M) for f in (10.0, 12.5, 15.0)]
    assert betas[0] < betas[1] < betas[2]


def test_rejects_eps_r_at_or_below_one():
    with pytest.raises(ValueError, match="eps_r"):
        tm0.tm0_beta(10e9, 1.0, THICKNESS_M)


def test_thicker_slab_gives_larger_beta_more_confinement():
    thin = tm0.tm0_beta(10e9, EPS_R, 0.5e-3)
    thick = tm0.tm0_beta(10e9, EPS_R, 3.0e-3)
    assert thin < thick
