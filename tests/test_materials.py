import math

import pytest

from src.constants import C0
from src.materials import Metal, Substrate


def rt5880(**kw):
    """Test fixture values (the solver itself has no built-in defaults)."""
    args = dict(name="test", eps_r=2.2, tan_delta=0.0009, thickness=1.575e-3)
    args.update(kw)
    return Substrate(**args)


def test_no_builtin_electrical_defaults():
    with pytest.raises(TypeError):
        Substrate()
    with pytest.raises(TypeError):
        Substrate(name="x")


def test_complex_permittivity_definition():
    s = rt5880()
    assert s.eps_complex.real == pytest.approx(2.2)
    assert s.eps_complex.imag == pytest.approx(-2.2 * 0.0009)
    # loss tangent is recovered as -Im/Re under exp(+jwt)
    assert -s.eps_complex.imag / s.eps_complex.real == pytest.approx(0.0009, rel=1e-14)


def test_lossless_is_purely_real():
    s = rt5880(tan_delta=0.0)
    assert s.eps_complex.imag == 0.0
    assert s.refractive_index == pytest.approx(math.sqrt(2.2), rel=1e-15)
    assert s.refractive_index.imag == 0.0


@pytest.mark.parametrize("eps_r, tan_d", [(2.2, 0.0009), (10.2, 0.0023), (4.4, 0.02), (1.0, 0.5)])
def test_refractive_index_squares_back_and_is_passive(eps_r, tan_d):
    s = rt5880(eps_r=eps_r, tan_delta=tan_d)
    n = s.refractive_index
    assert n * n == pytest.approx(s.eps_complex, rel=1e-14)
    assert n.real > 0.0
    assert n.imag < 0.0  # passive medium under exp(+jwt)


def test_substrate_wavenumber_is_k0_times_index():
    s = rt5880(tan_delta=0.0)
    f = 10e9
    k = s.wavenumber(f)
    assert k.imag == 0.0
    assert k.real == pytest.approx(2 * math.pi * f * math.sqrt(2.2) / C0, rel=1e-14)
    # wavelength in medium = lambda0 / sqrt(eps_r)
    assert 2 * math.pi / k.real == pytest.approx((C0 / f) / math.sqrt(2.2), rel=1e-14)


def test_substrate_wavenumber_lossy_attenuates():
    s = rt5880()
    k = s.wavenumber(10e9)
    assert k.imag < 0.0
    # small-loss approximation Im(k) ~ -k' * tan_delta / 2
    assert k.imag == pytest.approx(-k.real * 0.0009 / 2, rel=1e-3)


@pytest.mark.parametrize(
    "kw",
    [
        dict(eps_r=0.0),
        dict(eps_r=-2.2),
        dict(eps_r=float("nan")),
        dict(eps_r=float("inf")),
        dict(eps_r="2.2"),
        dict(eps_r=True),
        dict(eps_r=2.2 + 0j),
        dict(tan_delta=-1e-4),  # would be an active (gain) medium
        dict(tan_delta=float("nan")),
        dict(thickness=0.0),
        dict(thickness=-1e-3),
        dict(thickness=float("nan")),
    ],
)
def test_substrate_rejects_invalid(kw):
    with pytest.raises(ValueError):
        rt5880(**kw)


def test_metal_only_pec_in_phase1():
    assert Metal().model == "pec"
    assert Metal("pec").model == "pec"
    for model in ("copper", "finite_conductivity", ""):
        with pytest.raises(ValueError):
            Metal(model)
