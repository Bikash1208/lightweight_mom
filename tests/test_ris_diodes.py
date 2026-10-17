"""Phase 20 R2 item 1: PIN-diode equivalent circuits -- passivity, limits, resonance, provenance."""

import numpy as np
import pytest

from src.ris.diodes import (
    ASSUMED,
    NOT_AVAILABLE,
    PROVISIONAL,
    DiodeModel,
    bias_current_budget,
    smp1345_off,
    smp1345_on,
)

F = np.linspace(1e9, 20e9, 191)


@pytest.mark.parametrize(
    "m",
    [
        smp1345_on(),
        smp1345_on("1mA"),
        smp1345_off(),
        smp1345_off(3e3),
        smp1345_off(10e3),
        DiodeModel("g", 2.0, 1e-9, 0.2e-12, 5e3),
    ],
)
def test_passive_over_1_to_20_ghz(m):
    assert np.all(m.impedance(F).real >= 0.0)


def test_limits_give_short_and_open():
    short = DiodeModel("s", R_s=0.0, L_s=0.0).impedance(10e9)
    assert short == 0
    assert abs(DiodeModel("c", C_T=1.0).impedance(10e9)) < 1e-9  # C_T -> inf: short
    assert abs(DiodeModel("c", C_T=1e-21).impedance(10e9)) > 1e7  # C_T -> 0: open
    assert abs(smp1345_off().scaled(C_T=1e-6).impedance(10e9)) > 1e7


def test_off_series_resonance_equals_one_over_two_pi_sqrt_LC_and_is_zero_reactance():
    m = smp1345_off()
    f0 = 1 / (2 * np.pi * np.sqrt(0.7e-9 * 0.19e-12))
    assert m.series_resonance_hz == pytest.approx(f0, rel=1e-12)
    assert abs(m.impedance(f0).imag) < 1e-9 * abs(m.impedance(0.5 * f0))
    assert m.impedance(0.8 * f0).imag < 0 < m.impedance(1.2 * f0).imag  # capacitive below, inductive above
    assert smp1345_on().series_resonance_hz is None


def test_on_state_is_Rs_plus_jwLs():
    z = smp1345_on().impedance(10e9)
    assert z.real == 1.5 and z.imag == pytest.approx(2 * np.pi * 10e9 * 0.7e-9, rel=1e-12)
    assert smp1345_on("1mA").impedance(10e9).real == 3.5


def test_rp_in_parallel_with_C_and_lossless_default():
    m = smp1345_off(3e3)
    w = 2 * np.pi * 10e9
    want = 1j * w * 0.7e-9 + 1 / (1 / 3e3 + 1j * w * 0.19e-12)
    assert m.impedance(10e9) == pytest.approx(want, rel=1e-12)
    assert smp1345_off().impedance(10e9).real == 0.0
    with pytest.raises(ValueError):
        DiodeModel("x", R_p=1e3)


def test_provenance_labels_and_scaling():
    on, off = smp1345_on(), smp1345_off()
    assert all(p.status == PROVISIONAL for p in on.provenance.values())
    assert off.provenance["R_p"].status == NOT_AVAILABLE and off.provenance["R_p"].value is None
    assert smp1345_off(3e3).provenance["R_p"].status == ASSUMED
    assert "extrapolation" in on.provenance["R_s"].condition
    s = off.scaled(C_T=1.2)
    assert pytest.approx(0.19e-12 * 1.2) == s.C_T and s.provenance["C_T"].status == "PERTURBED_SENSITIVITY"
    assert off.C_T == 0.19e-12  # original untouched
    with pytest.raises(ValueError):
        off.scaled(R_p=2.0)


def test_input_validation():
    for bad in (dict(R_s=-1.0), dict(L_s=float("nan")), dict(C_T=0.0), dict(C_T=-1e-12)):
        with pytest.raises(ValueError):
            DiodeModel("b", **bad)
    with pytest.raises(ValueError):
        smp1345_on().impedance(0.0)
    with pytest.raises(ValueError):
        smp1345_on("5mA")


def test_bias_budget():
    i, p = bias_current_budget(5685, 0.5)
    assert i == pytest.approx(5685 * 0.5 * 10e-3) and p == pytest.approx(i * 0.89)
