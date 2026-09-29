"""Phase 8B Sec. 19-20: the two candidate modulation kernels, numerically tested independent of
whether either is validated for this project (that is `control_variable`/`library_compatibility`'s job)."""

import numpy as np
import pytest

from src.phase8b.holographic_modulation import (
    interference_modulation,
    phase_difference_hologram,
    phase_difference_hologram_from_phases_deg,
    reactance_from_modulation,
)


def test_interference_modulation_is_real_and_bounded_for_unit_amplitude_fields():
    e_ref = np.exp(1j * np.radians([0, 45, 90, 180, -90]))
    e_obj = np.exp(1j * np.radians([0, 0, 0, 0, 0]))
    psi = interference_modulation(e_ref, e_obj)
    assert np.all(np.isreal(psi))
    assert np.all(psi <= 1.0 + 1e-12) and np.all(psi >= -1.0 - 1e-12)


def test_interference_modulation_peaks_when_fields_are_in_phase():
    e_ref_aligned = np.array([1.0 + 0.0j])
    e_obj_aligned = np.array([1.0 + 0.0j])
    e_obj_opposed = np.array([-1.0 + 0.0j])
    assert interference_modulation(e_ref_aligned, e_obj_aligned)[0] == pytest.approx(1.0)
    assert interference_modulation(e_ref_aligned, e_obj_opposed)[0] == pytest.approx(-1.0)


def test_reactance_from_modulation_rejects_out_of_range_index():
    with pytest.raises(ValueError):
        reactance_from_modulation(np.array([0.5]), 100.0, modulation_index=1.5)
    with pytest.raises(ValueError):
        reactance_from_modulation(np.array([0.5]), 100.0, modulation_index=0.0)


def test_reactance_from_modulation_is_linear_in_psi():
    psi = np.array([-1.0, 0.0, 1.0])
    x = reactance_from_modulation(psi, x_avg_ohm=200.0, modulation_index=0.3)
    assert x[1] == pytest.approx(200.0)
    assert x[0] == pytest.approx(200.0 * (1 - 0.3))
    assert x[2] == pytest.approx(200.0 * (1 + 0.3))


def test_phase_difference_hologram_recovers_a_known_phase_shift():
    e_ref = np.exp(1j * np.radians([10.0, -170.0]))
    e_obj = np.exp(1j * np.radians([40.0, 170.0]))
    required = phase_difference_hologram(e_ref, e_obj)
    assert required[0] == pytest.approx(30.0, abs=1e-6)
    # -170 -> wrap(170-(-170)) = wrap(340) = -20
    assert required[1] == pytest.approx(-20.0, abs=1e-6)


def test_phase_difference_hologram_from_phases_matches_complex_version():
    ref_deg = np.array([12.3, -45.6, 179.9])
    obj_deg = np.array([88.8, 45.6, -179.9])
    a = phase_difference_hologram_from_phases_deg(ref_deg, obj_deg)
    b = phase_difference_hologram(np.exp(1j * np.radians(ref_deg)), np.exp(1j * np.radians(obj_deg)))
    assert np.allclose(a, b, atol=1e-6)


def test_phase_difference_hologram_result_is_always_wrapped():
    ref_deg = np.array([0.0])
    obj_deg = np.array([540.0])   # 540 wraps to 180
    required = phase_difference_hologram_from_phases_deg(ref_deg, obj_deg)
    assert -180.0 < required[0] <= 180.0
