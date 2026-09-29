"""Phase 8B Sec. 8-10, 43-44: reference-wave input types and the monopole/guided-surface-wave
field, including the honest refusal to substitute a missing `k_swg`."""

import numpy as np
import pytest

from src.phase8b.reference_wave import (
    ReferenceField,
    ReferenceFieldKind,
    monopole_reference_field,
    phase_only_reference_requires_amplitude_assumption,
    surface_impedance_reference_sufficiency,
)


def test_reference_field_requires_matching_shapes():
    with pytest.raises(ValueError):
        ReferenceField(
            ReferenceFieldKind.PHASE_ONLY, values=np.array([1.0, 2.0]), x_m=np.array([0.0]), y_m=np.array([0.0])
        )


def test_reference_field_kinds_are_distinguished():
    rf = ReferenceField(
        ReferenceFieldKind.SURFACE_IMPEDANCE, values=np.array([1.0 + 2.0j]), x_m=np.array([0.0]), y_m=np.array([0.0])
    )
    assert rf.kind is ReferenceFieldKind.SURFACE_IMPEDANCE


def test_monopole_reference_field_requires_k_swg():
    with pytest.raises(ValueError, match="k_swg"):
        monopole_reference_field(np.array([0.01]), np.array([0.0]), k_swg=None)


def test_monopole_reference_field_has_cylindrical_1_over_sqrt_rho_falloff():
    k_swg = 2 * np.pi * 12.5e9 / 3e8  # an arbitrary illustrative value, not claimed physical
    rho = np.array([0.01, 0.04])  # rho2 = 4 * rho1
    field = monopole_reference_field(rho, np.zeros_like(rho), k_swg=k_swg)
    # amplitude ~ rho^-1/2, so |field(rho1)| / |field(rho2)| = sqrt(rho2/rho1) = 2
    assert abs(field[0]) / abs(field[1]) == pytest.approx(2.0, rel=1e-6)


def test_monopole_reference_field_phase_matches_exp_minus_jkr():
    k_swg = 50.0
    rho = np.array([0.02])
    field = monopole_reference_field(rho, np.zeros_like(rho), k_swg=k_swg)
    assert np.angle(field[0]) == pytest.approx(float(np.angle(np.exp(-1j * k_swg * rho[0]))), abs=1e-9)


def test_phase_only_reference_sufficiency_differs_by_formulation():
    result = phase_only_reference_requires_amplitude_assumption()
    assert "UNSUPPORTED" in result["formulation_A_fong2010_interference_kernel"]
    assert "sufficient" in result["formulation_B_phase_difference_hologram"]


def test_surface_impedance_reference_sufficiency_is_explicit():
    note = surface_impedance_reference_sufficiency()
    assert "REQUIRES_ADDITIONAL_DATA" in note
