"""Sec. 9-10: the exact complex quantity C_req(x,y) each candidate formulation derives from
(E_ref, E_obj), and the measurable unit-cell quantity each one asks to equal it. Candidates A/B/C
reuse `phase8b.holographic_modulation` unmodified; Candidate D is new (Sec. 5's generalized
vector/matrix operator).
"""

from dataclasses import dataclass

import numpy as np

from ..phase8b.holographic_modulation import (
    interference_modulation,
    phase_difference_hologram,
    reactance_from_modulation,
)


@dataclass(frozen=True)
class ControlResponseResult:
    formulation: str                  # "A" | "B" | "C" | "D"
    value_kind: str                    # "complex_reactance" | "phase_deg" | "transfer_operator"
    values: object                     # array (or dict of arrays for D)
    required_quantity: str
    mathematical_derivation: str
    physical_meaning: str
    dimensionality: str
    phase_convention: str
    required_physical_surface_response: str
    assumptions: tuple
    limitations: tuple


def required_reactance_a(e_ref, e_obj, x_avg_ohm, modulation_index):
    """Candidate A: X(x,y) = X_avg (1 + M * Re{E_ref . conj(E_obj)}) -- `interference_modulation`/
    `reactance_from_modulation` (Phase 8B, unmodified)."""
    psi = interference_modulation(e_ref, e_obj)
    values = reactance_from_modulation(psi, x_avg_ohm, modulation_index)
    return ControlResponseResult(
        formulation="A", value_kind="complex_reactance", values=values,
        required_quantity="surface reactance X(x,y) of a grounded modulated-impedance boundary",
        mathematical_derivation="Psi = Re{E_ref . conj(E_obj)}; X = X_avg(1 + M*Psi) (Fong et al. 2010)",
        physical_meaning="local reactive boundary-condition coefficient of a grounded surface",
        dimensionality="real-valued [ohm]",
        phase_convention="exp(+j*omega*t); Psi built from complex E_ref, E_obj under exp(-j*k.r)",
        required_physical_surface_response="grounded surface reactance (not available, see architecture_compatibility.py)",
        assumptions=("grounded structure", "|Gamma|=1 (purely reactive boundary)", "M, X_avg known"),
        limitations=("M/X_avg cannot be calibrated from this project's ungrounded library",),
    )


def required_phase_b(e_ref, e_obj):
    """Candidate B: required reflection phase = wrap(arg(E_obj) - arg(E_ref))."""
    values = phase_difference_hologram(e_ref, e_obj)
    return ControlResponseResult(
        formulation="B", value_kind="phase_deg", values=values,
        required_quantity="normal-incidence co-polarized reflection phase (arg(S11xx))",
        mathematical_derivation="required_phase = wrap(arg(E_obj) - arg(E_ref))",
        physical_meaning="phase shift the reflected wave must add to turn E_ref into E_obj locally",
        dimensionality="real-valued [deg], phase only",
        phase_convention="exp(+j*omega*t); exp(-j*k.r); principal range (-180, 180] (phase.wrap_deg)",
        required_physical_surface_response="arg(S11xx) at the design frequency, normal incidence, x-pol",
        assumptions=("reflective architecture", "local periodic approximation", "amplitude of E_ref does not matter"),
        limitations=("phase-only: does not control |S11xx|; amplitude mismatch is not represented",),
    )


def required_phase_c(e_ref, e_obj):
    """Candidate C: the identical kernel as B, claimed against the transmission phase instead."""
    values = phase_difference_hologram(e_ref, e_obj)
    return ControlResponseResult(
        formulation="C", value_kind="phase_deg", values=values,
        required_quantity="normal-incidence co-polarized transmission phase (arg(S21xx))",
        mathematical_derivation="required_phase = wrap(arg(E_obj) - arg(E_ref))",
        physical_meaning="phase shift the transmitted wave must add to turn E_ref into E_obj locally",
        dimensionality="real-valued [deg], phase only",
        phase_convention="exp(+j*omega*t); exp(-j*k.r); principal range (-180, 180] (phase.wrap_deg)",
        required_physical_surface_response="arg(S21xx) at the design frequency, normal incidence, x-pol",
        assumptions=("transmissive architecture", "local periodic approximation", "amplitude of E_ref does not matter"),
        limitations=("phase-only: does not control |S21xx|; amplitude mismatch is not represented",),
    )


def generalized_transfer_operator(ex_ref, ey_ref, ex_obj, ey_obj):
    """Candidate D (Sec. 5): the minimum-Frobenius-norm 2x2 operator T(x,y) solving the single
    local vector equation `E_obj = T . E_ref` at every point:

        T = outer(E_obj, conj(E_ref)) / (E_ref^H . E_ref)

    (the standard minimum-norm solution of an underdetermined linear system -- direct linear
    algebra, not a physical law). Reduces to a scalar ratio `Ex_obj/Ex_ref` when the reference has
    no y-component (this project's default synthetic reference fields) and Txy=Tyx=0 only if the
    object wave also has no y-component; a nonzero Txy/Tyx is exactly the *requirement* for
    polarization conversion (Sec. 5's off-diagonal terms), reported, never assumed zero."""
    ex_ref, ey_ref = np.asarray(ex_ref, dtype=complex), np.asarray(ey_ref, dtype=complex)
    ex_obj, ey_obj = np.asarray(ex_obj, dtype=complex), np.asarray(ey_obj, dtype=complex)
    denom = np.abs(ex_ref) ** 2 + np.abs(ey_ref) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        Txx = ex_obj * np.conj(ex_ref) / denom
        Txy = ex_obj * np.conj(ey_ref) / denom
        Tyx = ey_obj * np.conj(ex_ref) / denom
        Tyy = ey_obj * np.conj(ey_ref) / denom
    undefined = denom == 0.0   # E_ref vanishes: T is not defined there (never fabricated)
    values = {"Txx": Txx, "Txy": Txy, "Tyx": Tyx, "Tyy": Tyy, "undefined_mask": undefined}
    return ControlResponseResult(
        formulation="D", value_kind="transfer_operator", values=values,
        required_quantity="local 2x2 linear-basis transfer operator [[Sxx,Sxy],[Syx,Syy]]",
        mathematical_derivation="T = outer(E_obj, conj(E_ref)) / (E_ref^H . E_ref) (minimum-norm solution of E_obj = T.E_ref)",
        physical_meaning="local scattering operator relating incident (reference) to desired (object) field vectors",
        dimensionality="complex-valued 2x2 matrix, dimensionless",
        phase_convention="exp(+j*omega*t); exp(-j*k.r)",
        required_physical_surface_response="full linear-basis 2x2 S-block (Sxx, Sxy, Syx, Syy)",
        assumptions=(
            "the minimum-norm T is the physically intended one (not unique in general)",
            "a single (E_ref, E_obj) sample per point determines the whole 2x2 operator",
        ),
        limitations=(
            "not unique unless E_ref genuinely has two independent components at every point",
            "undefined wherever E_ref = 0 (undefined_mask)",
        ),
    )
