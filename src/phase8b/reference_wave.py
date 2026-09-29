"""Phase 8B Sec. 8-10, 43-44: the reference-wave definition.

This project has no monopole-feed implementation, no surface-wave/guided-mode solver, and no
dispersion-relation code anywhere in `src/` (grep-verified, `results/phase8a_feasibility/`'s own
audit and Sec. 6's repository audit below). A "monopole reference wave" is standard terminology in
the holographic-impedance-surface antenna literature this phase's structure matches (Fong, Colburn,
Ottusch, Visher, Sievenpiper, "Scalar and Tensor Holographic Artificial Impedance Surfaces," IEEE
Trans. Antennas Propag., vol. 58, no. 10, pp. 3212-3221, 2010 -- `source_registry.json` entry
`fong2010`): a monopole at the centre of a grounded modulated-impedance surface launches a
cylindrically-spreading **guided surface wave** along the surface, which serves as the reference
wave. Reproducing that field requires the guided surface wave's propagation constant
`k_swg` -- generally *not* `k0`, and only obtainable from a surface-wave dispersion analysis of the
periodic modulated boundary, which this project's solver does not perform (it solves only normal-
incidence plane-wave scattering, `CONVENTIONS.md` "Ports"). `monopole_reference_field` therefore
requires `k_swg` as an explicit, non-defaulted argument and documents, rather than hides, why: this
project cannot compute it.

Sec. 9's three input types are supported as explicit, distinguished cases (never silently
interchanged, Sec. 10):
"""

from dataclasses import dataclass
from enum import Enum

import numpy as np


class ReferenceFieldKind(Enum):
    COMPLEX_E_FIELD = "complex_E_field"          # Sec. 9.A: Ex(x,y), Ey(x,y), (Ez optional)
    PHASE_ONLY = "phase_only"                      # Sec. 9.B: phi_ref(x,y), no amplitude
    SURFACE_IMPEDANCE = "surface_impedance"        # Sec. 9.C: Zs(x,y)


@dataclass(frozen=True)
class ReferenceField:
    """One reference-field input, tagged with exactly which of Sec. 9's three types it is --
    the holographic-modulation step (`holographic_modulation.py`) branches on ``kind`` explicitly
    rather than assuming any two of these are interchangeable (Sec. 9-10)."""

    kind: ReferenceFieldKind
    values: np.ndarray            # complex E (A), real phase [deg] (B), or complex Zs [ohm] (C)
    x_m: np.ndarray
    y_m: np.ndarray

    def __post_init__(self):
        values = np.asarray(self.values)
        x_m, y_m = np.asarray(self.x_m, dtype=float), np.asarray(self.y_m, dtype=float)
        if values.shape != x_m.shape or values.shape != y_m.shape:
            raise ValueError("values/x_m/y_m must share one shape")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "x_m", x_m)
        object.__setattr__(self, "y_m", y_m)


def phase_only_reference_requires_amplitude_assumption():
    """Sec. 43: does this project's formulation support a phase-only reference field (no
    amplitude) without inventing one? Answer: **no** -- `holographic_modulation.interference_modulation`
    (Formulation A, the Fong et al. interference kernel `Re{E_ref . conj(E_obj)}`) is a product of
    two *complex* fields; using a phase-only reference there requires assuming a specific amplitude
    profile (uniform, or the standard 1/sqrt(rho) cylindrical-spreading law) that Sec. 43 forbids
    inventing silently. `holographic_modulation.phase_difference_hologram` (Formulation B, a pure
    phase-difference hologram) is the one case where a phase-only reference is sufficient by
    construction (it only ever uses arg(E_ref), never |E_ref|) -- documented per-formulation
    rather than as one project-wide answer, since the two candidate formulations differ here."""
    return {
        "formulation_A_fong2010_interference_kernel": "PHASE_ONLY_REFERENCE_UNSUPPORTED "
        "(the interference kernel needs |E_ref|; assuming an amplitude profile is not permitted "
        "without it being separately, explicitly justified)",
        "formulation_B_phase_difference_hologram": "phase-only reference is sufficient (the "
        "kernel only ever uses arg(E_ref))",
    }


def surface_impedance_reference_sufficiency():
    """Sec. 44: is a supplied `Zs(x,y)` alone sufficient? Only if the formulation's control
    variable *is* surface impedance (Formulation A) -- and even then, `Zs(x,y)` as a *reference*
    field (rather than a computed unit-cell *response*) is not a quantity either candidate
    formulation in this project actually asks for as an input; see
    `control_variable.py`/Sec. 24 for the separate question of whether the *unit-cell library*
    can *produce* a Zs quantity."""
    return "SURFACE_IMPEDANCE_INPUT_REQUIRES_ADDITIONAL_DATA (not a defined reference-wave input for either candidate formulation in this project; see control_variable.py)"


def monopole_reference_field(x_m, y_m, k_swg, amplitude_at_1m=1.0):
    """The standard cylindrically-spreading guided-surface-wave field of a monopole-fed
    holographic impedance surface (Fong et al. 2010; asymptotic cylindrical-wave form, standard in
    the surface-wave/leaky-wave antenna literature): `E_ref(rho) ~ rho^(-1/2) exp(-j k_swg rho)`
    (Sec. 19's plane-wave phase convention `exp(-j k.r)` applied radially). ``k_swg`` (the guided
    surface-wave propagation constant) has **no default**: this project's solver never computes a
    surface-wave dispersion relation (`phase8b.conventions.SURFACE_WAVE_DISPERSION_IMPLEMENTED`
    is `False`), so a caller must supply it from elsewhere (e.g. an independent eigenmode solve,
    out of this project's scope) -- this function will not silently substitute `k0` for it."""
    if k_swg is None:
        raise ValueError(
            "monopole_reference_field: k_swg (guided surface-wave propagation constant) is "
            "required and has no default -- this project's solver cannot compute it (no "
            "surface-wave/dispersion solver exists in src/); see reference_wave.py's module "
            "docstring"
        )
    x_m, y_m = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    rho = np.sqrt(x_m**2 + y_m**2)
    with np.errstate(divide="ignore"):
        amplitude = amplitude_at_1m / np.sqrt(np.maximum(rho, np.finfo(float).eps))
    return amplitude * np.exp(-1j * k_swg * rho)
