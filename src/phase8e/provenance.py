"""Sec. 27: formulation provenance for every candidate (A-D). Candidates A and B/C's provenance
was already established in Phase 8B (`source_registry.json`, `fong2010`); this module restates it
(never re-derives it) and adds C/D, which Phase 8B did not need to distinguish.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class FormulationProvenance:
    formulation_name: str
    mathematical_reference: str
    assumptions: tuple
    architecture_requirements: tuple
    applicability: str
    incompatibilities: tuple
    implementation_status: str


def build_provenance_registry():
    return {
        "A": FormulationProvenance(
            formulation_name="A: interference-modulation index (grounded holographic impedance surface)",
            mathematical_reference=(
                "Fong, Colburn, Ottusch, Visher, Sievenpiper, 'Scalar and Tensor Holographic "
                "Artificial Impedance Surfaces,' IEEE Trans. Antennas Propag., vol. 58, no. 10, "
                "pp. 3212-3221, 2010 (source_registry.json 'fong2010', already cited in Phase 8B); "
                "the standard reflection-coefficient/impedance-boundary relation "
                "Zs = eta0(1+Gamma)/(1-Gamma) (textbook, 'impedance_boundary_relation')."
            ),
            assumptions=(
                "grounded modulated-impedance boundary (|Gamma|=1, purely reactive Zs)",
                "monopole-launched cylindrical guided surface wave as the reference field",
                "surface-wave propagation constant k_swg known from a dispersion analysis",
            ),
            architecture_requirements=("ground plane", "surface-wave-supporting periodic boundary"),
            applicability="grounded modulated-impedance leaky-wave holographic antennas",
            incompatibilities=(
                "this project's solver has no ground plane (phase8b.conventions.NO_GROUND_PLANE)",
                "this project's solver has no surface-wave dispersion solver "
                "(phase8b.conventions.SURFACE_WAVE_DISPERSION_IMPLEMENTED = False)",
            ),
            implementation_status="MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED "
            "(phase8b.holographic_modulation), ARCHITECTURE INCOMPATIBLE",
        ),
        "B": FormulationProvenance(
            formulation_name="B: phase-difference hologram, reflection correspondence",
            mathematical_reference=(
                "Standard reflectarray phase-hologram design equation (phase-conjugate/phase-"
                "difference hologram): required_phase = wrap(arg(E_obj) - arg(E_ref)); textbook "
                "reflectarray antenna design (e.g. Huang & Encinar, 'Reflectarray Antennas,' "
                "Wiley-IEEE Press, 2007 -- standard reference for phase-only reflectarray unit-"
                "cell synthesis; already implemented, unmodified, as "
                "phase8b.holographic_modulation.phase_difference_hologram)."
            ),
            assumptions=(
                "each unit cell's normal-incidence reflection phase alone reproduces the locally "
                "required phase shift (local periodic approximation)",
                "the reflected field is the physically radiating one (reflective architecture)",
            ),
            architecture_requirements=("normal-incidence co-polarized reflection response (S11)",),
            applicability="passive reflectarray-style phase holograms",
            incompatibilities=(),
            implementation_status="MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED, "
            "ARCHITECTURE COMPATIBILITY UNVALIDATED (Phase 8B/8C)",
        ),
        "C": FormulationProvenance(
            formulation_name="C: phase-difference hologram, transmission correspondence",
            mathematical_reference=(
                "The identical kernel as B (phase8b.holographic_modulation.phase_difference_hologram "
                "uses only arg(E_ref)/arg(E_obj), independent of which physical wave carries the "
                "radiated field), applied to a transmitarray-style architecture (e.g. Huang & "
                "Encinar, op. cit., Ch. 8, 'Transmitarray Antennas'): the reference wave transmits "
                "through the surface and the far side carries the shaped object wave."
            ),
            assumptions=(
                "each unit cell's normal-incidence transmission phase alone reproduces the locally "
                "required phase shift (local periodic approximation)",
                "the transmitted field is the physically radiating one (transmissive architecture)",
            ),
            architecture_requirements=("normal-incidence co-polarized transmission response (S21)",),
            applicability="passive transmitarray-style phase holograms",
            incompatibilities=(),
            implementation_status="MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED, "
            "ARCHITECTURE COMPATIBILITY UNVALIDATED (Phase 8B/8C)",
        ),
        "D": FormulationProvenance(
            formulation_name="D: generalized (vector/matrix) local transfer operator",
            mathematical_reference=(
                "Minimum-Frobenius-norm solution of the underdetermined local relation "
                "E_obj(x,y) = T(x,y) . E_ref(x,y) for a 2x2 (or scalar) T: "
                "T = outer(E_obj, conj(E_ref)) / |E_ref|^2 (direct linear-algebra derivation, "
                "documented in control_response.py; the same tensor-surface generalization Fong "
                "et al. 2010's title already names ('Scalar AND TENSOR Holographic Artificial "
                "Impedance Surfaces'), applied here to a scattering operator rather than a "
                "grounded impedance tensor)."
            ),
            assumptions=(
                "the required local transfer operator is exactly the minimum-norm one consistent "
                "with the single observed (E_ref, E_obj) pair at each point -- NOT unique unless "
                "E_ref is a genuine 2-component vector at every point",
                "the unit cell's own 2x2 linear-basis S-block (Sxx, Sxy, Syx, Syy) is the "
                "measurable quantity that would need to realize T(x,y)",
            ),
            architecture_requirements=("polarization-converting (non-diagonal) scattering response",),
            applicability="polarization-converting / tensor holographic surfaces",
            incompatibilities=(
                "Phase 8C measured negligible cross-polarization (<1e-6 of incident power) across "
                "every one of the 324 cells in the only registered library "
                "(phase8c.power_balance.library_power_summary) -- this specific unit-cell family "
                "cannot realize a non-diagonal T(x,y), independent of whether the formulation "
                "itself is otherwise valid",
            ),
            implementation_status="MATHEMATICALLY DEFINED, NUMERICALLY IMPLEMENTED "
            "(control_response.generalized_transfer_operator), LIBRARY DATA INSUFFICIENT",
        ),
    }
