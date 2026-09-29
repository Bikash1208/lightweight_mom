"""Phase 8B Sec. 21-25: the physical unit-cell control-variable investigation.

Every conclusion below is one of `VALIDATED` / `INVALID` / `INSUFFICIENT_EVIDENCE` -- a forced
binary answer is not used when the evidence does not support one (Sec. 22's explicit instruction).
Every numeric claim is either read from this project's own files (cited by path) or a standard,
cited piece of EM theory (`source_registry.json`), never invented.
"""

from dataclasses import dataclass

CASE01_S11XX_DB_RANGE = (-15.89, -12.79)  # results/validation_cases/case01_non_resonant_patch/mom_results.json
CASE01_S21XX_DB_RANGE = (-0.24, -0.12)  # same file, recomputed directly from stored |S| (see
# tests/test_phase8b_control_variable.py, which re-derives
# these two ranges from the file itself rather than trusting
# the literal numbers here)


@dataclass(frozen=True)
class ControlVariable:
    name: str
    physical_definition: str
    complex_or_real: str  # "complex" | "real"
    phase_relevant: bool
    magnitude_relevant: bool
    library_quantity: str | None  # the Phase 7 database field that would supply it, or None
    validated: str  # "VALIDATED" | "INVALID" | "INSUFFICIENT_EVIDENCE"
    justification: str


def investigate_s11xx():
    """Sec. 22: S11xx = the co-polarized (x-in, x-out) complex reflection coefficient at port 1
    (the top reference plane), for an x-polarized plane wave normally incident from +z
    (`CONVENTIONS.md` "Ports", `S_PARAMETERS.md`). Answers Sec. 22's 11 questions in the
    justification text below."""
    justification = (
        "(1) Physically: the complex ratio of the reflected to incident tangential-E wave "
        "amplitude at port 1 (top reference plane), co-polarized x-to-x (`S_PARAMETERS.md`). "
        "(2) Incident field: a normally-incident (theta=0), x-polarized plane wave, unit "
        "amplitude referenced at the pattern plane z=0 (`excitation.py`). (3) Reflected field: "
        "the x-component of the scattered wave travelling +z, referenced at z_top "
        "(CONVENTIONS.md 'Ports'). (4) Reference plane: z_top = +top_air_thickness, in air, above "
        "the pattern (not at the metal itself). (5) Polarization: linear x, co-polarized. "
        "(6) The unit cell is NOT grounded (README.md 'Assumptions' item 2, unchanged since "
        "Phase 1: air | pattern | substrate | air, no ground plane). (7) Whether the intended MTS "
        "effect is reflective or transmissive is an architecture decision this project has never "
        "made (`UNIT_CELL_LIBRARY.md` Sec. 4, `phase_extraction.py`); Phase 6's case01 measured "
        f"the structure as *mostly transmissive*: |S11xx| in {CASE01_S11XX_DB_RANGE} dB vs. "
        f"|S21xx| in {CASE01_S21XX_DB_RANGE} dB (recomputed directly from "
        "`results/validation_cases/case01_non_resonant_patch/mom_results.json`), the opposite of "
        "a typical grounded reflectarray element where |S11|~=0 dB (near-total reflection). "
        "(8) Whether the holographic modulation is represented by a reflection coefficient at "
        "all depends on which candidate formulation is used: Formulation B "
        "(`holographic_modulation.phase_difference_hologram`) can use *either* S11xx or S21xx "
        "phase as its control variable in principle (both are normal-incidence co-pol phases of "
        "the right *kind*); Formulation A (Fong et al. 2010) needs a **surface reactance** of a "
        "**grounded** structure, which S11xx of this ungrounded, partially-transmissive cell is "
        "not (see `investigate_surface_impedance`). (9) S11xx phase alone is not obviously "
        "sufficient without a magnitude/usability floor -- Phase 7 already applies one "
        "(`library_schema.PHASE_USABLE_MIN_MAG_DB`). (10) Only Formulation B's kernel could use "
        "reflection-coefficient phase directly; Formulation A cannot (structural reason above). "
        "(11) Phase 7 stores S11xx as `arg(S)` on (-180, 180], a fixed, consistent reference "
        "(`phase.py`, `database.CONVENTIONS`)."
    )
    return ControlVariable(
        name="S11_xx",
        physical_definition="complex co-polarized (x-in, x-out) reflection coefficient at the "
        "top reference plane, normal incidence, x-polarized",
        complex_or_real="complex",
        phase_relevant=True,
        magnitude_relevant=True,
        library_quantity="S11 (family 'S11', component 'xx')",
        validated="INSUFFICIENT_EVIDENCE",
        justification=justification
        + " Conclusion: S11xx is a well-defined, correctly-stored physical quantity, but whether "
        "it is *the* intended holographic control variable depends on an architecture decision "
        "(reflective-vs-transmissive; which candidate formulation) this project has not made -- "
        "neither VALIDATED nor INVALID can be honestly asserted; see "
        "`library_compatibility.py` for the per-formulation resolution this ambiguity allows.",
    )


def investigate_s21xx():
    """Sec. 23: S21xx = the co-polarized complex transmission coefficient (port 1 -> port 2)."""
    justification = (
        "S21xx is the complex ratio of the transmitted (port 2, travelling -z) to incident "
        "(port 1, travelling -z) tangential-E wave amplitude, co-polarized x-to-x. This is NOT a "
        "grounded structure (no ground plane), so unlike a typical reflectarray/AMC unit cell, "
        f"S21xx is *not* negligible here: |S21xx| in {CASE01_S21XX_DB_RANGE} dB (near-total "
        "transmission) vs. |S11xx| in "
        f"{CASE01_S11XX_DB_RANGE} dB for Phase 6's case01 (same source as `investigate_s11xx`). "
        "This makes S21xx phase, not S11xx phase, the *larger*-magnitude, better-conditioned "
        "response for this specific ungrounded geometry family -- but a larger magnitude alone "
        "is not sufficient justification to prefer it (Sec. 23's explicit prohibition: 'do not "
        "select S21 merely because it offers a larger phase span/better conditioning'). Whether "
        "S21xx is the physically meaningful quantity depends on whether the intended MTS "
        "architecture is transmissive (a transmit-array-style hologram, e.g. the object wave "
        "radiates from the side opposite the feed) or reflective -- an architecture decision this "
        "project has never made."
    )
    return ControlVariable(
        name="S21_xx",
        physical_definition="complex co-polarized (x-in, x-out) transmission coefficient, "
        "port 1 to port 2, normal incidence, x-polarized",
        complex_or_real="complex",
        phase_relevant=True,
        magnitude_relevant=True,
        library_quantity="S21 (family 'S21', component 'xx')",
        validated="INSUFFICIENT_EVIDENCE",
        justification=justification,
    )


def investigate_surface_impedance():
    """Sec. 24: can this project's solver produce a surface impedance (or equivalent) quantity?"""
    justification = (
        "Formulation A (Fong et al. 2010) needs the surface reactance X(x,y) of a **grounded** "
        "modulated-impedance boundary, most commonly obtained in the literature from the "
        "closed-form normal-incidence relation Zs = eta0 (1+Gamma)/(1-Gamma) (standard "
        "transmission-line/impedance-boundary reflection-coefficient relation, "
        "`source_registry.json` `impedance_boundary_relation`), where Gamma is the reflection "
        "coefficient of a **grounded** unit cell (|Gamma|=1 for a lossless grounded boundary, so "
        "Zs is purely reactive). Phase 1-7's solver stack has **no ground plane** "
        "(`phase8b.conventions.NO_GROUND_PLANE`) and Phase 6's case01 shows the structure is "
        f"mostly transmissive (|S11xx| in {CASE01_S11XX_DB_RANGE} dB, far from |Gamma|=1=0 dB). "
        "Applying Zs=eta0(1+S11xx)/(1-S11xx) to this data anyway would not fail numerically, but "
        "it would not represent a physical surface impedance: the formula's derivation assumes "
        "total reflection from a boundary condition of the form E_t = Zs (n x H_t), which this "
        "partially-transmissive, ungrounded two-port structure does not satisfy (some power "
        "genuinely leaves through the back, which the single-boundary impedance model has no way "
        "to represent -- it would appear as a large, physically meaningless 'loss' term in a "
        "structure whose actual dielectric loss is tan(delta)=0.0009, negligible). Separately, "
        "the formulation's reference wave (`reference_wave.monopole_reference_field`) is itself "
        "a *guided surface wave*, whose propagation constant k_swg this project's solver cannot "
        "compute (`phase8b.conventions.SURFACE_WAVE_DISPERSION_IMPLEMENTED = False`: no "
        "eigenmode/dispersion solver exists in `src/`)."
    )
    return ControlVariable(
        name="surface_impedance_Zs",
        physical_definition="surface reactance/impedance boundary condition of a grounded, "
        "surface-wave-supporting periodic structure",
        complex_or_real="complex",
        phase_relevant=True,
        magnitude_relevant=True,
        library_quantity=None,
        validated="INVALID",
        justification=justification
        + " Conclusion: SURFACE_IMPEDANCE_NOT_AVAILABLE -- two independent, structural reasons "
        "(no ground plane invalidates the standard Zs<->Gamma relation; no surface-wave "
        "dispersion solver exists to characterize the reference wave itself), neither of which "
        "is a missing feature that could be worked around without changing the frozen solver.",
    )


def investigate_reflection_coefficient():
    """Sec. 25: is S11xx actually 'the' reflection coefficient Gamma a holographic formulation
    needing Gamma(x,y) would want?"""
    return ControlVariable(
        name="reflection_coefficient_Gamma",
        physical_definition="complex reflection coefficient at the top reference plane",
        complex_or_real="complex",
        phase_relevant=True,
        magnitude_relevant=True,
        library_quantity="S11 (family 'S11', component 'xx')",
        validated="VALIDATED",
        justification=(
            "S11xx *is* mathematically the reflection coefficient Gamma_xx at port 1 by "
            "construction (`S_PARAMETERS.md` Sec. 1: b1 = S11 a1 + S12 a2; with only port 1 "
            "excited, S11 = b1/a1 = Gamma). This equivalence (Gamma_xx = S11_xx) is physically "
            "justified and requires no additional assumption -- unlike the *surface-impedance* "
            "interpretation of that same Gamma (`investigate_surface_impedance`, which is INVALID "
            "for this ungrounded structure), 'S11xx is a reflection coefficient' is simply what "
            "S11xx already, definitionally, is (`database.py` module docstring, `S_PARAMETERS.md`)."
        ),
    )


ALL_CANDIDATES = (
    investigate_s11xx,
    investigate_s21xx,
    investigate_surface_impedance,
    investigate_reflection_coefficient,
)
