"""Phase 13 Part 2: the fixed validation-label vocabulary (PHASE13_WEBAPP_PROMPT.md Part 2) --
every result object the app returns carries a tuple of these, and the UI shows them next to the
result. Labels are grounded in this project's own recorded evidence (docs/phases/, `gate.py`
modules of earlier phases); nothing here invents a claim a phase didn't actually measure.

Some labels are plain (a fixed string); others are parameterized (a value baked into the label at
call time, e.g. `phase_coverage_shortfall(174.0)` -> `"PHASE_COVERAGE_SHORTFALL(174.0deg)"`) so the
UI can show the actual number without a separate lookup.
"""

from dataclasses import dataclass

# ---------------------------------------------------------------- plain (unparameterized) labels
HFSS_VALIDATED_OFF_RESONANCE = "HFSS_VALIDATED_OFF_RESONANCE"  # e.g. G1: ~0.3deg mean phase diff
MESH_CONVERGED = "MESH_CONVERGED"  # Richardson fit succeeded, within tol
MESH_NOT_CONVERGED = "MESH_NOT_CONVERGED"  # fixed mesh only, or fit failed/non-monotonic
LOCAL_PERIODICITY_APPROXIMATION = "LOCAL_PERIODICITY_APPROXIMATION"  # lookup-table forward model, not full-wave
ZERO_THICKNESS_METAL_MODEL = "ZERO_THICKNESS_METAL_MODEL"  # solver's PEC pattern has no thickness
EXPERIMENTAL = "EXPERIMENTAL"  # not yet validated against anything
DISPERSION_FROM_HFSS = "DISPERSION_FROM_HFSS"  # dispersion table imported from HFSS
DISPERSION_FROM_MOM = "DISPERSION_FROM_MOM"  # dispersion table from this project's own MoM eigenmode search (src/dispersion/mom_scaffold.py) -- a research item, see that module's own docstring
DISPERSION_FROM_SHEET_MODEL = "DISPERSION_FROM_SHEET_MODEL"  # homogenized patch-sheet-admittance + transverse resonance (src/dispersion/sheet_model.py), Luukkonen et al. 2008 style -- an approximation (angle-independent sheet admittance), not a rigorous eigenmode
REFERENCE_AMPLITUDE_ASSUMED_UNIFORM = "REFERENCE_AMPLITUDE_ASSUMED_UNIFORM"
FEED_REGION_MASKED = "FEED_REGION_MASKED"
GEOMETRY_INVERSION_LOCKED = "GEOMETRY_INVERSION_LOCKED"
STEP_EXPORT_LOCKED = "STEP_EXPORT_LOCKED"
VECTOR_PATTERN_MODEL = "VECTOR_PATTERN_MODEL"  # PHASE15 Part B: far field via the layered Green's
# function (reciprocity), not the scalar free-space Fourier transform
SCALAR_PATTERN_MODEL_LEGACY = "SCALAR_PATTERN_MODEL_LEGACY"  # PHASE15 Part B: the pre-Phase-15
# scalar model, kept for comparison/regression -- no polarization/vortex information
MULTI_PORT_SHARED_METHOD = "MULTI_PORT_SHARED_METHOD"  # PHASE15 Part C: every port of one
# multi-port design shares the same formulation and dispersion table / phase library -- a
# documented scope reduction, not a genuinely independent per-port method choice
PORT_ISOLATION_NOT_FULL_WAVE = "PORT_ISOLATION_NOT_FULL_WAVE"  # PHASE15 Part C: the inter-port
# isolation estimate is a first-order edge-distance-vs-surface-wavelength heuristic, not a
# full-wave (HFSS) S-parameter result
CODED_DESIGN_METHOD = "CODED_DESIGN_METHOD"  # PHASE15 Part E: 1-bit/N-bit phase quantization
# (Mishra TVT 2026), rather than a continuous reflection-phase or impedance-modulation hologram
ANNULAR_SLOT_PATCH_GEOMETRY = "ANNULAR_SLOT_PATCH_GEOMETRY"  # PHASE15 Part E: the second unit-cell
# geometry type (a square/rectangular patch with a concentric annular slot), solved by this
# project's own MoM the same way the square patch is


def hfss_residual_at_resonance(residual_deg):
    """G4's own measured ~3-4deg unexplained residual after the mesh/thickness budget
    (PHASE12D_ADDENDUM.md, `src/phase13/hfss_resonance_ingest.py`)."""
    return f"HFSS_RESIDUAL_AT_RESONANCE(~{residual_deg:.1f}deg)"


def phase_coverage_shortfall(worst_error_deg):
    """`worst_error_deg`: the largest per-cell achieved-vs-required phase error among the cells a
    design flagged as not reachable (PHASE13B A4 -- reused by `reflectarray.designer` and
    `hologram.hmsa_synthesis`, both against `PHASE_ERROR_UNREACHABLE_THRESHOLD_DEG`)."""
    return f"PHASE_COVERAGE_SHORTFALL({worst_error_deg:.0f}deg)"


def steep_gradient(threshold_deg_per_cell):
    """A per-cell phase gradient exceeds the measured local-periodicity step limit (12D:
    5.8 deg/cell off-resonance)."""
    return f"STEEP_GRADIENT(>{threshold_deg_per_cell:.1f}deg/cell)"


def reference_wave_rescaled_to_beta_bar(beta_bar_rad_per_m):
    """PHASE13B B3: the imported bare-slab reference phase was rescaled to the guided wavenumber
    beta_bar implied by the hologram's own mean surface reactance (dispersion library), instead of
    being used as-is."""
    return f"REFERENCE_WAVE_RESCALED_TO_BETA_BAR({beta_bar_rad_per_m:.1f}rad/m)"


def reference_wave_not_rescaled(worst_case_error_deg):
    """PHASE13B B3: the raw bare-slab reference phase was used directly; `worst_case_error_deg`
    is the accumulated phase error ~ |beta_bar - beta_slab|*r_max this could cost at the aperture
    edge."""
    return f"REFERENCE_WAVE_NOT_RESCALED(~{worst_case_error_deg:.0f}deg accumulated error possible)"


def isolation_edge_distance(port_a, port_b, gap_lambda_sw):
    """PHASE15 Part C: two ports' sub-apertures sit closer than the recommended guard-row
    distance (`hologram.ports.isolation_estimate`'s `min_lambda_sw`), reported in units of the
    board's own surface-wave wavelength."""
    return f"PORT_EDGE_DISTANCE_TIGHT(port{port_a}-port{port_b},{gap_lambda_sw:.2f}lambda_sw)"


def quantization_loss_db(loss_db, n_bits, textbook_db):
    """PHASE15 Part E: main-beam PEAK power lost by quantizing the continuous coded-method phase to
    `n_bits` bits, at the same aperture power, measured on this project's own vector model. The
    textbook value for uniform phase error (`coded_synthesis.textbook_quantization_loss_db`) is
    carried alongside so the two can be compared, never one substituted for the other."""
    return f"QUANTIZATION_LOSS({loss_db:.2f}dB,{n_bits}bit,textbook {textbook_db:.2f}dB)"


def two_state_surface_wave_check(beta_bar_rad_per_m, spread_fraction, spread_limit):
    """PHASE15 Part E: the pair's own TM0 surface-wave check -- the mean beta the coded phase uses,
    and how far the two states' betas disagree relative to the allowed limit."""
    return (
        f"TWO_STATE_SURFACE_WAVE(beta_bar={beta_bar_rad_per_m:.1f}rad/m,"
        f"spread={100.0 * spread_fraction:.2f}%,limit={100.0 * spread_limit:.0f}%)"
    )


def two_state_phase_separation(achieved_deg, n_bits=1):
    """PHASE15 Part E: the achieved reflection-phase separation between this design's two (or
    2^n_bits) coded states, as actually measured from the chosen cell pair/library -- reported
    honestly even when far from the ideal 360/2^n_bits deg spacing (never assumed)."""
    return f"TWO_STATE_PHASE_SEPARATION({achieved_deg:.1f}deg,{n_bits}bit)"


def feed_position_auto_detected(x_mm, y_mm):
    """PHASE15 Part A: the feed position used for polarization compensation was taken from the
    uploaded reference field's own NaN/singularity location, not entered by the user."""
    return f"FEED_POSITION_AUTO_DETECTED({x_mm:.1f}mm,{y_mm:.1f}mm)"


def reference_field_frequency_mismatch(upload_frequency_hz, synthesis_frequency_hz):
    """PHASE13B A3: the frequency the reference field was uploaded/analyzed at (e.g. its TM0
    check) differs from the frequency synthesis actually ran at -- silently accepting this
    produced a 100%-unreachable design with no explanation in the live review."""
    return (
        f"REFERENCE_FIELD_FREQUENCY_MISMATCH(upload={upload_frequency_hz / 1e9:g}GHz,"
        f"synthesis={synthesis_frequency_hz / 1e9:g}GHz)"
    )


@dataclass(frozen=True)
class LabeledResult:
    """A thin wrapper any app service can return: the underlying result value plus its labels
    and any warnings (plain strings, human-readable -- distinct from the fixed label vocabulary
    above, which is meant to be matched on / filtered by the UI)."""

    value: object
    labels: tuple = ()
    warnings: tuple = ()

    def with_label(self, label):
        return LabeledResult(self.value, self.labels + (label,), self.warnings)

    def with_warning(self, warning):
        return LabeledResult(self.value, self.labels, self.warnings + (warning,))

    def has(self, label):
        return label in self.labels or any(str(l).startswith(label.split("(")[0]) for l in self.labels)
