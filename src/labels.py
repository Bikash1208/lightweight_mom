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
