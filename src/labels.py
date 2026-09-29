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
HFSS_VALIDATED_OFF_RESONANCE = "HFSS_VALIDATED_OFF_RESONANCE"          # e.g. G1: ~0.3deg mean phase diff
MESH_CONVERGED = "MESH_CONVERGED"                                       # Richardson fit succeeded, within tol
MESH_NOT_CONVERGED = "MESH_NOT_CONVERGED"                               # fixed mesh only, or fit failed/non-monotonic
LOCAL_PERIODICITY_APPROXIMATION = "LOCAL_PERIODICITY_APPROXIMATION"     # lookup-table forward model, not full-wave
ZERO_THICKNESS_METAL_MODEL = "ZERO_THICKNESS_METAL_MODEL"               # solver's PEC pattern has no thickness
EXPERIMENTAL = "EXPERIMENTAL"                                            # not yet validated against anything
DISPERSION_FROM_HFSS = "DISPERSION_FROM_HFSS"                            # dispersion table imported from HFSS
DISPERSION_FROM_MOM = "DISPERSION_FROM_MOM"                              # (scaffold only -- see src/dispersion/)
REFERENCE_AMPLITUDE_ASSUMED_UNIFORM = "REFERENCE_AMPLITUDE_ASSUMED_UNIFORM"
FEED_REGION_MASKED = "FEED_REGION_MASKED"
GEOMETRY_INVERSION_LOCKED = "GEOMETRY_INVERSION_LOCKED"
STEP_EXPORT_LOCKED = "STEP_EXPORT_LOCKED"


def hfss_residual_at_resonance(residual_deg):
    """G4's own measured ~3-4deg unexplained residual after the mesh/thickness budget
    (PHASE12D_ADDENDUM.md, `src/phase13/hfss_resonance_ingest.py`)."""
    return f"HFSS_RESIDUAL_AT_RESONANCE(~{residual_deg:.1f}deg)"


def phase_coverage_shortfall(achieved_span_deg, target_span_deg=360.0):
    shortfall = target_span_deg - achieved_span_deg
    return f"PHASE_COVERAGE_SHORTFALL({shortfall:.0f}deg)"


def steep_gradient(threshold_deg_per_cell):
    """A per-cell phase gradient exceeds the measured local-periodicity step limit (12D:
    5.8 deg/cell off-resonance)."""
    return f"STEEP_GRADIENT(>{threshold_deg_per_cell:.1f}deg/cell)"


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
