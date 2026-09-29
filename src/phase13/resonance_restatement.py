"""Part 0 task 2 (PHASE12D_ADDENDUM.md): re-state the Phase 12B/12C/12D 8 mm results' own
resonance-region mesh error, without editing any of their outputs.

Three different meshes were actually used for the 8 mm grounded design point across phases:
  - the design-point CONFIG's own single-cell mesh (`config_phase12c_grounded_8mm.yaml`):
    40 cells, dx = 0.2 mm (`phase12c.design_point.BASE_MESH_CELLS`) -- the addendum's own "0.2 mm
    cells" reference;
  - `phase12c.locality_study.MESH_CELLS_PER_SUBCELL`: 20 cells, dx = 0.4 mm;
  - `phase12d.step_curve.EMBEDDED_MESH_CELLS_PER_SUBCELL`: 24 cells, dx = 0.333 mm.
Each is recomputed from the actual live constant, never hard-coded, so this module cannot silently
drift out of sync with those phases' real meshes.

**Method.** The addendum's own G4 Richardson fit (`phase13.convergence`, reproduced and checked
against the addendum's numbers in `tests/test_phase13_convergence.py`) gives a real, measured
power law `phi(dx) = phi_inf + a*dx^p` at resonance (p ~= 0.96). Re-solving a full 7-point
Richardson family at every mesh above (for every phase's own resonant L) would cost several more
hours of embedded-cell-equivalent compute for a side investigation the brief itself only asks to
"quantify" -- so, following the same cost-honesty discipline as Phase 12D's own step-size study,
this module instead SCALES that one already-measured, already-validated power law to each mesh's
own dx via `a*dx^p` (justified because the underlying cause -- an unresolved edge-current
singularity at the patch corners -- is a local mesh-resolution effect, not something that depends
on which nearby L was used to measure it). This is explicitly labelled a **scaled estimate**, not
an independent re-fit, everywhere it is used.

The step-size (embedded-vs-isolated) cancellation question is reported as genuinely UNRESOLVED:
answering it rigorously needs an embedded-supercell Richardson study, and Phase 12D's own cost
measurements already put a single embedded solve at the near/at-resonance sizes at this project's
practical ceiling even at ONE mesh resolution -- a multi-resolution embedded sweep is not attempted
here, and this module says so rather than guessing.
"""

from dataclasses import dataclass

from ..phase12c import design_point as dp
from ..phase12c.locality_study import MESH_CELLS_PER_SUBCELL as PHASE12C_LOCALITY_MESH_CELLS
from ..phase12d.step_curve import EMBEDDED_MESH_CELLS_PER_SUBCELL as PHASE12D_STEP_MESH_CELLS
from . import convergence as conv

RESONANT_L_MM = 6.0          # G4 -- the L value Phases 12C/12D/13 all treat as "at resonance"
RESONANT_FREQUENCY_HZ = 12.5e9

MESHES_TO_RESTATE = {
    "phase12c_design_point_config (single-cell library mesh)": dp.BASE_MESH_CELLS,
    "phase12c.locality_study (3x3/5x5 LPA suite mesh)": PHASE12C_LOCALITY_MESH_CELLS,
    "phase12d.step_curve (Task A step-size study mesh)": PHASE12D_STEP_MESH_CELLS,
}


@dataclass(frozen=True)
class RestatedMeshError:
    label: str
    cells: int
    dx_mm: float
    scaled_correction_deg: float   # |phi_inf - phi(this dx)|, scaled from the measured G4 power law


def restate_resonance_region_error(reference_fit=None):
    """Scales the measured G4 power law to every mesh in `MESHES_TO_RESTATE`. `reference_fit`
    lets a caller pass in an already-computed `phase13.convergence.FitResult` (from
    `hfss_resonance_ingest.converged_mom_curve`'s own G4-at-12.5GHz point, or a fresh
    `richardson_extrapolate_s11xx` call) to avoid re-solving; if omitted, a fresh 7-point
    Richardson family is run (the expensive path, only for a standalone call)."""
    if reference_fit is None:
        unit_cell = dp.unit_cell_at(RESONANT_L_MM * 1e-3)
        result = conv.richardson_extrapolate_s11xx(unit_cell, RESONANT_FREQUENCY_HZ, dp.REFERENCE_PLANES)
        reference_fit = result.phase
    if reference_fit.status != conv.STATUS_CONVERGED or reference_fit.coefficient is None:
        raise ValueError("restate_resonance_region_error needs a CONVERGED reference fit with a coefficient")

    period_mm = dp.PERIOD_M * 1e3
    out = []
    for label, cells in MESHES_TO_RESTATE.items():
        dx_mm = period_mm / cells
        dx_m = dx_mm * 1e-3
        correction_deg = abs(reference_fit.coefficient) * dx_m ** reference_fit.order
        out.append(RestatedMeshError(label, cells, dx_mm, correction_deg))
    return tuple(out)


STEP_SIZE_CANCELLATION_FINDING = (
    "UNRESOLVED: whether Phase 12C/12D's embedded-vs-isolated step-size errors (measured on the "
    "SAME mesh for both) partly cancel the mesh-convergence error found here requires an embedded-"
    "supercell Richardson study (multiple mesh resolutions of a 3x3+ supercell at near/at-resonance "
    "L values). Phase 12D's own cost measurements already place a SINGLE embedded solve at those "
    "sizes at this project's practical ceiling (~7-14 min real CPU time per case at one mesh); a "
    "multi-resolution sweep of the same case was not attempted in this phase. Reported as an open "
    "question, not answered by assumption."
)
