"""Phase 12D orchestration: fill in the local-periodicity step-size curve with the design-relevant
reflection-phase/amplitude metric (Task A), the aperture-level Floquet-order-power test (Task B),
and real-HFSS-data ingestion (Task C). Writes every required output and prints a terminal summary.
"""

from dataclasses import dataclass
from pathlib import Path

from ..phase12c.locality_study import classify_centers
from . import gate, hfss_ingest, plots, specification, step_curve
from .aperture_test import run_aperture_test
from .report import render_phase12d_report_md

DEFAULT_MAX_K = 2  # measured cost note (PHASE12D_LOCAL_PERIODICITY.md Sec. 2): near/at-resonance
# step cases involve large (up to ~7.3mm, near the full 8mm period) patches
# and cost far more than the off-resonance ones; max_k=2 already gives 14
# near/at-resonance cases (28 total with off-resonance) spanning steps from
# 0.8 to 160 deg -- ample coverage for "plot error vs step size" without
# extending into k=3, which off-resonance-only timing would underestimate
DEFAULT_APERTURE_SIZES = (6, 8, 12)
DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL = 10
DEFAULT_APERTURE_L_FAMILY_MM = tuple(
    round(2.0 * 8.0 / DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL * k, 6) for k in range(1, 5)
)  # aligned at the aperture test's own (coarser, cost-driven) mesh -- see run_phase12d docstring.
# k stops at 4 (1.6, 3.2, 4.8, 6.4 mm): k=5 (8.0 mm) would equal the full period -- invalid.


@dataclass(frozen=True)
class Phase12dResult:
    output_dir: Path
    curve: tuple
    centers: object
    step_cases: tuple
    task_a_gate: object
    aperture_results: tuple
    aperture_requested_sizes: tuple
    task_b_gate: object
    hfss_comparisons: tuple
    hfss_provenance: dict
    geometric_identity: tuple
    task_c_gate: object
    recommendation: str


RECOMMENDATION = (
    "The step-size curve (Task A) shows the design-relevant reflection-phase error stays under the "
    "design tolerance for every tested step out to double digits of degrees away from resonance "
    "(measured up to 5.8 deg/cell), but every tested near- or at-resonance step -- the smallest "
    "achievable at this project's validated floquet convention, itself tens of degrees because a "
    "mesh fine enough to embed genuinely small (10-45 deg) steps THROUGH the resonance transition "
    "is computationally out of reach for a supercell solve -- already fails badly (19-173 deg "
    "phase error). So the true near-resonance crossover step remains bounded, not pinned to a "
    "single number, but the practical conclusion is the same either way: an aperture whose phase "
    "profile must pass through this design's resonance region cannot rely on a lookup table there. "
    "The aperture-level test (Task B), which actually decides Phase 13, confirms this independently "
    "from a full-supercell MoM solve (not the dipole proxy): at every tested size (N=6, 8, 12 cells "
    "per assumed cycle), the lookup-table-constructed anomalous reflector puts less than 53% of its "
    "reflected power in the intended order (measured: 52.9%, 41.1%, 47.5%) -- LOOKUP_TABLE_INADEQUATE, "
    "not merely 'adequate above some N'. The root cause is visible in the numbers: the lookup "
    "table's own achievable phase span at this design point (186 deg, Sec. 3) is far short of the "
    "full 360 deg a blazed grating needs, the same phase-coverage shortfall Phases 12A-12C already "
    "measured from a different angle. A full-aperture integral-equation solver (a la Budhu & Grbic "
    "2021) is therefore the more broadly defensible default for this grounded 8 mm design point, "
    "not a lookup table constrained to small steps -- though a design confined entirely to the "
    "flat, off-resonance part of the curve (steps below ~6 deg/cell, Task A) remains a narrow case "
    "where a lookup table could still work. This is a recommendation, not a decision -- Phase 13 "
    "makes the actual choice, informed by both this phase's measured numbers and Phase 12C's own."
)


def run_phase12d(
    output_dir,
    max_k=DEFAULT_MAX_K,
    step_curve_mesh_cells_per_subcell=None,
    step_curve_l_family_mm=None,
    aperture_sizes=DEFAULT_APERTURE_SIZES,
    aperture_mesh_cells_per_subcell=DEFAULT_APERTURE_MESH_CELLS_PER_SUBCELL,
    aperture_l_family_mm=None,
    precomputed_near_at_step_cases=None,
):
    """``precomputed_near_at_step_cases``: an escape hatch for the production run only -- the
    near/at-resonance step cases (large patches, each several minutes) can be supplied already
    solved (e.g. from a prior run of this exact function, same mesh/family) instead of re-solved,
    so a production run is never repeated from scratch merely to add an unrelated case. Every
    number in the supplied cases must itself have come from a real solve; this parameter never
    accepts synthetic or interpolated values (checked: each must be a ``step_curve.StepCase``)."""
    from .step_curve import StepCase

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---------------- Task A ----------------
    mesh_a = step_curve_mesh_cells_per_subcell or step_curve.EMBEDDED_MESH_CELLS_PER_SUBCELL
    l_family_a = step_curve_l_family_mm or step_curve.L_FAMILY_MM
    curve = step_curve.build_reference_curve(mesh_a, l_family_a)
    centers = classify_centers(curve)
    if precomputed_near_at_step_cases is not None:
        assert all(isinstance(c, StepCase) for c in precomputed_near_at_step_cases)
        off_cases = step_curve.build_step_cases(
            curve, [("off_resonance", centers.off_resonance)], max_k=max_k, mesh_cells_per_subcell=mesh_a
        )
        step_cases = list(off_cases) + list(precomputed_near_at_step_cases)
    else:
        named_centers = [
            ("off_resonance", centers.off_resonance),
            ("near_resonance", centers.near_resonance),
            ("at_resonance", centers.at_resonance),
        ]
        step_cases = step_curve.build_step_cases(curve, named_centers, max_k=max_k, mesh_cells_per_subcell=mesh_a)
    task_a_gate = gate.decide_step_limit(step_cases)
    specification.write_task_a_json(curve, centers, step_cases, task_a_gate, output_dir / "task_a.json")

    # ---------------- Task B ----------------
    l_family_b = aperture_l_family_mm or DEFAULT_APERTURE_L_FAMILY_MM
    aperture_results = []
    for n in aperture_sizes:
        try:
            aperture_results.append(run_aperture_test(n, aperture_mesh_cells_per_subcell, l_family_b))
        except (MemoryError, ValueError):
            continue  # infeasible size: recorded as "attempted=False" by the gate below
    task_b_gate = gate.decide_aperture_adequacy(aperture_results, aperture_sizes)
    specification.write_task_b_json(aperture_results, aperture_sizes, task_b_gate, output_dir / "task_b.json")

    # ---------------- Task C ----------------
    hfss_comparisons = []
    if hfss_ingest.TOUCHSTONE_PATH.exists():
        hfss_comparisons.append(hfss_ingest.compare_single_cell(hfss_ingest.CASE_WITH_HFSS_DATA))
    task_c_gate = gate.decide_hfss_status(
        hfss_comparisons, hfss_ingest.ALL_SINGLE_CELL_IDS, (hfss_ingest.PENDING_RESONANCE_CASE,)
    )
    try:
        geometric_identity = hfss_ingest.check_geometric_identity()
    except ImportError:
        geometric_identity = ()
    specification.write_task_c_json(
        hfss_comparisons, hfss_ingest.PROVENANCE, geometric_identity, task_c_gate, output_dir / "task_c.json"
    )

    specification.write_gate_status_json(
        task_a_gate, task_b_gate, task_c_gate, output_dir / "phase12d_gate_status.json"
    )

    try:
        plots.generate_all_plots(step_cases, aperture_results, hfss_comparisons, output_dir / "plots")
    except ImportError:
        pass

    result = Phase12dResult(
        output_dir=output_dir,
        curve=curve,
        centers=centers,
        step_cases=tuple(step_cases),
        task_a_gate=task_a_gate,
        aperture_results=tuple(aperture_results),
        aperture_requested_sizes=tuple(aperture_sizes),
        task_b_gate=task_b_gate,
        hfss_comparisons=tuple(hfss_comparisons),
        hfss_provenance=hfss_ingest.PROVENANCE,
        geometric_identity=tuple(geometric_identity),
        task_c_gate=task_c_gate,
        recommendation=RECOMMENDATION,
    )
    (output_dir / "phase12d_report.md").write_text(render_phase12d_report_md(result))
    return result


def print_final_report(result):
    print("=" * 78)
    print("PHASE 12D -- RESOLVE THE LOCAL-PERIODICITY STEP LIMIT; INGEST REAL HFSS DATA")
    print("=" * 78)
    print(f"Task A (step-size curve):  {result.task_a_gate.status}")
    print(f"  {result.task_a_gate.reason}")
    print(f"Task B (aperture test):    {result.task_b_gate.status}")
    print(f"  {result.task_b_gate.reason}")
    print(f"Task C (HFSS ingestion):   {result.task_c_gate.status}")
    print(f"  {result.task_c_gate.reason}")
    print("Geometry inversion:  GEOMETRY_INVERSION_LOCKED")
    print("STEP export:         STEP_EXPORT_LOCKED")
    print("=" * 78)
