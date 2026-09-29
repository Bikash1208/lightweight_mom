"""Phase 12D gates: explicit named states, never a bare pass/fail (project-wide rule)."""

from dataclasses import dataclass

# Two tolerances, both always reported side by side (brief's own explicit request: "use a
# design-relevant one, e.g. 10deg/1dB, and also show 5deg/0.5dB").
DESIGN_TOLERANCE_DEG = 10.0
DESIGN_TOLERANCE_DB = 1.0
STRICT_TOLERANCE_DEG = 5.0
STRICT_TOLERANCE_DB = 0.5

LOCAL_PERIODICITY_STEP_LIMIT_MEASURED = "LOCAL_PERIODICITY_STEP_LIMIT_MEASURED"
LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED = "LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED"

LOOKUP_TABLE_ADEQUATE = "LOOKUP_TABLE_ADEQUATE"
LOOKUP_TABLE_INADEQUATE = "LOOKUP_TABLE_INADEQUATE"
LOOKUP_TABLE_UNTESTED = "LOOKUP_TABLE_UNTESTED"

HFSS_SINGLE_CELL_VALIDATED_PARTIAL = "HFSS_SINGLE_CELL_VALIDATED_PARTIAL"
HFSS_SINGLE_CELL_NOT_VALIDATED = "HFSS_SINGLE_CELL_NOT_VALIDATED"
HFSS_SUPERCELL_NOT_AVAILABLE = "HFSS_SUPERCELL_NOT_AVAILABLE"


@dataclass(frozen=True)
class StepLimitGateResult:
    status: str
    reason: str
    max_step_within_design_tolerance_deg: "float | None"
    max_step_within_strict_tolerance_deg: "float | None"
    per_case: tuple  # (name, achieved_step_deg, phase_error_deg, mag_error_db, within_design, within_strict)


def decide_step_limit(step_curve_cases):
    """``step_curve_cases``: iterable of objects with ``.name``, ``.achieved_step_deg``,
    ``.reflection_phase_error_deg``, ``.reflection_magnitude_error_db`` (Task A's own reflection-
    coefficient metric, not the raw dipole one). Reports the largest tested step that stays within
    each tolerance -- never invents a step size that was not actually tested (brief: "report the
    achieved steps honestly")."""
    per_case = []
    for c in step_curve_cases:
        within_design = (
            c.reflection_phase_error_deg <= DESIGN_TOLERANCE_DEG
            and c.reflection_magnitude_error_db <= DESIGN_TOLERANCE_DB
        )
        within_strict = (
            c.reflection_phase_error_deg <= STRICT_TOLERANCE_DEG
            and c.reflection_magnitude_error_db <= STRICT_TOLERANCE_DB
        )
        per_case.append(
            (
                c.name,
                c.achieved_step_deg,
                c.reflection_phase_error_deg,
                c.reflection_magnitude_error_db,
                within_design,
                within_strict,
            )
        )

    design_ok = [c for c in per_case if c[4]]
    strict_ok = [c for c in per_case if c[5]]
    max_design = max((c[1] for c in design_ok), default=None)
    max_strict = max((c[1] for c in strict_ok), default=None)

    if not per_case:
        return StepLimitGateResult(
            LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED,
            "no step-curve case was measured",
            None,
            None,
            (),
        )
    if not design_ok:
        status = LOCAL_PERIODICITY_STEP_LIMIT_NOT_RESOLVED
        reason = (
            f"no tested case (achieved steps {min(c[1] for c in per_case):.1f}-"
            f"{max(c[1] for c in per_case):.1f} deg) stayed within the {DESIGN_TOLERANCE_DEG} "
            f"deg / {DESIGN_TOLERANCE_DB} dB design tolerance; the true crossover step (if any "
            f"smaller step exists) was not computationally reachable at this mesh -- see the "
            f"report's cost note"
        )
    else:
        status = LOCAL_PERIODICITY_STEP_LIMIT_MEASURED
        reason = (
            f"local periodicity holds within the design tolerance up to {max_design:.1f} "
            f"deg/cell (largest tested step that stayed within it)"
            + (
                f"; within the stricter {STRICT_TOLERANCE_DEG} deg/{STRICT_TOLERANCE_DB} dB "
                f"tolerance up to {max_strict:.1f} deg/cell"
                if strict_ok
                else f"; no tested step stayed within the stricter {STRICT_TOLERANCE_DEG} deg/"
                f"{STRICT_TOLERANCE_DB} dB tolerance"
            )
        )
    return StepLimitGateResult(status, reason, max_design, max_strict, tuple(per_case))


@dataclass(frozen=True)
class ApertureGateResult:
    status: str
    reason: str
    per_size: tuple  # (n_cells, fraction_intended, max_spurious_fraction, adequate, attempted)


def decide_aperture_adequacy(aperture_results, requested_sizes, adequacy_threshold=0.9):
    """``adequate`` iff the intended order carries at least ``adequacy_threshold`` of the total
    propagating reflected power (a design-relevant, stated threshold, not an invented pass mark
    hidden in the code -- reported alongside the raw fraction so the reader can apply a different
    one). Sizes in ``requested_sizes`` that never produced a result (infeasible; brief: "if N=12 is
    infeasible, say so") are recorded with ``attempted=False``, never silently dropped."""
    attempted = {r.n_cells: r for r in aperture_results}
    per_size = []
    for n in requested_sizes:
        r = attempted.get(n)
        if r is None:
            per_size.append((n, None, None, None, False))
            continue
        max_spurious = max((f for _, f in r.spurious_orders), default=0.0)
        adequate = r.fraction_intended >= adequacy_threshold
        per_size.append((n, r.fraction_intended, max_spurious, adequate, True))

    tested = [p for p in per_size if p[4]]
    if not tested:
        return ApertureGateResult(LOOKUP_TABLE_UNTESTED, "no requested aperture size was feasible", tuple(per_size))
    adequate_ns = sorted(n for n, _, _, ok, _ in tested if ok)
    if adequate_ns:
        status = LOOKUP_TABLE_ADEQUATE
        reason = (
            f"the lookup-table-constructed aperture puts >= {adequacy_threshold * 100:.0f}% of "
            f"reflected power in the intended order for N >= {min(adequate_ns)} (tested: "
            f"{sorted(n for n, *_ in tested)})"
        )
    else:
        status = LOOKUP_TABLE_INADEQUATE
        reason = (
            f"no tested, feasible N reached {adequacy_threshold * 100:.0f}% of reflected power "
            f"in the intended order (measured: " + ", ".join(f"N={n}: {f:.3f}" for n, f, *_ in tested) + ")"
        )
    return ApertureGateResult(status, reason, tuple(per_size))


@dataclass(frozen=True)
class HfssGateResult:
    status: str
    reason: str
    cases_with_data: tuple
    cases_pending: tuple


def decide_hfss_status(comparisons, all_case_ids, pending_ids):
    """``comparisons``: iterable of objects with ``.case_id`` and ``.meets_tolerance``. Never
    claims validation for a case with no HFSS file on disk; a case with data but outside tolerance
    is still listed (in ``cases_with_data``), just not counted toward the PARTIAL status without
    also checking ``.meets_tolerance``."""
    by_id = {c.case_id: c for c in comparisons}
    validated = tuple(cid for cid, c in by_id.items() if c.meets_tolerance)
    pending = tuple(cid for cid in all_case_ids if cid not in by_id)
    if not by_id:
        return HfssGateResult(
            "HFSS_SINGLE_CELL_NOT_VALIDATED",
            "no HFSS Touchstone data found for any single-cell case",
            (),
            tuple(all_case_ids),
        )
    status = HFSS_SINGLE_CELL_VALIDATED_PARTIAL if validated else "HFSS_SINGLE_CELL_NOT_VALIDATED"
    reason = (
        f"validated within tolerance: {validated or '(none)'}; cases with HFSS data but "
        f"outside tolerance: {tuple(cid for cid in by_id if cid not in validated) or '(none)'}; "
        f"still pending (no HFSS file yet): {pending}"
    )
    return HfssGateResult(status, reason, validated, pending)
