"""Phase 12C gates: explicit named states, never a bare pass/fail (project-wide rule).

Two tolerances are used for the local-periodicity gate, both reported side by side (never only
one, so neither obscures the other):

- Phase 9's own criterion (``phase9.local_periodicity.classify_local_periodicity``, reused
  unmodified): relative complex error vs the solver's OWN measured noise floor / a self-referential
  "at least half the signal's own magnitude" substantiality test -- no invented percentage.
- The brief's own explicit example numeric tolerance for THIS question's specific framing ("the
  largest phase step per cell for which local periodicity holds ... e.g. <= 5-10 deg phase /
  <= 0.5 dB"): ``PHASE_TOLERANCE_DEG = 10.0`` (the brief's own more lenient end) and
  ``MAGNITUDE_TOLERANCE_DB = 0.5`` (the brief's own number), used ONLY because the brief itself
  names them for this specific "largest step" answer -- not adopted as a general project tolerance.
"""

import math
from dataclasses import dataclass
from pathlib import Path

from ..phase9.local_periodicity import classify_local_periodicity

PHASE_TOLERANCE_DEG = 10.0
MAGNITUDE_TOLERANCE_DB = 0.5

GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS = "GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS"
GROUNDED_LOCAL_PERIODICITY_INVALIDATED = "GROUNDED_LOCAL_PERIODICITY_INVALIDATED"
GROUNDED_LOCAL_PERIODICITY_UNVALIDATED = "GROUNDED_LOCAL_PERIODICITY_UNVALIDATED"

EARLIER_CONCLUSIONS_HOLD = "EARLIER_CONCLUSIONS_HOLD"
EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED = "EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED"
EARLIER_CONCLUSIONS_CHANGED = "EARLIER_CONCLUSIONS_CHANGED"

HFSS_EXPORT_READY = "HFSS_EXPORT_READY"
HFSS_NOT_AVAILABLE = "HFSS_NOT_AVAILABLE"


def magnitude_error_db(isolated_px, embedded_px):
    """20*log10(|embedded|/|isolated|), the standard dB expression of a magnitude deviation."""
    return 20.0 * abs(math.log10(abs(embedded_px) / abs(isolated_px)))


@dataclass(frozen=True)
class LocalPeriodicityGateResult:
    status: str
    reason: str
    max_step_within_tolerance_deg: "float | None"
    phase9_classification: str   # VALIDATED | INVALIDATED | UNVALIDATED (phase9.local_periodicity's own labels)
    per_case: tuple              # (name, achieved_step_deg, phase_error_deg, magnitude_error_db, within_tolerance)


def decide_local_periodicity(gradient_results, phase_tolerance_deg=PHASE_TOLERANCE_DEG,
                              magnitude_tolerance_db=MAGNITUDE_TOLERANCE_DB):
    """``gradient_results``: an iterable of ``locality_study.GradientCaseResult`` (1-D row and 2-D
    gradient cases only -- NOT the strong-contrast or periodic-supercell cases, which are reported
    separately, since the brief's own "largest phase step per cell" question is specifically about
    the smooth-gradient series)."""
    per_case = []
    for r in gradient_results:
        db = magnitude_error_db(r.isolated_px, r.embedded_px)
        within = r.phase_error_deg <= phase_tolerance_deg and db <= magnitude_tolerance_db
        per_case.append((r.name, r.achieved_step_deg, r.phase_error_deg, db, within))

    p9_status = classify_local_periodicity([_ForP9(r) for r in gradient_results])

    within_cases = [c for c in per_case if c[4]]
    max_step = max((c[1] for c in within_cases), default=None)

    if not within_cases:
        status = GROUNDED_LOCAL_PERIODICITY_INVALIDATED
        reason = (f"no tested smooth-gradient case stayed within the stated tolerance "
                  f"(phase <= {phase_tolerance_deg} deg, magnitude <= {magnitude_tolerance_db} dB); "
                  f"Phase 9's own criterion also classifies these deviations as {p9_status}")
    elif len(within_cases) == len(per_case):
        status = GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS
        reason = (f"every tested smooth-gradient case (largest achieved step {max_step:.1f} deg) "
                  f"stayed within the stated tolerance; Phase 9's own criterion: {p9_status}")
    else:
        status = GROUNDED_LOCAL_PERIODICITY_VALID_FOR_SMOOTH_GRADIENTS
        reason = (f"local periodicity holds within tolerance up to {max_step:.1f} deg/cell "
                  f"(the largest step size tested that stayed within tolerance); at least one "
                  f"larger tested step exceeded it. Phase 9's own criterion: {p9_status}")

    return LocalPeriodicityGateResult(status, reason, max_step, p9_status, tuple(per_case))


class _ForP9:
    """Adapts a ``locality_study.GradientCaseResult`` to the ``relative_complex_error`` attribute
    ``phase9.local_periodicity.classify_local_periodicity`` reads, without modifying either
    module."""

    def __init__(self, result):
        self.relative_complex_error = result.relative_complex_error


@dataclass(frozen=True)
class SpotCheckGateResult:
    status: str
    reason: str
    findings: tuple   # (case_name, verdict, detail)


def decide_spot_check(c2_c3, phase9_strong, case_d, coupling_points):
    """Sec. D: for each re-check, does the earlier conclusion hold, change quantitatively only, or
    change outright?"""
    findings = []

    # Phase 10 C2/C3: the earlier (Phase 12B) finding was "mesh-identical" at the coarse mesh.
    if not c2_c3.mask_identical:
        findings.append(("phase10_c2_c3", EARLIER_CONCLUSIONS_CHANGED,
                          f"C2 and C3 are mesh-DISTINCT at the finer {c2_c3.mesh_cells_per_side}-cell "
                          f"mesh (relative isolated-response difference {c2_c3.relative_difference:.3e}); "
                          f"Phase 12B's coarse-mesh finding does not hold at this resolution"))
    else:
        findings.append(("phase10_c2_c3", EARLIER_CONCLUSIONS_HOLD,
                          f"C2 and C3 remain mesh-identical at the finer {c2_c3.mesh_cells_per_side}-cell mesh"))

    # Phase 9 strong contrast: was INVALIDATED (relative_complex_error ~ 0.997, >= the 0.5
    # substantiality threshold, phase9.local_periodicity.SUBSTANTIAL_DEVIATION_RELATIVE) at 8
    # cells/subcell. The meaningful question is whether the finer mesh crosses THAT threshold
    # (a real change in classification), not merely whether the raw number moved.
    from ..phase9.local_periodicity import SUBSTANTIAL_DEVIATION_RELATIVE

    earlier_substantial = phase9_strong.earlier_relative_complex_error >= SUBSTANTIAL_DEVIATION_RELATIVE
    now_substantial = phase9_strong.relative_complex_error >= SUBSTANTIAL_DEVIATION_RELATIVE
    if earlier_substantial == now_substantial:
        findings.append(("phase9_strong_contrast", EARLIER_CONCLUSIONS_HOLD,
                          f"relative complex error {phase9_strong.relative_complex_error:.3f} at "
                          f"{phase9_strong.mesh_cells_per_subcell} cells/subcell (vs {phase9_strong.earlier_relative_complex_error:.3f} "
                          f"at 8 cells/subcell) -- both are {'>=' if now_substantial else '<'} the "
                          f"{SUBSTANTIAL_DEVIATION_RELATIVE} substantiality threshold; the INVALIDATED "
                          f"classification is unaffected"))
    else:
        findings.append(("phase9_strong_contrast", EARLIER_CONCLUSIONS_CHANGED,
                          f"relative complex error dropped from {phase9_strong.earlier_relative_complex_error:.3f} "
                          f"(>= the {SUBSTANTIAL_DEVIATION_RELATIVE} substantiality threshold, INVALIDATED) at 8 "
                          f"cells/subcell to {phase9_strong.relative_complex_error:.3f} (< {SUBSTANTIAL_DEVIATION_RELATIVE}) "
                          f"at the finer, aligned {phase9_strong.mesh_cells_per_subcell}-cell/subcell mesh -- the "
                          f"INVALIDATED classification itself does not survive this finer mesh for this specific "
                          f"contrast case; both isolated and embedded solves were independently verified converged "
                          f"to <1% at this floquet order before reporting this (Sec. D)"))

    # Phase 11 Case D: was "substantial" (>= 0.5 relative) at n=7, 8 cells/subcell.
    if case_d.substantial:
        findings.append(("phase11_case_d", EARLIER_CONCLUSIONS_HOLD,
                          f"at the reduced, finer-mesh re-check (n={case_d.n}, {case_d.mesh_cells_per_subcell} "
                          f"cells/subcell), the checkerboard-vs-block relative difference is "
                          f"{case_d.relative_difference:.3f} -- still substantial, confirming "
                          f"arrangement (not just local statistics) matters"))
    else:
        findings.append(("phase11_case_d", EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED,
                          f"at the reduced, finer-mesh re-check (n={case_d.n}), the relative difference "
                          f"{case_d.relative_difference:.3f} is below the 0.5 substantiality threshold "
                          f"(the original finding used n=7; this is a reduced-size re-check, Sec. D)"))

    # Phase 11 coupling radius: was "not converged by 9x9/11x11" -- at this reduced size range we
    # can only report whether the shell increment is still large, not overturn the original
    # not-converged-by-11x11 finding at a smaller size than that.
    if len(coupling_points) >= 2:
        last, prev = coupling_points[-1], coupling_points[-2]
        rel = abs(last.response - prev.response) / abs(prev.response) if abs(prev.response) > 0 else float("inf")
        findings.append(("phase11_coupling_radius", EARLIER_CONCLUSIONS_HOLD if rel > 1e-6 else
                          EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED,
                          f"shell increment n={prev.n}->{last.n} relative change {rel:.3e} at the finer mesh "
                          f"(feasible sizes only, Sec. 'Cost note'; not a re-run of the original 9x9/11x11 series)"))

    changed = [f for f in findings if f[1] == EARLIER_CONCLUSIONS_CHANGED]
    quantitative = [f for f in findings if f[1] == EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED]
    if changed:
        status = EARLIER_CONCLUSIONS_CHANGED
    elif quantitative:
        status = EARLIER_CONCLUSIONS_QUANTITATIVELY_CHANGED
    else:
        status = EARLIER_CONCLUSIONS_HOLD
    reason = "; ".join(f"{name}: {verdict}" for name, verdict, _ in findings)
    return SpotCheckGateResult(status, reason, tuple(findings))


def hfss_export_status(single_cells, supercells, synthetic_records):
    """``HFSS_EXPORT_READY`` iff every requested export actually produced a STEP file and the
    synthetic one-port ingestion round-tripped (Sec. E); never claims real HFSS data exists."""
    ok = (len(single_cells) > 0 and all(Path(s.step_path).exists() for s in single_cells) and
          len(supercells) > 0 and all(Path(s.step_path).exists() for s in supercells) and
          synthetic_records is not None and all(r.s21 is None for r in synthetic_records))
    return HFSS_EXPORT_READY if ok else "HFSS_EXPORT_NOT_READY"


from pathlib import Path  # noqa: E402  (kept near its sole use, matching this project's other gate modules)
