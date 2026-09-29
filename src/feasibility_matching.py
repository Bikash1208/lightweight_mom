"""Phase 8A Sec. 19-26: matching a required phase distribution against actual, already-simulated
library geometries. Every match is either an existing (L, W) cell within the configured tolerance,
or UNMATCHED -- never invented, interpolated, extrapolated, or rescued with an arbitrary phase
offset (Sec. 20, 21, 25, 26).
"""

from dataclasses import dataclass

from .phase import wrapped_difference_deg


@dataclass(frozen=True)
class MatchResult:
    required_phase_deg: float
    matched: bool
    nearest_library_phase_deg: float | None
    phase_error_deg: float | None  # signed, wrapped_difference_deg(nearest, required)
    candidate_L_mm: float | None
    candidate_W_mm: float | None
    candidate_validity: str | None  # "VALID" | None (no candidate at all)


def match_required_phase(required_phase_deg_array, table, tolerance_deg):
    """Sec. 21: for each required phase sample, the single *nearest* actual library candidate
    (correct circular metric, Sec. 26), kept as a match only if within ``tolerance_deg`` -- an
    explicit, caller-supplied parameter (Sec. 20 forbids inventing one silently)."""
    if tolerance_deg <= 0.0:
        raise ValueError(f"tolerance_deg must be > 0, got {tolerance_deg!r}")
    if not table:
        raise ValueError("match_required_phase: library phase table is empty")
    results = []
    for required in required_phase_deg_array:
        best_err, best_row = None, None
        for row in table:
            err = float(wrapped_difference_deg(row["phase_deg"], required))
            if best_err is None or abs(err) < abs(best_err):
                best_err, best_row = err, row
        matched = abs(best_err) <= tolerance_deg
        results.append(
            MatchResult(
                required_phase_deg=float(required),
                matched=matched,
                nearest_library_phase_deg=best_row["phase_deg"],
                phase_error_deg=best_err,
                candidate_L_mm=best_row["L_mm"],
                candidate_W_mm=best_row["W_mm"],
                candidate_validity="VALID",
            )
        )
    return results


@dataclass(frozen=True)
class CoverageSufficiencyResult:
    sufficient: bool
    n_required: int
    n_matched: int
    n_unmatched: int
    mean_error_deg: float | None
    max_error_deg: float | None
    unmatched_samples: tuple


def is_phase_coverage_sufficient(required_phase_deg_array, table, tolerance_deg):
    """Sec. 22: sufficient iff *every* required sample has a matched library candidate. This
    compares against the actual per-sample required phase distribution, never a naive
    min_required >= min_library / max_required <= max_library range check (Sec. 22's explicit
    prohibition), because phase is periodic and the library's usable phase need not be contiguous
    or centred anywhere in particular."""
    matches = match_required_phase(required_phase_deg_array, table, tolerance_deg)
    matched = [m for m in matches if m.matched]
    unmatched = [m for m in matches if not m.matched]
    errors = [abs(m.phase_error_deg) for m in matched]
    return CoverageSufficiencyResult(
        sufficient=len(matches) > 0 and len(unmatched) == 0,
        n_required=len(matches),
        n_matched=len(matched),
        n_unmatched=len(unmatched),
        mean_error_deg=(sum(errors) / len(errors)) if errors else None,
        max_error_deg=max(errors) if errors else None,
        unmatched_samples=tuple(m.required_phase_deg for m in unmatched),
    )
