"""Sec. 16: the convergence criterion, declared BEFORE any Phase 11 result was evaluated -- these
are not new numbers chosen to fit the coupling-radius data: they are ``phase9.local_periodicity``'s
own already-sourced constants (a measured solver self-consistency floor, and a self-referential
"at least half the signal's own magnitude" substantiality test), the same reuse
``phase10.neighborhood_decision`` already made for exactly this reason. Reusing a tolerance that
predates this phase is the most defensible way to satisfy "declare the tolerance before evaluating
the final result."
"""

from dataclasses import dataclass

from ..phase9.local_periodicity import NOISE_FLOOR_RELATIVE, SUBSTANTIAL_DEVIATION_RELATIVE

CONVERGED = "CONVERGED"
NOT_CONVERGED = "NOT_CONVERGED"
STILL_SUBSTANTIALLY_CHANGING = "STILL_SUBSTANTIALLY_CHANGING"


@dataclass(frozen=True)
class ConvergenceResult:
    center_name: str
    status: str
    last_increment_relative: float
    max_increment_relative: float
    converged_at_n: object  # the smallest n whose OWN incoming shell increment was <= the floor, or None


def classify_convergence(
    increments, center_name="", noise_floor=NOISE_FLOOR_RELATIVE, substantial=SUBSTANTIAL_DEVIATION_RELATIVE
):
    """Sec. 16: `R_(N+2) ~= R_N` only when the shell increment is at or below the solver's own
    measured self-consistency floor -- a solver-precision fact, not a judgement call. If ANY
    increment (not just the last) is still at least "substantial" relative to its own base value,
    the series has not settled down by the largest size tested."""
    if not increments:
        raise ValueError("classify_convergence: no shell increments given")
    last = increments[-1]
    max_rel = max(inc.delta_relative for inc in increments)
    converged_at_n = next((inc.n_to for inc in increments if inc.delta_relative <= noise_floor), None)

    if last.delta_relative <= noise_floor:
        status = CONVERGED
    elif any(inc.delta_relative >= substantial for inc in increments):
        status = STILL_SUBSTANTIALLY_CHANGING
    else:
        status = NOT_CONVERGED

    return ConvergenceResult(center_name, status, last.delta_relative, max_rel, converged_at_n)
