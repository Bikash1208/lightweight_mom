"""Sec. 19 (Phase 11O): computational reduction is only *considered* after the coupling physics is
established (Sec. 19's own ordering), and only when the classification actually supports it.

This project already applies one genuine, mathematically EXACT reduction throughout Phase 10-11:
``reference_solver.solve_grid``/``phase10.grid_solver.solve_grid`` solve the full heterogeneous
system ONCE and read out as many sub-cells' responses as needed from that single solve (`dipole_at`/
`response_at`) -- "matrix reuse" and "cached... reuse" (Sec. 19's own list) in their most literal
form, with zero approximation (it is the same self-consistent current solution, not a re-solve).
Nothing further (block-Toeplitz exploitation, hierarchical matrices, domain decomposition, etc.) is
implemented here: Sec. 19 requires any reduction to be validated against the full direct MoM
calculation before being trusted, and building one of those without first knowing whether the
underlying physics is even local enough to make such a reduction meaningful would be exactly the
"trying to make the problem look easier than it is" Sec. 36 warns against.
"""

from dataclasses import dataclass

from .classification import CASE_A_LOCAL

REDUCTION_VALIDATED = "REDUCTION_VALIDATED"
REDUCTION_NOT_ATTEMPTED = "REDUCTION_NOT_ATTEMPTED"


@dataclass(frozen=True)
class ReductionDecision:
    status: str
    rationale: str


def decide_reduction(overall_case):
    if overall_case == CASE_A_LOCAL:
        return ReductionDecision(
            REDUCTION_VALIDATED,
            "CASE_A_LOCAL: the already-implemented one-solve/many-read-out reuse in "
            "`reference_solver.solve_grid` is itself an exact, validated reduction relative to "
            "solving a separate system per sub-cell -- no further, approximate reduction is needed "
            "or attempted.",
        )
    return ReductionDecision(
        REDUCTION_NOT_ATTEMPTED,
        f"overall classification is {overall_case}, not CASE_A_LOCAL: building an approximate "
        "spatial-reduction scheme (hierarchical matrices, block-Toeplitz reuse, domain "
        "decomposition, etc.) before knowing the coupling is local enough for such a scheme to be "
        "meaningful would not be a physically justified reduction (Sec. 19's own requirement); none "
        "is implemented in this phase.",
    )
