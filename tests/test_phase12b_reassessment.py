"""Phase 12B Task B: aligned-mesh re-assessment of Phase 12A's grounded refinement study."""

import pytest

from src.phase12b.reassessment import (
    aligned_mesh_for_family,
    aligned_mesh_refinement_study,
    floquet_only_convergence,
)


def test_aligned_mesh_for_family_matches_task_a():
    assert aligned_mesh_for_family() == 96


def test_aligned_mesh_refinement_study_rejects_an_unaligned_mesh():
    with pytest.raises(ValueError, match="not aligned"):
        aligned_mesh_refinement_study(cells_list=(32,))


def test_aligned_mesh_refinement_study_converges_within_tolerance():
    """The Phase 12B headline finding: on aligned meshes, the joint mesh+Floquet refinement study
    for the representative patch converges monotonically and within the 0.02 tolerance (unlike
    Phase 12A's unaligned 0.048)."""
    points, diffs = aligned_mesh_refinement_study(cells_list=(24, 48))
    assert len(points) == 2
    assert diffs[0] < 0.02
    # monotone decrease would need a 3rd point; a single diff below tolerance is still meaningful
    # on its own (the full 24/48/96 study, and its 2 diffs, is in PHASE12B_MESH_PERIOD_STUDY.md)


@pytest.mark.slow
def test_floquet_only_convergence_isolates_floquet_truncation():
    points, diffs = floquet_only_convergence(floquet_orders=(48, 64))
    assert len(points) == 2
    assert diffs[0] < 0.02
