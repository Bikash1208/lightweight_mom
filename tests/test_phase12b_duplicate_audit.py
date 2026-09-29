"""Phase 12B Task C: read-only duplicate-geometry audit."""

import pytest

from src.phase8c.data_source import DEFAULT_CHECKPOINT_JSON
from src.phase12b.duplicate_audit import (
    phase7_delivered_lookup_duplicate_count,
    phase7_library_duplicate_report,
    phase9_strong_contrast_check,
    phase10_center_duplicate_report,
    phase10_center_perturbation_duplicate_report,
    phase10_holdout_leakage_check,
    phase11_pattern_endpoints_check,
)


def test_phase7_library_duplicate_report_matches_the_brief_own_documented_count():
    """The brief itself measured 169 distinct / 105 duplicate groups by equal *solved S*; this
    independently re-derives it from the discretized *metal mask* (frequency-independent) and
    they must agree exactly, since identical masks solve to identical S at every frequency."""
    r = phase7_library_duplicate_report()
    assert r.n_geometries == 324
    assert r.n_distinct == 169
    assert r.n_duplicate_groups == 105


@pytest.mark.skipif(not DEFAULT_CHECKPOINT_JSON.exists(), reason="Phase 7 checkpoint not present")
@pytest.mark.slow
def test_delivered_lookup_cross_check_matches_the_geometric_count():
    result = phase7_delivered_lookup_duplicate_count()
    assert result["n_cells"] == 324
    assert result["n_distinct_solved"] == 169


def test_phase9_strong_contrast_geometries_are_not_mesh_degenerate():
    """Phase 9's LOCAL_PERIODICITY_INVALIDATED finding depends on the center and the strong-
    contrast neighbors being genuinely distinct geometries; confirm they are (no duplicate
    groups) at Phase 9's own solve mesh."""
    assert phase9_strong_contrast_check() == ()


def test_phase10_c2_and_c3_are_mesh_identical_at_the_supercell_mesh():
    """A real, measured finding (Phase 12B): C2=(2.0,2.8)mm and C3=(2.8,2.8)mm discretize to the
    exact same metal mask at Phase 9-11's 8-cell/side mesh (dx=0.6mm)."""
    groups = phase10_center_duplicate_report()
    assert groups == (("C2", "C3"),)


def test_phase10_holdout_center_c5_has_no_leakage_into_training_centers():
    """The blind-test/train split (Phase 10's own geometry-based holdout by center name) is not
    compromised by mesh duplication: none of C1-C4 (training) collides with C5 (held out)."""
    result = phase10_holdout_leakage_check()
    assert result["leaks"] == ()


def test_phase10_cross_center_perturbation_collisions_are_confined_to_the_training_centers():
    result = phase10_center_perturbation_duplicate_report()
    involved = {name.split("_L")[0] for group in result["cross_center_groups"] for name in group}
    assert involved <= {"C2", "C3", "C4"}
    assert "C1" not in involved and "C5" not in involved


def test_phase11_checkerboard_endpoints_are_not_mesh_identical():
    assert phase11_pattern_endpoints_check() is False
