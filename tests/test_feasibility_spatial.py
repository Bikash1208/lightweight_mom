"""Phase 8A Sec. 29-31: sampling compatibility and spatial-mapping availability (brief Sec. 44
Test 12: the actual 2 mm reference sampling / 4.8 mm MTS period case)."""

from src.feasibility_spatial import check_mapping_method_defined, check_sampling_compatibility, estimate_max_aperture


def test_mismatched_2mm_vs_4p8mm_grid_is_flagged_non_identical():
    result = check_sampling_compatibility(2.0, 2.0, 4.8, 4.8)
    assert result.grids_identical is False
    assert "non-identical" in result.note
    assert "validated spatial resampling/mapping method is required" in result.note


def test_identical_grid_is_recognized():
    result = check_sampling_compatibility(4.8, 4.8, 4.8, 4.8)
    assert result.grids_identical is True


def test_missing_reference_spacing_is_not_identical():
    result = check_sampling_compatibility(None, None, 4.8, 4.8)
    assert result.grids_identical is False


def test_mapping_method_is_not_defined_in_this_project():
    result = check_mapping_method_defined()
    assert result.defined is False
    assert "No validated" in result.note
    assert "Nearest-neighbour" in result.note


def test_estimate_max_aperture_is_a_floor_division_upper_bound():
    est = estimate_max_aperture(50.0, 50.0, 4.8, 4.8)
    assert est.max_nx == 10  # floor(50/4.8) = 10
    assert est.max_ny == 10


def test_estimate_max_aperture_reports_unknown_without_extent():
    est = estimate_max_aperture(None, None, 4.8, 4.8)
    assert est.max_nx is None and est.max_ny is None
    assert "unknown" in est.note
