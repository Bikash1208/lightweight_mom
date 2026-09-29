"""Phase 12D Task C: the Touchstone reader (phase9.hfss_reference.read_touchstone) and the
single-cell comparison it feeds must reproduce the real HFSS data on disk and the already-measured
G1 result the brief cites. Skips cleanly if the (large, real-data) export directory is absent."""

import pytest

from src.phase12d import hfss_ingest as hi

pytestmark = [
    pytest.mark.hfss_data,
    pytest.mark.skipif(not hi.TOUCHSTONE_PATH.exists(), reason="real HFSS export not present in this checkout"),
]


def test_compare_single_cell_g1_reproduces_recorded_comparison():
    c = hi.compare_single_cell("G1")
    assert c.l_mm == 2.0
    assert hi.reproduces_recorded_g1_comparison(c)


def test_g1_mean_and_max_phase_deviation_match_brief():
    c = hi.compare_single_cell("G1")
    assert abs(c.mean_dphase_deg - (-0.29)) < 0.05
    assert 0.1 < c.max_abs_dphase_deg < 1.0


def test_g1_magnitude_and_crosspol_are_tiny():
    c = hi.compare_single_cell("G1")
    assert c.max_dmag_db < 0.01
    assert c.max_crosspol_mag < 1e-4


def test_g1_top_to_bottom_transmission_is_zero_grounded():
    c = hi.compare_single_cell("G1")
    assert c.max_top_to_bottom_mag == 0.0


def test_g1_meets_the_declared_tolerance():
    c = hi.compare_single_cell("G1")
    assert c.meets_tolerance


def test_manifest_loads_and_lists_expected_cases():
    manifest = hi.load_manifest()
    ids = {s["id"] for s in manifest["singles"]}
    assert ids == set(hi.ALL_SINGLE_CELL_IDS)
