"""Phase 15 Part C: multi-feed reference-field construction (superposition and segmentation), and
its plumbing into the existing single-reference-field hologram synthesis functions unmodified."""

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.hologram import hmsa_synthesis as hs
from src.hologram.multi_feed import (
    FeedSpec,
    feed_owner_index,
    feed_xy_m_for_beams,
    segmented_reference_field_and_regions,
    superposed_reference_field,
)
from src.phase8b.object_wave import Beam
from src.phase8b.reference_wave import monopole_reference_field
from src.phase8e.reference_field_input import phase_only_reference_field

FREQ = 10e9
BETA = tm0_beta(FREQ, 2.2, 1.575e-3)


def _grid(n=9, extent_m=0.02):
    coords = np.linspace(-extent_m, extent_m, n)
    return np.meshgrid(coords, coords, indexing="xy")


# =========================================================================== FeedSpec


def test_feed_spec_rejects_non_finite_position():
    with pytest.raises(ValueError):
        FeedSpec(x_m=float("nan"), y_m=0.0)


# ============================================================== superposed_reference_field


def test_superposition_of_one_feed_at_origin_matches_monopole_reference_field():
    x, y = _grid()
    feed = FeedSpec(x_m=0.0, y_m=0.0)
    combined = superposed_reference_field([feed], x, y, k_swg=BETA)
    expected = monopole_reference_field(x, y, BETA)
    assert np.allclose(combined.Ex, expected)


def test_superposition_of_two_feeds_is_a_complex_sum():
    x, y = _grid()
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0), FeedSpec(x_m=0.01, y_m=0.0)]
    combined = superposed_reference_field(feeds, x, y, k_swg=BETA)
    expected = monopole_reference_field(x - feeds[0].x_m, y - feeds[0].y_m, BETA) + monopole_reference_field(
        x - feeds[1].x_m, y - feeds[1].y_m, BETA
    )
    assert np.allclose(combined.Ex, expected)
    assert combined.amplitude_available


def test_superposition_applies_each_feeds_own_amplitude():
    x, y = _grid()
    feeds_unit = [FeedSpec(x_m=0.0, y_m=0.0, amplitude=1.0)]
    feeds_scaled = [FeedSpec(x_m=0.0, y_m=0.0, amplitude=2.0 + 1.0j)]
    unit = superposed_reference_field(feeds_unit, x, y, k_swg=BETA)
    scaled = superposed_reference_field(feeds_scaled, x, y, k_swg=BETA)
    assert np.allclose(scaled.Ex, (2.0 + 1.0j) * unit.Ex)


def test_superposition_requires_k_swg_for_analytic_feeds():
    x, y = _grid()
    with pytest.raises(ValueError, match="k_swg"):
        superposed_reference_field([FeedSpec(x_m=0.0, y_m=0.0)], x, y, k_swg=None)


def test_superposition_can_mix_an_uploaded_field_with_an_analytic_one():
    x, y = _grid()
    uploaded = phase_only_reference_field(x, y, np.degrees(np.arctan2(y, x)))
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0, reference_field=uploaded), FeedSpec(x_m=0.01, y_m=0.0)]
    combined = superposed_reference_field(feeds, x, y, k_swg=BETA)
    assert combined.Ex.shape == x.shape
    assert np.all(np.isfinite(combined.Ex))


def test_superposition_of_empty_feeds_raises():
    x, y = _grid()
    with pytest.raises(ValueError):
        superposed_reference_field([], x, y, k_swg=BETA)


# =============================================================================== segmentation


def test_feed_owner_index_partitions_by_nearest_feed():
    x, y = np.array([-1.0, 1.0]), np.array([0.0, 0.0])
    feeds = [FeedSpec(x_m=-2.0, y_m=0.0), FeedSpec(x_m=2.0, y_m=0.0)]
    owner = feed_owner_index(feeds, x, y)
    assert list(owner) == [0, 1]


def test_segmentation_is_piecewise_not_a_sum():
    x, y = _grid()
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0), FeedSpec(x_m=0.01, y_m=0.0)]
    grid, owner = segmented_reference_field_and_regions(feeds, x, y, k_swg=BETA)
    field0 = monopole_reference_field(x - feeds[0].x_m, y - feeds[0].y_m, BETA)
    field1 = monopole_reference_field(x - feeds[1].x_m, y - feeds[1].y_m, BETA)
    expected = np.where(owner == 0, field0, field1)
    assert np.allclose(grid.Ex, expected)
    # NOT a sum: verify at least one cell differs from field0+field1
    assert not np.allclose(grid.Ex, field0 + field1)


def test_segmentation_of_empty_feeds_raises():
    x, y = _grid()
    with pytest.raises(ValueError):
        segmented_reference_field_and_regions([], x, y, k_swg=BETA)


# =============================================================================== feed_xy_m_for_beams


def test_feed_xy_m_for_beams_builds_the_sequence_format():
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0), FeedSpec(x_m=0.01, y_m=0.005)]
    xs, ys = feed_xy_m_for_beams(feeds, [0, 1, 0])
    assert xs == [-0.01, 0.01, -0.01]
    assert ys == [0.0, 0.005, 0.0]


# ===================================================== end-to-end through hmsa_synthesis (unmodified)


def _candidates():
    return [(f"c{i}", 1.0, 1.0, "lib", np.exp(1j * np.radians(i))) for i in range(-180, 180)]


def test_synthesize_formulation_b_accepts_a_superposed_multi_feed_reference():
    x, y = _grid()
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0), FeedSpec(x_m=0.01, y_m=0.0)]
    ref = superposed_reference_field(feeds, x, y, k_swg=BETA)
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=15.0, phi_deg=180.0)]
    feed_xy = feed_xy_m_for_beams(feeds, [0, 1])
    design = hs.synthesize_formulation_b(ref, beams, FREQ, _candidates(), period_m=8e-3, feed_xy_m=feed_xy)
    assert design.formulation == "B_phase_difference"
    assert len(design.cells) == design.x_m.size > 0
    assert all(c.reachable for c in design.cells)


def test_superposition_vs_segmentation_produce_different_designs():
    x, y = _grid()
    feeds = [FeedSpec(x_m=-0.01, y_m=0.0), FeedSpec(x_m=0.01, y_m=0.0)]
    beams = [Beam(theta_deg=20.0, phi_deg=0.0), Beam(theta_deg=15.0, phi_deg=180.0)]
    feed_xy = feed_xy_m_for_beams(feeds, [0, 1])

    sup_ref = superposed_reference_field(feeds, x, y, k_swg=BETA)
    seg_ref, _owner = segmented_reference_field_and_regions(feeds, x, y, k_swg=BETA)
    sup_design = hs.synthesize_formulation_b(sup_ref, beams, FREQ, _candidates(), period_m=8e-3, feed_xy_m=feed_xy)
    seg_design = hs.synthesize_formulation_b(seg_ref, beams, FREQ, _candidates(), period_m=8e-3, feed_xy_m=feed_xy)
    assert not np.allclose(sup_design.control_map, seg_design.control_map)


def test_single_feed_multi_feed_pipeline_matches_the_single_feed_path():
    """A degenerate "multi-feed" design of exactly one feed at the origin must reproduce the
    ordinary single-reference-field path exactly -- multi_feed.py adds no new physics for N=1.
    Both paths must start from the SAME complex field (amplitude included): comparing against a
    phase-only reference would also change how `resample_to_cell_grid` interpolates (unit-magnitude
    phasor vs amplitude-weighted complex interpolation), a real but unrelated difference this test
    is not about."""
    from dataclasses import replace

    x, y = _grid()
    beam = [Beam(theta_deg=20.0, phi_deg=0.0)]
    plain_field = monopole_reference_field(x, y, BETA)
    plain_ref = replace(
        phase_only_reference_field(x, y, np.zeros_like(x)),
        Ex=plain_field,
        Ey=np.zeros_like(plain_field),
        amplitude_available=True,
        phase_deg=np.degrees(np.angle(plain_field)),
    )
    single_design = hs.synthesize_formulation_b(plain_ref, beam, FREQ, _candidates(), period_m=8e-3)

    feeds = [FeedSpec(x_m=0.0, y_m=0.0)]
    multi_ref = superposed_reference_field(feeds, x, y, k_swg=BETA)
    multi_design = hs.synthesize_formulation_b(multi_ref, beam, FREQ, _candidates(), period_m=8e-3)
    assert np.array_equal(single_design.control_map, multi_design.control_map)
