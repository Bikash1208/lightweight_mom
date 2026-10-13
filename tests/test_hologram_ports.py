"""Phase 15 Part C (gap-fill): sub-aperture port layouts and the first-order isolation estimate."""

import math

import pytest

from src.hologram.ports import (
    IsolationFlag,
    PortRegion,
    custom_layout,
    isolation_estimate,
    layout_preset_1xn,
    layout_preset_2x2_plus_1_large,
)


def test_port_region_rejects_degenerate_rectangle():
    with pytest.raises(ValueError):
        PortRegion(0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0)


def test_port_region_area_and_contains():
    r = PortRegion(0, -1.0, 1.0, -2.0, 2.0, 0.0, 0.0)
    assert r.width_m == pytest.approx(2.0)
    assert r.height_m == pytest.approx(4.0)
    assert r.area_m2 == pytest.approx(8.0)
    assert r.contains(0.5, 1.0)
    assert not r.contains(5.0, 0.0)


# ======================================================================= layout presets


def test_2x2_plus_1_large_has_five_ports_tiling_without_gap():
    s = 0.05
    regions = layout_preset_2x2_plus_1_large(s, gap_m=0.0)
    assert len(regions) == 5
    small = regions[:4]
    large = regions[4]
    # the 2x2 block occupies x in [-2s, 0], y in [-s, s]; the large port x in [0, 2s], y in [-s, s]
    assert min(p.x_min_m for p in small) == pytest.approx(-2 * s)
    assert max(p.x_max_m for p in small) == pytest.approx(0.0)
    assert large.x_min_m == pytest.approx(0.0)
    assert large.x_max_m == pytest.approx(2 * s)
    assert large.area_m2 == pytest.approx((2 * s) ** 2)
    for p in small:
        assert p.area_m2 == pytest.approx(s * s)
    # total board tiles exactly (no gaps, no overlaps): total area == sum of port areas
    total_area = (4 * s) * (2 * s)
    assert sum(p.area_m2 for p in regions) == pytest.approx(total_area)


def test_2x2_plus_1_large_feed_defaults_to_each_ports_own_centre():
    regions = layout_preset_2x2_plus_1_large(0.05)
    for p in regions:
        assert p.feed_x_m == pytest.approx((p.x_min_m + p.x_max_m) / 2.0)
        assert p.feed_y_m == pytest.approx((p.y_min_m + p.y_max_m) / 2.0)


def test_large_port_feed_stays_centred_on_the_2x2_block_with_a_gap():
    """Regression (PHASE16 follow-ups): port 5 (the large port) used to inherit the small ports'
    own row-gap shift, pulling its feed off-centre by gap/2 even though it already spans the full
    height of the 2x2 block centred on the origin."""
    s, gap = 0.05, 0.0279
    regions = layout_preset_2x2_plus_1_large(s, gap_m=gap)
    large = regions[4]
    assert large.feed_y_m == pytest.approx(0.0)
    assert large.feed_x_m == pytest.approx((large.x_min_m + large.x_max_m) / 2.0)
    small_top_y = max(p.y_max_m for p in regions[:4])
    small_bottom_y = min(p.y_min_m for p in regions[:4])
    assert (small_top_y + small_bottom_y) / 2.0 == pytest.approx(0.0)


def test_gap_keeps_each_active_side_and_adds_the_gap_outside_it():
    s = 0.05
    gap = 0.004
    no_gap = layout_preset_2x2_plus_1_large(s, gap_m=0.0)
    with_gap = layout_preset_2x2_plus_1_large(s, gap_m=gap)
    for a, b in zip(no_gap, with_gap, strict=True):
        assert b.width_m == pytest.approx(a.width_m)
        assert b.height_m == pytest.approx(a.height_m)
    # the two columns of small ports and the large port are each separated by exactly one gap
    small_left, small_right = sorted(with_gap[:2], key=lambda p: p.x_min_m)[0], with_gap[1]
    assert with_gap[4].x_min_m - small_right.x_max_m == pytest.approx(gap)
    assert small_right.x_min_m - small_left.x_max_m == pytest.approx(gap)
    # the board is centred on the origin
    assert min(p.x_min_m for p in with_gap) == pytest.approx(-max(p.x_max_m for p in with_gap))


def test_1xn_gap_is_outside_each_active_side():
    regions = layout_preset_1xn(3, 0.12, gap_m=0.028)
    for p in regions:
        assert p.width_m == pytest.approx(0.12)
    for a, b in zip(regions, regions[1:], strict=False):
        assert b.x_min_m - a.x_max_m == pytest.approx(0.028)


def test_1xn_preset_builds_n_equal_ports_in_a_row():
    regions = layout_preset_1xn(3, 0.02)
    assert len(regions) == 3
    for p in regions:
        assert p.width_m == pytest.approx(0.02)
        assert p.height_m == pytest.approx(0.02)
    # contiguous, centred on the origin
    assert regions[0].x_min_m == pytest.approx(-0.03)
    assert regions[-1].x_max_m == pytest.approx(0.03)
    assert regions[0].x_max_m == pytest.approx(regions[1].x_min_m)


def test_custom_layout_defaults_feed_to_centre_and_accepts_explicit_feed():
    regions = custom_layout([(-1.0, 1.0, -1.0, 1.0), (0.0, 2.0, 0.0, 2.0, 0.3, 0.3)])
    assert regions[0].feed_x_m == pytest.approx(0.0)
    assert regions[1].feed_x_m == pytest.approx(0.3)
    assert regions[1].feed_y_m == pytest.approx(0.3)


def test_custom_layout_rejects_wrong_tuple_length():
    with pytest.raises(ValueError):
        custom_layout([(0.0, 1.0, 0.0)])


# ======================================================================= isolation estimate


def test_isolation_estimate_flags_touching_ports_and_not_distant_ones():
    k_swg = 2.0 * math.pi / 0.01  # lambda_sw = 10 mm
    touching = (
        PortRegion(0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
        PortRegion(1, 1.0, 2.0, 0.0, 1.0, 1.5, 0.5),  # shares the x=1 edge: gap = 0
    )
    flags = isolation_estimate(touching, k_swg, min_lambda_sw=1.0)
    assert len(flags) == 1
    assert isinstance(flags[0], IsolationFlag)
    assert flags[0].gap_m == pytest.approx(0.0)
    assert flags[0].flagged
    assert flags[0].guard_row_suggestion_m == pytest.approx(0.01)  # needs a full lambda_sw of gap

    far = (
        PortRegion(0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
        PortRegion(1, 2.0, 3.0, 0.0, 1.0, 2.5, 0.5),  # 1 m gap >> 10 mm lambda_sw
    )
    flags_far = isolation_estimate(far, k_swg, min_lambda_sw=1.0)
    assert not flags_far[0].flagged
    assert flags_far[0].guard_row_suggestion_m is None


def test_isolation_estimate_does_not_flag_a_gap_at_exactly_the_minimum():
    """Regression (PHASE16 follow-ups): the app's own default gap is exactly 1.00 lambda_sw
    (rounded for display), and a strict `<` comparison flagged it anyway on floating-point
    rounding in the lambda_sw conversion."""
    k_swg = 2.0 * math.pi / 0.01  # lambda_sw = 10 mm
    lambda_sw = 2.0 * math.pi / k_swg
    exact = (
        PortRegion(0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
        PortRegion(1, 1.0 + lambda_sw, 2.0 + lambda_sw, 0.0, 1.0, 1.5 + lambda_sw, 0.5),
    )
    flags = isolation_estimate(exact, k_swg, min_lambda_sw=1.0)
    assert flags[0].gap_lambda_sw == pytest.approx(1.0)
    assert not flags[0].flagged

    just_under = (
        PortRegion(0, 0.0, 1.0, 0.0, 1.0, 0.5, 0.5),
        PortRegion(1, 1.0 + 0.5 * lambda_sw, 2.0 + 0.5 * lambda_sw, 0.0, 1.0, 1.5 + 0.5 * lambda_sw, 0.5),
    )
    flags_under = isolation_estimate(just_under, k_swg, min_lambda_sw=1.0)
    assert flags_under[0].flagged


def test_isolation_estimate_requires_positive_k_swg():
    regions = layout_preset_1xn(2, 0.02)
    with pytest.raises(ValueError):
        isolation_estimate(regions, None)
    with pytest.raises(ValueError):
        isolation_estimate(regions, 0.0)


# ======================================================================= Phase 21: grid, snap, overlap


def test_grid_2x2_matches_the_reference_board_port_order_and_geometry():
    from src.hologram.ports import layout_preset_grid

    s, g = 0.3008, 0.028
    regions = layout_preset_grid(2, 2, s, g)
    centre = (s + g) / 2.0
    assert [(r.feed_x_m > 0, r.feed_y_m > 0) for r in regions] == [(True, True), (False, True), (False, False), (True, False)]
    for r in regions:
        assert abs(r.feed_x_m) == pytest.approx(centre) and abs(r.feed_y_m) == pytest.approx(centre)
        assert r.width_m == pytest.approx(s) and r.height_m == pytest.approx(s)
        assert r.feed_x_m == pytest.approx((r.x_min_m + r.x_max_m) / 2.0)
    assert regions[0].x_min_m == pytest.approx(g / 2.0)  # 14 mm for a 28 mm gap
    assert [r.port_id for r in regions] == [0, 1, 2, 3]


def test_grid_is_centred_serpentine_and_gap_free_tiles():
    from src.hologram.ports import layout_preset_grid

    regions = layout_preset_grid(2, 3, 0.1, 0.0)
    assert len(regions) == 6
    assert min(r.x_min_m for r in regions) == pytest.approx(-0.15) and max(r.x_max_m for r in regions) == pytest.approx(0.15)
    assert min(r.y_min_m for r in regions) == pytest.approx(-0.1) and max(r.y_max_m for r in regions) == pytest.approx(0.1)
    # top row right->left, second row left->right
    assert [round(r.feed_x_m, 3) for r in regions] == [0.1, 0.0, -0.1, -0.1, 0.0, 0.1]
    assert [round(r.feed_y_m, 3) for r in regions] == [0.05] * 3 + [-0.05] * 3


def test_grid_1x1_and_validation():
    from src.hologram.ports import layout_preset_grid

    assert len(layout_preset_grid(1, 1, 0.2)) == 1
    for bad in ((0, 2, 0.1, 0.0), (2, 2, 0.0, 0.0), (2, 2, 0.1, -1.0), (2.5, 2, 0.1, 0.0)):
        with pytest.raises(ValueError):
            layout_preset_grid(*bad)


def test_snap_side_to_cells():
    from src.hologram.ports import snap_side_to_cells

    side, n = snap_side_to_cells(0.300, 0.0064)
    assert n == 47 and side == pytest.approx(0.3008)
    assert snap_side_to_cells(0.0001, 0.0064)[1] == 1


def test_check_no_overlap():
    from src.hologram.ports import check_no_overlap

    check_no_overlap(custom_layout([(0, 1, 0, 1), (1, 2, 0, 1)]))  # touching is fine
    with pytest.raises(ValueError, match="ports 1 and 2 overlap"):
        check_no_overlap(custom_layout([(0, 1, 0, 1), (0.5, 2, 0.5, 2)]))


def test_build_port_regions_grid_snap_and_custom_validation():
    from src.app.routers.hmsa import _build_port_regions

    base = {"preset": "grid", "rows": 2, "cols": 2, "cell_side_mm": 300.0, "gap_mm": 0.0, "snap_side_to_cells": True,
            "small_side_mm": None, "n": None, "regions_mm": None}
    snapped = _build_port_regions(base, 6.4)
    assert snapped[0].width_m == pytest.approx(0.3008)
    unsnapped = _build_port_regions({**base, "snap_side_to_cells": False}, 6.4)
    assert unsnapped[0].width_m == pytest.approx(0.3)
    overlap = {**base, "preset": "custom", "regions_mm": [[0, 10, 0, 10], [5, 15, 5, 15]]}
    with pytest.raises(Exception, match="overlap"):
        _build_port_regions(overlap, 6.4)
    with pytest.raises(Exception, match="rows, cols"):
        _build_port_regions({**base, "rows": None}, 6.4)
