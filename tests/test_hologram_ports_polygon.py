"""Phase 23 A: polygon (hexagon) port regions -- geometry against analytic values, unique cell
ownership, isolation/overlap on polygons, and that a plain rectangle `PortRegion` is untouched."""

import math

import numpy as np
import pytest

from src.hologram.multi_port_synthesis import mask_design_to_region, region_cell_mask
from src.hologram.ports import (
    PortRegion,
    check_no_overlap,
    isolation_estimate,
    layout_preset_grid,
    layout_preset_hex4,
    overlap_area_m2,
    polygon_area_m2,
    region_gap_m,
    regular_hexagon_vertices,
)

SIDE = 0.150
PERIOD = 0.0064


def test_rectangle_region_is_unchanged():
    r = PortRegion(3, -0.1, 0.2, -0.05, 0.05, 0.0, 0.0)
    assert r.vertices_m is None and not r.is_polygon
    assert r.area_m2 == r.width_m * r.height_m == pytest.approx(0.3 * 0.1)
    assert r.contains(0.2, 0.05) and not r.contains(0.2000001, 0.0)
    assert r.outline_m == ((-0.1, -0.05), (0.2, -0.05), (0.2, 0.05), (-0.1, 0.05))
    # presets still build plain rectangles
    assert all(not p.is_polygon for p in layout_preset_grid(2, 2, 0.3, gap_m=0.0))


def test_hexagon_area_and_extents_are_analytic():
    verts = regular_hexagon_vertices(0.0, 0.0, SIDE, "pointy")
    assert polygon_area_m2(verts) == pytest.approx(1.5 * math.sqrt(3.0) * SIDE**2, rel=1e-12)
    assert polygon_area_m2(verts) * 1e6 == pytest.approx(58456.7, abs=0.1)  # the brief's 58,457 mm^2
    xs, ys = [v[0] for v in verts], [v[1] for v in verts]
    assert max(xs) - min(xs) == pytest.approx(math.sqrt(3.0) * SIDE)  # flat-to-flat 259.8 mm
    assert max(ys) - min(ys) == pytest.approx(2.0 * SIDE)  # corner-to-corner 300 mm
    flat = regular_hexagon_vertices(0.0, 0.0, SIDE, "flat")
    assert max(v[0] for v in flat) - min(v[0] for v in flat) == pytest.approx(2.0 * SIDE)
    assert max(v[1] for v in flat) - min(v[1] for v in flat) == pytest.approx(math.sqrt(3.0) * SIDE)


def test_hex_membership_and_cell_count_match_area():
    region = layout_preset_hex4(side_m=SIDE)[0]
    cx, cy = region.feed_x_m, region.feed_y_m
    assert region.contains(cx, cy)
    assert region.contains(cx, cy + SIDE)  # a corner is on the boundary (inclusive)
    assert region.contains(cx + 0.5 * math.sqrt(3.0) * SIDE, cy)  # mid-edge
    assert not region.contains(cx + 0.5 * math.sqrt(3.0) * SIDE + 1e-4, cy)
    assert not region.contains(region.x_max_m, region.y_max_m)  # a corner of the bounding box
    # cell centres of a 6.4 mm lattice inside the hexagon ~ area / period^2
    xs = np.arange(region.x_min_m + PERIOD / 2, region.x_max_m, PERIOD)
    ys = np.arange(region.y_min_m + PERIOD / 2, region.y_max_m, PERIOD)
    gx, gy = np.meshgrid(xs, ys)
    n = int(region.contains_array(gx, gy).sum())
    assert n == pytest.approx(region.area_m2 / PERIOD**2, rel=0.01)
    # the array form agrees with the scalar form
    assert [bool(region.contains(x, y)) for x, y in zip(gx.ravel()[:50], gy.ravel()[:50], strict=True)] == list(
        region.contains_array(gx, gy).ravel()[:50]
    )


def test_hex4_cluster_geometry():
    regions = layout_preset_hex4(side_m=SIDE)
    assert len(regions) == 4
    # centred on the origin, four equal hexagons
    assert sum(r.feed_x_m for r in regions) == pytest.approx(0.0, abs=1e-12)
    assert sum(r.feed_y_m for r in regions) == pytest.approx(0.0, abs=1e-12)
    assert {round(r.area_m2, 12) for r in regions} == {round(1.5 * math.sqrt(3.0) * SIDE**2, 12)}
    # default beam directions point away from the centre: port k's own direction is within 15 deg of 45+90k
    for k, r in enumerate(regions):
        direction = math.degrees(math.atan2(r.feed_y_m, r.feed_x_m)) % 360.0
        assert abs((direction - (45.0 + 90.0 * k) + 180.0) % 360.0 - 180.0) <= 15.0 + 1e-9
    check_no_overlap(regions)  # touching edges are allowed
    # every touching pair has gap 0; only the two ends of the rows face corner to corner: the centres are
    # 3*side apart along a corner direction, minus two circumradii = one side
    gaps = sorted(f.gap_m for f in isolation_estimate(regions, 100.0))
    assert gaps[:5] == [0.0] * 5
    assert gaps[5] == pytest.approx(SIDE)
    assert all(f.gap_lambda_sw == 0.0 for f in isolation_estimate(regions, 100.0) if f.gap_m == 0.0)


def test_hex4_alternative_size_input_gap_and_orientation():
    a = layout_preset_hex4(side_m=SIDE)
    b = layout_preset_hex4(flat_to_flat_m=math.sqrt(3.0) * SIDE)
    assert [r.vertices_m for r in a] == pytest.approx([r.vertices_m for r in b], abs=1e-12)
    with pytest.raises(ValueError):
        layout_preset_hex4()
    with pytest.raises(ValueError):
        layout_preset_hex4(side_m=SIDE, flat_to_flat_m=0.26)
    gapped = layout_preset_hex4(side_m=SIDE, gap_m=0.010)
    check_no_overlap(gapped)
    assert min(f.gap_m for f in isolation_estimate(gapped, 100.0)) == pytest.approx(0.010, abs=1e-12)
    flat = layout_preset_hex4(side_m=SIDE, orientation="flat")
    check_no_overlap(flat)
    assert max(v[0] for v in flat[0].vertices_m) - min(v[0] for v in flat[0].vertices_m) == pytest.approx(2 * SIDE)
    assert min(f.gap_m for f in isolation_estimate(flat, 100.0)) == 0.0
    for k, r in enumerate(flat):  # port 1 is still the one nearest 45 deg
        direction = math.degrees(math.atan2(r.feed_y_m, r.feed_x_m)) % 360.0
        assert abs((direction - (45.0 + 90.0 * k) + 180.0) % 360.0 - 180.0) <= 15.0 + 1e-9


def test_lattice_snapped_boxes_share_one_lattice():
    regions = layout_preset_hex4(side_m=SIDE, period_m=PERIOD)
    for r in regions:
        for v in (r.x_min_m, r.x_max_m, r.y_min_m, r.y_max_m):
            assert v / PERIOD == pytest.approx(round(v / PERIOD), abs=1e-9)
        assert all(r.x_min_m <= x <= r.x_max_m and r.y_min_m <= y <= r.y_max_m for x, y in r.vertices_m)
    # the polygon itself is unchanged by the snapping
    plain = layout_preset_hex4(side_m=SIDE)
    assert [r.vertices_m for r in regions] == [r.vertices_m for r in plain]


def test_overlap_is_detected_and_area_is_exact():
    a, b = layout_preset_hex4(side_m=SIDE)[:2]
    assert overlap_area_m2(a, b) == pytest.approx(0.0, abs=1e-12)
    shifted = PortRegion(
        9,
        b.x_min_m + 0.05,
        b.x_max_m + 0.05,
        b.y_min_m,
        b.y_max_m,
        b.feed_x_m + 0.05,
        b.feed_y_m,
        vertices_m=tuple((x + 0.05, y) for x, y in b.vertices_m),
    )
    assert overlap_area_m2(a, shifted) > 1e-3
    assert region_gap_m(a, shifted) == 0.0
    with pytest.raises(ValueError, match="overlap"):
        check_no_overlap([a, shifted])
    # a hexagon against itself
    assert overlap_area_m2(a, a) == pytest.approx(a.area_m2, rel=1e-9)


def test_invalid_polygons_are_rejected():
    with pytest.raises(ValueError, match="convex"):
        PortRegion(0, 0, 1, 0, 1, 0.5, 0.5, vertices_m=((0, 0), (0, 1), (1, 1), (1, 0)))  # clockwise
    with pytest.raises(ValueError, match="convex"):
        PortRegion(0, 0, 1, 0, 1, 0.5, 0.5, vertices_m=((0, 0), (1, 0), (0.2, 0.2), (1, 1), (0, 1)))  # dent
    with pytest.raises(ValueError, match="bounding box"):
        PortRegion(0, 0, 1, 0, 1, 0.5, 0.5, vertices_m=((0, 0), (2, 0), (1, 1)))
    with pytest.raises(ValueError):
        PortRegion(0, 0, 1, 0, 1, 0.5, 0.5, vertices_m=((0, 0), (1, 0)))


def test_each_cell_has_at_most_one_owner_and_a_gapless_cluster_is_fully_covered():
    regions = layout_preset_hex4(side_m=SIDE, period_m=PERIOD)
    x0, x1 = min(r.x_min_m for r in regions), max(r.x_max_m for r in regions)
    y0, y1 = min(r.y_min_m for r in regions), max(r.y_max_m for r in regions)
    xs = x0 + (np.arange(round((x1 - x0) / PERIOD)) + 0.5) * PERIOD
    ys = y0 + (np.arange(round((y1 - y0) / PERIOD)) + 0.5) * PERIOD
    gx, gy = np.meshgrid(xs, ys)
    owners = np.zeros(gx.shape, dtype=int)
    for k, r in enumerate(regions):
        owners += region_cell_mask(r, gx, gy, earlier_regions=regions[:k])
    assert owners.max() == 1
    # a cell centre inside the union of the hexagons is owned (the union has no gap at gap_m = 0)
    union = np.zeros(gx.shape, dtype=bool)
    for r in regions:
        union |= r.contains_array(gx, gy)
    assert np.array_equal(owners.astype(bool), union)


def test_mask_design_to_region_zeroes_and_drops_outside_cells():
    from dataclasses import dataclass

    from src.hologram.hmsa_synthesis import GeometryCell

    @dataclass(frozen=True)
    class _Design:
        cells: tuple
        aperture_field: object

    ny, nx = 3, 4
    cells = tuple(
        GeometryCell(i, j, float(j), float(i), 0.0, 0.0, 0.0, 0.0, True) for i in range(ny) for j in range(nx)
    )
    d = _Design(cells, np.ones((ny, nx), dtype=complex))
    mask = np.zeros((ny, nx), dtype=bool)
    mask[1, 1:3] = True
    out = mask_design_to_region(d, mask)
    assert [(c.i, c.j) for c in out.cells] == [(1, 1), (1, 2)]
    assert out.aperture_field.sum() == 2 and out.aperture_field[0, 0] == 0
