"""The STEP writer centres the board on the origin, while ports and the combined cell grid are indexed
from the board's own corner. Feed markers and port rectangles must be expressed in the STEP frame, and
each feed must lie inside its own port (the check applied to every export that writes feed markers)."""

import numpy as np

from src.dispersion.analytic_tm0 import tm0_beta
from src.hologram import hmsa_synthesis as hs
from src.hologram.multi_port_synthesis import (
    combined_board_lw_grid_m,
    export_frame_feeds_and_rects_mm,
    synthesize_ports,
)
from src.hologram.ports import custom_layout
from src.phase8b.object_wave import Beam

FREQ = 10e9
PERIOD_M = 4e-3
K_SWG = tm0_beta(FREQ, 2.2, 1.575e-3)
CANDIDATES = [(f"c{i}", 2.0, 2.0, "lib", np.exp(1j * np.radians(i * 5.0))) for i in range(72)]


def _per_port(ref, beams, feed_xy_m):
    return hs.synthesize_formulation_b(ref, beams, FREQ, CANDIDATES, period_m=PERIOD_M, feed_xy_m=feed_xy_m)


def _cell_m(cell):
    if cell.geometry_param is None:
        return None
    return cell.geometry_param * 1e-3, (cell.geometry_param_w or cell.geometry_param) * 1e-3


def _board(rects):
    regions = custom_layout(rects)
    beams = {r.port_id: [Beam(theta_deg=0.0, phi_deg=0.0)] for r in regions}
    multi = synthesize_ports(regions, beams, _per_port, PERIOD_M, k_swg=K_SWG)
    return multi, combined_board_lw_grid_m(multi, _cell_m)


def _inside(feed, rect):
    x, y = feed
    x0, x1, y0, y1 = rect
    return x0 <= x <= x1 and y0 <= y <= y1


def test_every_feed_marker_lies_inside_its_own_port_in_the_step_frame():
    # a board that does not start at the origin: its corner frame and its centred STEP frame differ
    multi, board = _board([(0.0, 0.040, 0.0, 0.040), (0.048, 0.088, 0.0, 0.040), (0.0, 0.088, 0.048, 0.088)])
    feeds, rects = export_frame_feeds_and_rects_mm(multi, board)
    for feed, rect in zip(feeds, rects, strict=True):
        assert _inside(feed, rect)


def test_corner_frame_feeds_would_fall_outside_the_centred_board():
    multi, board = _board([(0.0, 0.040, 0.0, 0.040), (0.048, 0.088, 0.0, 0.040)])
    corner_feeds = [(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports]
    _feeds, rects = export_frame_feeds_and_rects_mm(multi, board)
    assert not all(_inside(f, r) for f, r in zip(corner_feeds, rects, strict=True))


def test_step_frame_is_centred_on_the_origin():
    multi, board = _board([(0.0, 0.040, 0.0, 0.040), (0.048, 0.088, 0.0, 0.040)])
    _feeds, rects = export_frame_feeds_and_rects_mm(multi, board)
    x_min = min(r[0] for r in rects)
    x_max = max(r[1] for r in rects)
    assert abs(0.5 * (x_min + x_max)) < PERIOD_M * 1e3 / 2.0 + 1e-6
