"""Shared default design inputs for Phase 19 (D = 150 mm octagon, p = 5 mm, F = 1.2 D, 28 GHz, +-30 deg)."""

from __future__ import annotations

import math

from .feed import q_for_edge_taper
from .layout import LensLayout

C0 = 299_792_458.0
FREQ_HZ = 28e9
D_FLATS_M = 0.150
PITCH_M = 0.005
F_OVER_D = 1.2
BEAM_DEG = 30.0
EDGE_TAPER_DB = 10.0
WINDOW_LAMBDA_OVER_D = 0.6


def make_layout() -> LensLayout:
    return LensLayout(D_FLATS_M, PITCH_M, F_OVER_D * D_FLATS_M)


def q_feed(layout: LensLayout, edge_taper_db: float = EDGE_TAPER_DB) -> float:
    return q_for_edge_taper(math.degrees(math.atan(layout.d_flats_m / 2 / layout.focal_m)), edge_taper_db)


def dual_beam_phase_map():
    """(x, y, theta_cell_deg, illumination weight, target T phase base in deg, without psi0) for the default lens.

    Same synthesis as the final design (symmetric, pointing-selected GS): arg T_c = psi_c + k0 r_c + psi0.
    """
    from .pipeline import DesignConfig, synthesize

    lay = make_layout()
    cfg = DesignConfig(FREQ_HZ, lay, BEAM_DEG, EDGE_TAPER_DB, None, 0.6e-3, 230.0, 0.8044, "")
    syn = synthesize(cfg)
    return syn["x"], syn["y"], syn["theta"], syn["weight"], syn["base_phase_deg"]
