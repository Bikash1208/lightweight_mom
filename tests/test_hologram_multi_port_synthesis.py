"""Phase 15 Part C (gap-fill): multi-port sub-aperture synthesis -- N independent single-aperture
syntheses, one per `ports.PortRegion`, reusing `hologram.hmsa_synthesis.synthesize_formulation_b`
unmodified."""

import math

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.hologram import hmsa_synthesis as hs
from src.hologram.multi_port_synthesis import (
    combined_board_lw_grid_m,
    port_native_grid,
    port_reference_field,
    synthesize_ports,
)
from src.hologram.ports import layout_preset_1xn, layout_preset_2x2_plus_1_large
from src.phase8b.object_wave import Beam
from src.phase8b.reference_wave import monopole_reference_field
from src.phase8e.reference_field_input import phase_only_reference_field

FREQ = 10e9
K_SWG = tm0_beta(FREQ, 2.2, 1.575e-3)
PERIOD_M = 2e-3

# A dense phase-only candidate library (5deg steps) so Formulation B can match any required phase
# to within ~2.5deg -- good enough for a direction/steering regression test.
_CANDIDATES = [(f"c{i}", 2.0, 2.0, "lib", np.exp(1j * math.radians(i * 5.0))) for i in range(72)]


def _per_port(ref, beams, feed_xy_m):
    return hs.synthesize_formulation_b(ref, beams, FREQ, _CANDIDATES, period_m=PERIOD_M, feed_xy_m=feed_xy_m)


# ============================================================================ port_native_grid


def test_port_native_grid_spans_the_region():
    regions = layout_preset_1xn(1, 0.02)
    x_m, y_m = port_native_grid(regions[0], resolution_m=0.002)
    assert x_m.min() == pytest.approx(regions[0].x_min_m, abs=1e-9)
    assert x_m.max() == pytest.approx(regions[0].x_max_m, abs=1e-9)
    assert y_m.min() == pytest.approx(regions[0].y_min_m, abs=1e-9)
    assert y_m.max() == pytest.approx(regions[0].y_max_m, abs=1e-9)


# ========================================================================= port_reference_field


def test_port_reference_field_with_no_upload_is_the_analytic_monopole_phase_about_the_feed():
    from dataclasses import replace

    region = replace(layout_preset_1xn(1, 0.02)[0], feed_x_m=0.003)  # an off-centre feed
    ref = port_reference_field(region, resolution_m=0.002, k_swg=K_SWG)
    expected = monopole_reference_field(ref.x_m - region.feed_x_m, ref.y_m - region.feed_y_m, K_SWG)
    assert np.allclose(ref.Ex, np.exp(1j * np.angle(expected)))


def test_analytic_single_port_reproduces_a_phase_only_feed_design():
    """Regression: the analytic monopole's rho^-1/2 amplitude made Formulation A report 28.5% of a
    full-board port unreachable, while a phase-only field carrying the same phase reports 0%
    (the uploaded .fld is phase-only). The analytic port must match that phase-only design."""
    from src.dispersion.hfss_table import DispersionSample, DispersionTable

    table = DispersionTable(
        source_path="synthetic", frequency_hz=FREQ,
        samples=tuple(DispersionSample(geometry_param=g, beta_rad_per_m=None, x_s_ohm=x)
                      for g, x in [(1.0, 50.0), (2.0, 100.0), (3.0, 150.0), (4.0, 200.0)]),
    )
    region = layout_preset_1xn(1, 0.296)[0]  # a 37 x 37 cell board at 8 mm
    beams = [Beam(theta_deg=20.0, phi_deg=0.0, weight=1.0 + 0j)]
    analytic = port_reference_field(region, resolution_m=0.008, k_swg=K_SWG)
    phase_only = phase_only_reference_field(
        analytic.x_m, analytic.y_m, np.degrees(np.angle(monopole_reference_field(analytic.x_m, analytic.y_m, K_SWG)))
    )
    design_analytic = hs.synthesize_formulation_a(
        analytic, beams, FREQ, 100.0, 0.3, table, 0.008, rescale_to_beta_bar=False, feed_xy_m=(0.0, 0.0)
    )
    design_phase_only = hs.synthesize_formulation_a(
        phase_only, beams, FREQ, 100.0, 0.3, table, 0.008, rescale_to_beta_bar=False, feed_xy_m=(0.0, 0.0)
    )
    assert design_analytic.unreachable_fraction <= design_phase_only.unreachable_fraction + 0.02
    assert design_analytic.unreachable_fraction < 0.05


# =============================================================================== synthesize_ports


def test_synthesize_ports_runs_one_independent_design_per_port():
    regions = layout_preset_1xn(2, 0.02)
    beams = {
        0: [Beam(theta_deg=20.0, phi_deg=0.0)],
        1: [Beam(theta_deg=30.0, phi_deg=90.0)],
    }
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)
    assert len(multi.ports) == 2
    for port_result, region in zip(multi.ports, regions, strict=True):
        assert port_result.region is region
        x_m, y_m = port_result.design.x_m, port_result.design.y_m
        assert np.all(x_m >= region.x_min_m - 1e-9)
        assert np.all(x_m <= region.x_max_m + 1e-9)
        assert np.all(y_m >= region.y_min_m - 1e-9)
        assert np.all(y_m <= region.y_max_m + 1e-9)


def test_synthesize_ports_requires_beams_for_every_port():
    regions = layout_preset_1xn(2, 0.02)
    with pytest.raises(ValueError):
        synthesize_ports(regions, {0: [Beam(theta_deg=0.0, phi_deg=0.0)]}, _per_port, period_m=PERIOD_M, k_swg=K_SWG)


def test_synthesize_ports_computes_isolation_flags_when_k_swg_given():
    regions = layout_preset_2x2_plus_1_large(0.02)
    beams = {i: [Beam(theta_deg=10.0 * (i + 1), phi_deg=0.0)] for i in range(5)}
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)
    # 5 ports -> C(5,2) = 10 pairwise isolation flags
    assert len(multi.isolation_flags) == 10


def test_synthesize_ports_skips_isolation_without_k_swg():
    from src.phase8e.reference_field_input import phase_only_reference_field

    regions = layout_preset_1xn(2, 0.02)
    beams = {0: [Beam(theta_deg=0.0, phi_deg=0.0)], 1: [Beam(theta_deg=0.0, phi_deg=0.0)]}
    native = {i: port_native_grid(r, resolution_m=0.002) for i, r in enumerate(regions)}
    uploaded = {i: phase_only_reference_field(x, y, np.zeros_like(x)) for i, (x, y) in native.items()}
    multi = synthesize_ports(
        regions, beams, _per_port, period_m=PERIOD_M, k_swg=None, uploaded_fields=uploaded
    )
    assert multi.isolation_flags == ()


# ======================================================================= combined_board_lw_grid_m


def _cell_m(cell):
    if cell.geometry_param is None:
        return None
    return cell.geometry_param * 1e-3, (cell.geometry_param_w or cell.geometry_param) * 1e-3


def test_combined_board_lw_grid_covers_the_whole_board_rectangle():
    regions = layout_preset_1xn(2, 0.02)
    beams = {0: [Beam(theta_deg=0.0, phi_deg=0.0)], 1: [Beam(theta_deg=0.0, phi_deg=0.0)]}
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)
    board = combined_board_lw_grid_m(multi, _cell_m)
    assert board.width_m == pytest.approx(0.04)
    assert board.height_m == pytest.approx(0.02)
    n_y, n_x = len(board.lw_grid_m), len(board.lw_grid_m[0])
    assert n_x == round(board.width_m / PERIOD_M)
    assert n_y == round(board.height_m / PERIOD_M)
    # every cell is populated with a real (non-None) (length_m, width_m) pair -- the 1xN preset
    # tiles exactly, with gap_m=0, so no combined-grid cell should fall outside every port
    for row in board.lw_grid_m:
        for length_m, width_m_ in row:
            assert length_m > 0.0 and width_m_ > 0.0
    # cell_center_m reconstructs a cell's true board position, matching this module's own lattice
    cx, cy = board.cell_center_m(0, 0)
    assert cx == pytest.approx(board.x_min_m + board.period_x_m / 2.0)
    assert cy == pytest.approx(board.y_min_m + board.period_y_m / 2.0)


def test_combined_board_lw_grid_marks_gap_cells_as_no_patch():
    regions = layout_preset_1xn(2, 0.02, gap_m=0.006)  # 3 periods of guard gap between the ports
    beams = {0: [Beam(theta_deg=0.0, phi_deg=0.0)], 1: [Beam(theta_deg=0.0, phi_deg=0.0)]}
    multi = synthesize_ports(regions, beams, _per_port, period_m=PERIOD_M, k_swg=K_SWG)
    board = combined_board_lw_grid_m(multi, _cell_m)
    assert any(cell == (0.0, 0.0) for row in board.lw_grid_m for cell in row)
