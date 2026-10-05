"""The OAM charge sign rule for circularly polarized vortices with a monopole feed (this code's
convention: RHCP needs l<0, LHCP needs l>0), plus the app-facing warning for the wrong sign.

The physics test uses an ideal synthetic phase library (every phase available) and the analytic
monopole phase, so it runs in the fast suite without a MoM solve."""

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.green_function import LayeredGreenFunction
from src.hologram.hmsa_synthesis import synthesize_formulation_b
from src.hologram.polarization_compensation import oam_sign_mismatch
from src.labels import oam_sign_mismatch_label
from src.materials import Substrate
from src.pattern.vector_aperture_pattern import (
    _vector_far_field,
    aperture_current_from_scalar,
    compute_vector_pattern,
    measure_vortex_winding,
)
from src.phase8b.object_wave import Beam
from src.phase8b.reference_wave import monopole_reference_field
from src.phase8e.reference_field_input import phase_only_reference_field

FREQ = 10.54e9
PERIOD = 8e-3
SUB = Substrate("user", 2.2, 0.0009, 1.575e-3)
K_SWG = tm0_beta(FREQ, 2.2, 1.575e-3)


@pytest.mark.parametrize(
    ("pol", "charge", "flagged"),
    [("rhcp", -1, False), ("rhcp", -2, False), ("rhcp", 2, True), ("rhcp", 1, True),
     ("lhcp", 2, False), ("lhcp", -1, True), ("lhcp", 0, False), ("natural", 2, False)],
)
def test_sign_rule_flags_only_the_opposite_sign(pol, charge, flagged):
    assert oam_sign_mismatch(pol, charge) is flagged


def test_sign_mismatch_label_names_the_handedness_and_charge():
    assert oam_sign_mismatch_label("rhcp", 2) == "OAM_SIGN_MISMATCH(RHCP,l=+2)"


def _ideal_library():
    k = np.arange(200)
    phases = np.radians(360.0 * k / 200.0)
    lengths = 0.5 + 0.02 * k  # mm, a monotone label for each ideal phase
    return [(f"ideal{i}", float(lengths[i]), float(lengths[i]), "ideal", complex(np.cos(p), np.sin(p)))
            for i, p in enumerate(phases)]


def _ring_xpd_db(pat_fn, theta_b_deg, phi_b_deg, desired, radius_uv=0.02, n=48):
    """Co- over cross-polar power, averaged on a small (u,v) ring around the beam direction."""
    u_b = np.sin(np.radians(theta_b_deg)) * np.cos(np.radians(phi_b_deg))
    v_b = np.sin(np.radians(theta_b_deg)) * np.sin(np.radians(phi_b_deg))
    alpha = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    u = u_b + radius_uv * np.cos(alpha)
    v = v_b + radius_uv * np.sin(alpha)
    theta = np.degrees(np.arcsin(np.clip(np.hypot(u, v), 0.0, 1.0)))
    phi = np.degrees(np.arctan2(v, u)) % 360.0
    e_theta, e_phi = pat_fn(theta, phi)
    from src.pattern.vector_aperture_pattern import circular_components

    e_r, e_l = circular_components(e_theta, e_phi)
    co, cross = (e_l, e_r) if desired == "rhcp" else (e_r, e_l)
    return float(10.0 * np.log10(np.mean(np.abs(co) ** 2) / np.mean(np.abs(cross) ** 2)))


def _vortex_design_pattern(desired, charge):
    xs = (np.arange(37) - 18) * PERIOD + PERIOD / 2.0  # cell centres offset from the feed
    X, Y = np.meshgrid(xs, xs)
    reference = phase_only_reference_field(X, Y, np.degrees(np.angle(monopole_reference_field(X, Y, K_SWG))))
    beam = Beam(theta_deg=20.0, phi_deg=0.0, weight=1.0 + 0j, desired_polarization=desired, oam_charge=charge)
    design = synthesize_formulation_b(reference, [beam], FREQ, _ideal_library(), PERIOD, feed_xy_m=(0.0, 0.0))
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    dx = float(design.x_m[0, 1] - design.x_m[0, 0])
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=(0.0, 0.0))
    pat = compute_vector_pattern(design.x_m, design.y_m, jx, jy, FREQ, green, dx, dx)
    peak = float(np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2))

    def pat_fn(theta, phi):
        return _vector_far_field(theta, phi, design.x_m, design.y_m, jx, jy, FREQ, green, dx, dx)

    return design, green, jx, jy, peak, pat, pat_fn


def test_rhcp_with_negative_charge_gives_a_clean_vortex():
    design, green, jx, jy, peak, _pat, pat_fn = _vortex_design_pattern("rhcp", -1)
    vortex = measure_vortex_winding(
        design.x_m, design.y_m, jx, jy, FREQ, green, PERIOD, PERIOD, 20.0, 0.0, "rhcp", pattern_peak_power=peak,
    )
    assert vortex.winding_number == pytest.approx(-1.0, abs=0.05)
    # The ideal synthetic library gives about -24.6 dB; the real-library check (slow test) enforces -25 dB.
    assert vortex.null_depth_db <= -20.0
    xpd_ring = _ring_xpd_db(pat_fn, 20.0, 0.0, "rhcp")
    assert xpd_ring >= 10.0, xpd_ring


def test_rhcp_with_positive_charge_two_is_flagged_and_is_not_a_clean_rhcp_vortex():
    assert oam_sign_mismatch("rhcp", 2)
    design, green, jx, jy, peak, _pat, _pat_fn = _vortex_design_pattern("rhcp", 2)
    vortex = measure_vortex_winding(
        design.x_m, design.y_m, jx, jy, FREQ, green, PERIOD, PERIOD, 20.0, 0.0, "rhcp", pattern_peak_power=peak,
    )
    assert abs(vortex.winding_number - 2.0) > 0.5 or vortex.null_depth_db > -25.0


@pytest.mark.slow
def test_rhcp_l_minus_1_on_the_real_square_library_meets_the_spec(tmp_path):
    """The submission figure's own path: the tracked example .fld, the fast-mesh square-patch library
    solved here at 10.54 GHz, Formulation B, RHCP l=-1 -- winding -1, on-axis null <= -25 dB."""
    from pathlib import Path

    from src.app.library_query import library_candidates_at_frequency
    from src.app.routers.hmsa import AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY, _l_fractions_of_period, _mask_feed_region
    from src.fields.fld_parser import read_fld
    from src.fields.reference_field_grid import (
        feed_position_from_nan_mask,
        fld_to_reference_field_grid,
        resample_to_cell_grid,
    )
    from src.library_sweep import build_square_library_cells, write_square_library_json
    from src.unit_cell_library import UnitCellLibrary

    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY["fast"]
    l_values = tuple(round(PERIOD * 1e3 * f, 6) for f in _l_fractions_of_period(mesh_cells))
    base = {
        "period_x": PERIOD, "period_y": PERIOD,
        "substrate": {"name": "user", "eps_r": 2.2, "tan_delta": 0.0009, "thickness": 1.575e-3, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": l_values[0] * 1e-3, "width": l_values[0] * 1e-3},
        "excitation": {"frequency": FREQ},
        "mesh": {"cells_x": mesh_cells, "cells_y": mesh_cells, "floquet_order": mesh_cells},
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }
    cells = build_square_library_cells(base, l_values, FREQ)
    library_path = tmp_path / "square_library.json"
    write_square_library_json(cells, library_path)
    candidates = library_candidates_at_frequency(UnitCellLibrary.from_json(library_path), "real", FREQ)

    fld = Path(__file__).resolve().parents[1] / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
    raw = fld_to_reference_field_grid(read_fld(fld))
    nan = np.isnan(raw.phase())
    masked, _ = _mask_feed_region(raw)
    x0, y0 = feed_position_from_nan_mask(raw.x_m, raw.y_m, nan)
    reference = resample_to_cell_grid(masked, PERIOD)

    beam = Beam(theta_deg=20.0, phi_deg=0.0, weight=1.0 + 0j, desired_polarization="rhcp", oam_charge=-1)
    design = synthesize_formulation_b(reference, [beam], FREQ, candidates, PERIOD, feed_xy_m=(float(x0), float(y0)))
    green = LayeredGreenFunction(SUB, FREQ, backing="pec")
    jx, jy = aperture_current_from_scalar(design.x_m, design.y_m, design.aperture_field, feed_xy_m=(float(x0), float(y0)))
    pat = compute_vector_pattern(design.x_m, design.y_m, jx, jy, FREQ, green, PERIOD, PERIOD)
    peak = float(np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2))
    vortex = measure_vortex_winding(
        design.x_m, design.y_m, jx, jy, FREQ, green, PERIOD, PERIOD, 20.0, 0.0, "rhcp", pattern_peak_power=peak,
    )
    assert vortex.winding_number == pytest.approx(-1.0, abs=0.05)
    assert vortex.null_depth_db <= -25.0
