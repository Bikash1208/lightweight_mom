"""Phase 29 design point B3 (MODEL_ONLY): angle-dependent superposition (Frey et al. 2023 in the cell model) and the B3 evaluation helpers: known-answer checks."""

import sys
from pathlib import Path

import numpy as np
import pytest

from src.hologram import seven_hex_cluster as c
from src.hologram import shared_aperture_frey as fr
from src.hologram import shared_aperture_multiport as sa
from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/research"))


def test_feed_ports_centre_plus_ring_opposite_the_beam_with_jitters():
    p = fr.feed_ports(40.0, rot_deg=10.0, radial_jitter_mm=[1, 0, 0, 0, 0, -2], tilt_deg=13.0)
    assert p[0] == (0.0, 0.0, 0.0, 0.0) and len(p) == 7
    assert all(q[2] == 13.0 for q in p[1:]) and [q[3] for q in p[1:]] == [0, 60, 120, 180, 240, 300]
    r = [np.hypot(q[0], q[1]) for q in p[1:]]
    assert r == pytest.approx([41.0, 40.0, 40.0, 40.0, 40.0, 38.0])
    az = np.degrees(np.arctan2(p[1][1], p[1][0]))
    assert (az - (0.0 + 180.0 + 10.0) + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=1e-9)


@pytest.fixture(scope="module")
def table():
    return SupercellTable.from_files(t.TABLE_GLOB)


def test_single_feed_frey_and_aia_coincide_and_orientation_follows_the_required_angle(table):
    ports = [(0.0, 0.0, 0.0, 0.0)]
    a = fr.design_frey(table, 40.0, ports, [0.0], 1000.0, 1e5, mode="frey", n_iter=0, n_psi=90)
    b = fr.design_frey(table, 40.0, ports, [0.0], 1000.0, 1e5, mode="aia", n_iter=0, n_psi=90)
    assert np.array_equal(a["slot"], b["slot"]) and np.allclose(a["axis_rad"], b["axis_rad"])
    assert (
        a["mean_abs_orientation_error_rad_history"][0] < 1e-6
    )  # one feed: the shared slot realises its own sub-hologram exactly (orientation)


def test_offsets_rotate_the_slot_axes_by_half_the_offset(table):
    ports = fr.feed_ports(20.0, tilt_deg=13.0)[:1]
    a = fr.design_frey(table, 40.0, ports, [0.0], 1000.0, 1e5, n_iter=0, n_psi=90)
    b = fr.design_frey(table, 40.0, ports, [90.0], 1000.0, 1e5, n_iter=0, n_psi=90)
    d = np.mod(b["axis_rad"] - a["axis_rad"], np.pi)
    assert np.allclose(d, np.pi / 4)  # offset psi0 = 90 deg -> axis shift psi0 / 2 = 45 deg


def test_more_feeds_do_not_reduce_the_orientation_error(table):
    e = []
    for k in (1, 3, 7):
        ports = fr.feed_ports(40.0, tilt_deg=13.0)[:k]
        d = fr.design_frey(table, 60.0, ports, [0.0] * k, 1000.0, 1e5, n_iter=0, n_psi=90)
        e.append(d["mean_abs_orientation_error_rad_history"][0])
    assert e[0] < 1e-6 and e[1] > e[0] and e[2] >= e[1] - 0.05


def _beam(th, ph, th0, ph0, w):
    tg, pg = np.meshgrid(th, ph, indexing="ij")
    ang = np.degrees(np.arccos(np.clip(np.sum(c._unit(tg, pg) * c._unit(th0, ph0), axis=-1), -1, 1)))
    a = np.exp(-2 * np.log(2) * (ang / w) ** 2)  # amplitude; power -3 dB at ang = w / 2
    return a, 1j * a  # (E_theta, E_phi) = (a, j a): pure RHCP


def test_crossover_of_two_gaussian_beams_one_hpbw_apart_is_minus_3_db():
    import phase29_pointB3 as b3

    th, ph = np.arange(0.0, 90.1, 1.0), np.arange(0.0, 360.0, 1.0)
    F = [_beam(th, ph, 0.0, 0.0, 13.0)] + [_beam(th, ph, 13.0, 60.0 * k, 13.0) for k in range(6)]
    pm = [{"peak_theta_phi": (0.0, 0.0)}] + [{"peak_theta_phi": (13.0, 60.0 * k)} for k in range(6)]
    x = b3.crossovers(F, th, ph, pm)
    assert x["centre_outer"]["mean_db"] == pytest.approx(-3.0, abs=0.4)
    assert x["outer_outer"]["mean_db"] == pytest.approx(-3.0, abs=0.8)


def test_sum_and_difference_patterns_have_a_null_at_the_sum_peak_for_symmetric_beams():
    import phase29_pointB3 as b3

    th, ph = np.arange(0.0, 90.1, 1.0), np.arange(0.0, 360.0, 1.0)
    F = [_beam(th, ph, 0.0, 0.0, 20.0)] + [_beam(th, ph, 13.0, 60.0 * k, 20.0) for k in range(6)]
    r = b3.sum_diff(F, th, ph)
    assert r["sum"]["peak_theta_phi_deg"][0] < 2.0  # symmetric in-phase sum peaks at broadside
    for k in ("delta_x", "delta_y"):
        assert r[k]["null_depth_at_sum_peak_db"] < -25.0
