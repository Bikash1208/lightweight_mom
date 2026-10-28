"""Phase 29 B3 diagnostic at R 300 mm (MODEL_ONLY): helper checks (gene decoding, sidelobe diagnostics on synthetic patterns, the user's DoF ceiling formula)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/research"))
import phase29_B3_scale300 as s3  # noqa: E402


def test_decode_four_feed_ring_radius_from_adjacent_spacing():
    g = np.array([40.0, 0.0, 0, 0, 0, 0, 3.0])
    ports, off, dfo, q0 = s3.decode(g, 4)
    r = [np.hypot(p[0], p[1]) for p in ports]
    assert r == pytest.approx([40.0 / (2 * np.sin(np.pi / 4))] * 4)
    d = np.hypot(ports[0][0] - ports[1][0], ports[0][1] - ports[1][1])
    assert d == pytest.approx(40.0)  # adjacent feeds are one spacing apart
    assert [p[3] for p in ports] == [0, 90, 180, 270] and q0 == pytest.approx(1000.0) and not dfo.any()


def test_decode_seven_and_one_feed():
    p7, off, dfo, q0 = s3.decode(np.concatenate([[40.0, 0.0], np.zeros(6), np.arange(7.0), np.full(7, 2.0), [2.0]]), 7)
    assert (
        len(p7) == 7
        and p7[0][:2] == (0.0, 0.0)
        and all(p[2] == s3.TILT for p in p7[1:])
        and list(off) == list(range(7))
        and q0 == pytest.approx(100.0)
    )
    p1, *_ = s3.decode(np.array([3.0]), 1)
    assert p1 == [(0.0, 0.0, 0.0, 0.0)]


def test_dof_ceiling_formula_of_the_user():
    om = 2 * np.pi * (1 - np.cos(np.radians(11.0)))
    assert 10 * np.log10(4 * np.pi * 7 / om) == pytest.approx(28.8, abs=0.05)  # +-11 deg, 7 beams


def test_sidelobe_diagnostics_on_a_known_pattern():
    th, ph = np.arange(0.0, 90.1, 1.0), np.arange(0.0, 360.0, 1.0)
    tg, pg = np.meshgrid(th, ph, indexing="ij")
    a = np.exp(-2 * np.log(2) * (tg / 8.0) ** 2) + 0.1 * np.exp(
        -2 * np.log(2) * ((tg - 40.0) / 8.0) ** 2
    )  # beam at 0 deg + a -20 dB lobe at 40 deg
    f = (a.astype(complex), -1j * a)  # pure RHCP (project convention)
    pm = {"peak_theta_phi": (0.0, 0.0), "hpbw_deg": 8.0}
    d = s3.sidelobe_diag(f, th, ph, pm)
    assert d["worst_rhcp_sidelobe_db"] == pytest.approx(-20.0, abs=0.7) and d[
        "worst_sidelobe_deg_from_peak"
    ] == pytest.approx(40.0, abs=2.0)
    assert d["lhcp_power_fraction"] < 1e-6 and 0.7 < d["power_fraction_in_1p5hpbw"] < 0.95
