"""Side quest: switched parasitic-pin launcher network model (scripts/research/sidequest_launcher_model.py) — conservation / symmetry / limit checks."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
warnings.simplefilter("ignore")
import sidequest_launcher_model as L  # noqa: E402


def _on(n, back=True):
    return [(np.cos(np.radians((m + 0.5) * 360 / n)) < 0) == back for m in range(n)]


def test_bare_probe_is_the_calibrated_feed():
    r = L.solve(4, 10e-3, 6.0, [False] * 4, z_off=1e12)
    assert r["s11_db"] == pytest.approx(-15.0, abs=0.2)
    assert r["p_sw"] / r["p_acc"] == pytest.approx(L.pin(6.0)["launch"], abs=2e-3)  # = the PinLaunch launch fraction
    fb, _ = L.fb_2q(r)
    assert abs(fb) < 1e-6  # symmetric: no steering


@pytest.mark.parametrize("n,rho,top", [(12, 10.5e-3, 6.0), (8, 6e-3, 8.0), (4, 15e-3, 5.0)])
def test_power_balance_closes(n, rho, top):
    r = L.solve(n, rho, top, _on(n))
    assert (r["p_sw"] + r["p_space"] + r["p_diode"]) / r["p_acc"] == pytest.approx(1.0, abs=2e-3)  # surface wave + space radiation + diode loss = accepted


def test_mirror_symmetry_of_the_two_states():
    ra, rb = L.solve(12, 10.5e-3, 6.0, _on(12, True)), L.solve(12, 10.5e-3, 6.0, _on(12, False))
    fa, fb_ = L.fb_2q(ra)[0], L.fb_2q(rb)[0]
    assert fa == pytest.approx(-fb_, abs=1e-6) and ra["s11_db"] == pytest.approx(rb["s11_db"], abs=1e-6)


def test_nominal_design_meets_the_targets_and_a_post_inside_the_substrate_cannot_scatter():
    r = L.solve(12, 10.5e-3, 6.0, _on(12))
    assert L.fb_2q(r)[0] >= 10.0 and r["s11_db"] <= -10.0
    assert L.pin(0.0)["R_sw"] < 0.5  # a shorted post that stays inside the 1.575 mm substrate radiates ~0.26 ohm into the surface wave: a negligible scatterer
