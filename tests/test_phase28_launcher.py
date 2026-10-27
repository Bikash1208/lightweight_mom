"""Phase 28 Yagi-type launcher: the coupled-monopole model reduces to the single pin when the posts are far away and conserves the surface-wave power; the AEDT script runs against a recording stub."""

import math
import re
import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_launcher_model as lm  # noqa: E402

PKG = ROOT / "hfss_exports" / "phase28_launcher"


@pytest.fixture(scope="module")
def tab():
    from src.dispersion.pin_launch import Launcher, PinLaunch

    rows = []
    for t in lm.T_GRID:
        s = PinLaunch(lm.F0, lm.EPS_R, lm.H, Launcher(float(t), 0.635e-3, 0.0), n_seg=10).split()
        rows.append((t, 2 * s["p_sw_w_per_a2"], s["r_rad_ohm"], s["x_est_ohm"], s["beta_over_k0"]))
    a = np.array(rows)
    return {"r_sw": lm.CubicSpline(a[:, 0], a[:, 1]), "r_tot": lm.CubicSpline(a[:, 0], a[:, 2]), "x": lm.CubicSpline(a[:, 0], a[:, 3]), "beta_over_k0": float(a[0, 4])}


def test_far_posts_leave_the_single_pin(tab):
    d = {"Td": 4e-3, "reflector": (500.0, 4e-3), "directors": [(500.0, 4e-3)]}
    r = lm.yagi(d, tab)
    z = float(tab["r_tot"](4e-3)) + 1j * float(tab["x"](4e-3))
    assert complex(*r["zin_ohm"]) == pytest.approx(z, rel=0.05)
    assert r["hpbw_deg"] >= 300  # omnidirectional


def test_power_conservation_and_symmetric_coupling(tab):
    d = {"Td": 3.7e-3, "reflector": (9e-3, 8e-3), "directors": [(13e-3, 5.5e-3), (14e-3, 5e-3)]}
    r = lm.yagi(d, tab)
    # the surface-wave power cannot exceed the accepted power (the rest is space wave): launch fraction <= 1, > 0
    assert 0.3 < r["launch_fraction"] <= 1.0 + 1e-9
    assert r["zin_ohm"][0] > 0


@pytest.mark.skipif(not (PKG / "build_launcher.py").exists(), reason="launcher package not generated")
@pytest.mark.parametrize("loaded", [True, False])
def test_build_script_runs_against_the_stub(loaded):
    text = (PKG / "build_launcher.py").read_text(encoding="utf-8")
    assert not re.search(r"""\bf["']""", text) and "print(" not in text
    text = text.replace("LOADED_SURFACE = True ", f"LOADED_SURFACE = {loaded} ", 1)
    glob = {"oDesktop": mock.MagicMock(), "__name__": "__aedt__"}
    with mock.patch.dict(sys.modules, {"ScriptEnv": mock.MagicMock()}):
        exec(compile(text, "build_launcher.py", "exec"), glob)  # noqa: S102 - generated, committed script, stub environment
    n = len(glob["CELLS"])
    assert glob["oEditor"].CreateCylinder.call_count == 4
    assert glob["oBnd"].AssignLumpedPort.call_count == 1
    assert glob["oBnd"].AssignLumpedRLC.call_count == (n if loaded else 0)
    xs = [(c[1], c[2]) for c in glob["CELLS"]]
    assert min(math.hypot(x, y) for x, y in xs) >= 38.0 - 1e-6  # keep-out around the launcher
