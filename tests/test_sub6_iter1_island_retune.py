"""Geometry regression and passive-circuit checks; no expensive solve."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1]/"scripts/research/sub6_iter1_island_retune.py"
spec = importlib.util.spec_from_file_location("retune",PATH)
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def test_baseline_mask_matches_established_cell():
    import diag_full_tensor_lib as old
    old_geo = old.SlotDiscIsland(.005435,36,12,18,5,11.6/96)
    new_geo = R.Island(R.BASE,96)
    for offset in (0,.5):
        grid = (np.arange(96)-48+offset)*.0116/96
        x,y = np.meshgrid(grid,grid)
        np.testing.assert_array_equal(new_geo.in_feature(x,y),old_geo.in_feature(x,y))


@pytest.mark.parametrize("n,gap",[(96,8*11.6/96),(128,11*11.6/128)])
def test_one_bridge_and_clearances(n,gap):
    g = R.Island(R.BASE,n)
    assert g.gap == pytest.approx(gap)
    assert g.in_feature(g.x*.001,g.y*.001)
    assert not g.in_feature((g.x+.3)*.001,g.y*.001)
    assert not g.in_feature(0,(-g.sw+g.cb/2)*.001)
    assert g.in_feature(0,(g.bot+g.top)/2*.001)


def test_passive_circuit_and_nominal_impedance():
    f = 5.8e9
    on = R.impedance(f,.11,4000,"ON")
    off = R.impedance(f,.11,4000,"OFF")
    assert on == pytest.approx(3+1j*2*np.pi*f*.6e-9)
    assert off.real > 0 and off.imag < 0
    assert off.real*.35 == pytest.approx(5.42,abs=.01)
    assert off.imag*.35 == pytest.approx(-79.320,abs=.002)
    for _,c,rp in R.CASES:
        assert R.impedance(f,c,rp,"OFF").real > 0
