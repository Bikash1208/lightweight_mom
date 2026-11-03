"""Phase 30 library (scripts/research/phase30_library.py): interpolation properties and the stored held-out validation."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
LIB = ROOT / "results" / "phase30" / "library"

import phase30_library as L  # noqa: E402


def test_trig_interp_reproduces_samples_and_is_180_periodic():
    rng = np.random.default_rng(0)
    v = rng.normal(size=(12, 3))
    for i in (0, 3, 7):
        assert np.allclose(L.trig_interp(v, L.ANG[i], 0), v[i], atol=1e-12)
    a = 37.3
    assert np.allclose(L.trig_interp(v, a, 0), L.trig_interp(v, a + 180.0, 0), atol=1e-12)


def test_slot_sizes_are_even_mesh_cells_and_leave_a_bridge():
    for nl, nw in L.slots():
        assert nl % 2 == 0 and nw % 2 == 0
        if nl:
            assert (L.MAJOR - nl * L.HC) / 2 >= 0.3e-3
    assert all(n % 2 == 0 for n in L.HELD_NL)
    assert not set(L.HELD_NL) & set(L.NL)


@pytest.mark.skipif(not (LIB / "angle_lut_v1.npz").exists(), reason="library not generated")
def test_lut_is_exact_on_grid_c4_symmetric_and_reciprocal():
    d = np.load(LIB / "angle_lut_v1.npz")
    Y, sl = d["Y"], [tuple(map(int, s)) for s in d["slots_cells"]]
    i = sl.index((22, 6))
    assert np.allclose(L.lut_eval(Y, sl, 22, 6, L.ANG[3], L.ANG[5]), Y[i, 3, 5], rtol=1e-9, atol=1e-12)
    meta = json.loads((LIB / "angle_lut_v1.json").read_text())
    assert meta["c4_relative_error_max"] < 1e-9  # square lattice + fixed outline: exact C4 of the mask
    assert np.allclose(Y[..., 0, 1], Y[..., 1, 0], rtol=1e-6, atol=1e-9)  # reciprocity


@pytest.mark.skipif(not (LIB / "heldout.json").exists(), reason="held-out run not generated")
def test_heldout_validation_recorded_and_below_five_percent():
    h = json.loads((LIB / "heldout.json").read_text())
    assert h["n"] == 60 and len(h["rows"]) == 60
    assert h["rel_err_max"] < 0.05 and h["rel_err_median"] < 0.01
    assert all(r["nl"] in L.HELD_NL for r in h["rows"])  # none of the held-out lengths is a training length


@pytest.mark.skipif(not (LIB / "uniform_table.json").exists(), reason="library not generated")
def test_uniform_table_monotone_in_slot_length_and_no_slot_isotropic():
    rows = {(r["nl"], r["nw"]): r for r in json.loads((LIB / "uniform_table.json").read_text())["rows"]}
    assert rows[(0, 0)]["M_t"] < 1e-3
    for w in L.NW:
        m = [rows[(n, w)]["M_t"] for n in L.NL[1:]]
        assert all(b > a for a, b in zip(m, m[1:], strict=False))
