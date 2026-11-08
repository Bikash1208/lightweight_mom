import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sub6_slotdisc2_lib as T  # noqa: E402
import sub6_slotdisc_lib as B  # noqa: E402


def test_zero_second_gap_reproduces_single_notch():
    rng = np.random.default_rng(0)
    x, y = rng.uniform(-6e-3, 6e-3, 4000), rng.uniform(-6e-3, 6e-3, 4000)
    a = B.SlotDiscNotch(4.5e-3, 3.6e-3, 2.4e-3, 1e-3, bridge=True)
    b = T.SlotDiscNotch2(4.5e-3, 3.6e-3, 2.4e-3, 1e-3, bridge=True, notch2_width=0.0)
    assert np.array_equal(a.in_feature(x, y), b.in_feature(x, y))


def test_second_gap_splits_disc_and_is_one_cell_on_the_chosen_mesh():
    for g2, n in T.MESH_FOR_G2.items():
        k, w = T.gap_cells(11.6, 4.9, 7.8, 2.9, g2, n)
        assert k == {0.36: 3}.get(g2, 1) and abs(w - g2) / g2 < 0.25, (g2, n, k, w)
