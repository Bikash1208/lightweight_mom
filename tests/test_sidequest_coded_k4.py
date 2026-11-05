"""K-state engine (K = 2, 4): word bookkeeping and the K = 2 consistency with the earlier bit design."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_coded_k4 as K  # noqa: E402


def test_words_and_options():
    assert len(K.words(2)) == 2 and len(K.words(4)) == 14 and len(K.words(4, allw=True)) == 16
    opt = K.options(4, np.array([-1.0, 1.0]), [(1.0, -1.0, 1.0)])
    assert opt[0].shape == (4 + 14 * 2 * 2 + 14 * 14, 4)
    mx, my, nd = opt[0], opt[1], opt[2]
    assert set(np.unique(nd)) == {0.0, 1.0, 2.0} and np.all(np.abs(mx) <= 1) and np.all(np.abs(my) <= 1)


def test_line_classes_count_complements_once():
    bd = K.KBoard(*K.SETS["K4b"])
    W = K.words(4)
    ids = (np.array([0, 13, 1, -1]), np.array([-1, -1, -1, -1]))  # word 0 and 13 are complements (0001/1110)
    res = K.lines(bd, ids, None, W)
    assert res["x"] == {"words": 3, "classes": 2} and res["y"] == {"words": 0, "classes": 0}


def test_k2_bit_design_matches_the_two_state_engine():
    bd = K.KBoard(*K.SETS["K2"])
    r, _ = K.run(bd, K.options(2, np.array([-1.0, 1.0]), [(1.0, -1.0, 1.0)]), [0, 0], lam=0.555, iters=2)
    assert 700 <= r["diodes"] <= 800 and abs(r["states"][0]["gain"] - 20.56) < 0.15
    assert r["lines"]["x"]["classes"] == 1 and r["lines"]["x"]["words"] == 2
