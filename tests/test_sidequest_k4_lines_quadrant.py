"""Class-restricted options (line minimisation) and the quadrant (segmented 4-probe) board geometry."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sidequest_coded_k4 as K  # noqa: E402
import sidequest_quadrant_board as Q  # noqa: E402


def test_class_restriction_keeps_only_allowed_complement_classes():
    cls = {(0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0)}
    opt = K.options(4, np.array([-1.0, 1.0]), [(1.0, -1.0, 1.0)], cls_x=cls, cls_y=cls)
    W = opt[6]
    used = {K.canon(W[i]) for i in set(opt[4].tolist()) | set(opt[5].tolist()) if i >= 0}
    assert used == cls
    full = K.options(4, np.array([-1.0, 1.0]), [(1.0, -1.0, 1.0)])
    assert opt[0].shape[0] < full[0].shape[0]


def test_wedge_is_a_quarter_of_the_hexagon_and_cells_stay_inside():
    hexa_area = 58.457e-3
    for phi in (0.0, np.pi / 2):
        poly = Q.wedge_poly(phi)
        _, area = Q.centroid(poly)
        assert abs(area - hexa_area / 4) / (hexa_area / 4) < 0.06
    bd = Q.QBoard(0, "RHCP", feed_frac=0.75)
    ang = np.arctan2(bd.y, bd.x)
    assert np.all(np.abs(ang) <= np.radians(45) + 1e-6) and bd.feed[0] > 0.04 and abs(bd.feed[1]) < 1e-9
    assert np.all(bd.rho >= 6.5e-3 - 1e-12)
