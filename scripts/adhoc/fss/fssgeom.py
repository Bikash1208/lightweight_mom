import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))
from dataclasses import dataclass

from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice


@dataclass(frozen=True)
class SquareLoop:
    outer: float
    inner: float
    center_x: float = 0.0
    center_y: float = 0.0
    rotation: float = 0.0
    kind = "square_loop"

    @property
    def polarity(self):
        return "metal"

    @property
    def feature_area(self):
        return self.outer**2 - self.inner**2

    def corners(self):
        h = self.outer / 2
        return np.array([[h, h], [-h, h], [-h, -h], [h, -h]])

    def in_feature(self, x, y):
        ax = np.abs(np.asarray(x, float) - self.center_x)
        ay = np.abs(np.asarray(y, float) - self.center_y)
        return (np.maximum(ax, ay) <= self.outer / 2) & (np.maximum(ax, ay) >= self.inner / 2)


P = 10e-3
H = 1.0e-3


def unit_cell():
    return UnitCell(
        RectangularLattice(P, P),
        SquareLoop(8.8e-3, 6.6e-3),
        Substrate("RT5880", 2.2, 0.0009, H),
        Metal(),
        backing="air",
    )
