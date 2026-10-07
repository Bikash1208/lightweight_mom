import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from dataclasses import dataclass

from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.sparameters import compute_s_parameters


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
uc = UnitCell(
    RectangularLattice(P, P),
    SquareLoop(8.8e-3, 6.6e-3),
    Substrate("RT5880", 2.2, 0.0009, 1.0e-3),
    Metal(),
    backing="air",
)
planes = (
    ReferencePlanes(top_air_thickness=0.0, bottom_air_thickness=0.0, substrate_thickness=1.0e-3)
    if len(sys.argv) < 4
    else None
)
n = int(sys.argv[1])
fl = [float(f) * 1e9 for f in sys.argv[2].split(",")]
mesh = MeshSettings(cells_x=n, cells_y=n, floquet_order=int(sys.argv[3]) if len(sys.argv) > 3 else n)
for f in fl:
    t = time.time()
    s = compute_s_parameters(uc, f, mesh, planes, warn_grating=False)
    s11 = s.s11[0, 0]
    s21 = s.s21[0, 0]
    print(
        f"n={n} f={f / 1e9:.3f} S11 {20 * np.log10(abs(s11)):.3f} dB {np.degrees(np.angle(s11)):.2f} deg | S21 {20 * np.log10(abs(s21)):.3f} dB {np.degrees(np.angle(s21)):.2f} deg | yy S21 {20 * np.log10(abs(s.s21[1, 1])):.3f} | xpol {20 * np.log10(abs(s.s21[1, 0]) + 1e-15):.1f} | {time.time() - t:.1f}s",
        flush=True,
    )
