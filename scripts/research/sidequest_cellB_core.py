"""Architecture (b) back-switched mushroom, item 2 core.  MODEL_ONLY.
Patch grid: MoM Bloch sheet tensor (grounded slab, no via) at the actual transverse wavenumber (beta cos psi, beta sin psi).
Via + diode branch: lumped Z_cell = j w L_b + R / (1/(j w C)) per square cell (square cell: Z_sheet = Z_cell), added to the sheet admittance (Sievenpiper/Luukkonen mushroom model);
the via is NOT in the MoM (zero-thickness sheet).  Transverse resonance: Im[Y_air(beta) + Y_slab(beta) + Y_patch(beta, psi) + Y_branch] = 0.
    python sidequest_cellB_core.py test"""
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

from src.composite_geometry import AxisRect, CompositeGeometry
from src.geometry import UnitCell
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice

MM = 1e-3
C0 = 299792458.0
EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi

LAMINATES = {
    # name: (er, tand_10GHz) -- data sheet values, see LAMINATES.md for sources
    "RO3006": (6.15, 0.0020),
    "RT5880": (2.2, 0.0009),
    "RT5880_lo": (2.18, 0.0009),  # assumed Dk tolerance +-0.02 (Rogers RT/duroid 5880 datasheet value, not re-verified here)
    "RT5880_hi": (2.22, 0.0009),
    "RO3003": (3.00, 0.0010),
    "RO4003C": (3.38, 0.0027),
    "RO4003C_lo": (3.33, 0.0027),  # Dk -0.05 (Rogers spec +-0.05, process)
    "RO4003C_hi": (3.43, 0.0027),
}


def patch_cell(p, gap, er, tand, h, ssq=None):
    """Square patch (side p - gap unless ssq given) on a grounded slab, period p [mm]."""
    s = (p - gap) if ssq is None else ssq
    geo = CompositeGeometry(metal=(AxisRect(0.0, 0.0, s * MM, s * MM),))
    return UnitCell(RectangularLattice(p * MM, p * MM), geo, Substrate("lam", er, tand, h * MM), Metal(), backing="pec")


def patch_Y(args):
    """args = (p, gap, er, tand, h, f, mesh, beta_rad_m, psi_deg) -> 2x2 complex sheet tensor [S] of the patch grid at (beta cos psi, beta sin psi)."""
    import sidequest_cell_search as C

    p, gap, er, tand, h, f, mesh, beta, psi = args
    cell = patch_cell(p, gap, er, tand, h)
    y = C.bloch_loaded(cell, f, mesh, beta * np.cos(np.radians(psi)), beta * np.sin(np.radians(psi)))
    return y


if __name__ == "__main__":
    import time

    f = 10.5e9
    t0 = time.time()
    k0 = 2 * np.pi * f / C0
    for mesh in (48, 64):
        t0 = time.time()
        y = patch_Y((6.4, 0.4, 6.15, 0.002, 1.28, f, mesh, 1.25 * k0, 0.0))
        print(mesh, np.round(y * 1e3, 3).tolist(), "%.1f s" % (time.time() - t0), flush=True)
