"""Phase 28 (loop cell R4): 1-bit sampling limit at 9.0 mm pitch. For beams toward the feed (phi = 180) and away from it (phi = 0), and sideways (phi = 90), the fraction of cells whose modulation phase advance per lattice step along x or y exceeds pi,
using the cell's own beta(psi) (unrotated tensor, MADP L_s = 0). python scripts/research/phase28_loop_alias_table.py <out.json>"""
import json, math, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts" / "research"))
import phase28_loop_sector_model as lm
import phase28_sector_feed_model as sf
from src.hologram import coded_feasibility as cf
cf.SUB = lm.SUB; cf.set_frequency(10e9)
yon, yoff = lm.tensors("MEAS_10mA"); ybar = 0.5 * (yon + yoff)
poly, src = sf.source_xy()
x, y = cf.lattice(poly, lm.PITCH, keepout_m=0.0)
rho = np.hypot(x - src[0], y - src[1]); k = rho >= sf.RHO_MIN
x, y = x[k], y[k]; psi = np.arctan2(y - src[1], x - src[0])
grid = np.radians(np.linspace(-90, 90, 91))
bg = np.array([lm.beta_of_b((np.array([math.cos(p), math.sin(p)]) @ ybar @ np.array([math.cos(p), math.sin(p)])).imag * 1e3) for p in grid])
beta = np.interp(psi, grid, bg)
out = {"pitch_mm": 9.0, "beta_over_k0_at_psi0": float(bg[45] / cf.K0), "rows": []}
for phi in (180.0, 0.0, 90.0):
    for th in (15, 20, 25, 30, 35, 40, 45, 50, 60):
        kt = cf.K0 * math.sin(math.radians(th)); a = math.radians(phi)
        kx, ky = beta * np.cos(psi) - kt * math.cos(a), beta * np.sin(psi) - kt * math.sin(a)
        frac = float(np.mean((np.abs(kx) * lm.PITCH > math.pi) | (np.abs(ky) * lm.PITCH > math.pi)))
        out["rows"].append({"theta_deg": th, "phi_deg": phi, "aliased_cell_fraction": frac})
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
for r in out["rows"]: print(r["phi_deg"], r["theta_deg"], round(r["aliased_cell_fraction"], 2))
