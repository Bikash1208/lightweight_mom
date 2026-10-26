"""Phase 26 B2.4f: what conversion can the SOLVER_SUPERCELL table reach over the 150 mm hexagon?

    python scripts/research/phase26_table_reach.py <out.json>

Best case: every ray uses, at its local modulation wavenumber K = |beta rho_hat - k0 sin(theta) b_hat|
(beta = 1.065 k0, theta 30 deg, azimuth 225 deg), the LARGEST alpha_leak of the table (any depth), linearly
interpolated in K; conversion of a ray `1 - exp(-2 alpha rho_edge)` (alpha constant along the ray: the
launched wave's power per radian decays as exp(-2 alpha rho); cylindrical spreading does not change the
FRACTION). Collinear (gamma = 0) table applied to every ray: an approximation, labelled.
"""

import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.hologram.leaky_forward import polygon_ray_distance  # noqa: E402

pts = []
for f in sorted(glob.glob(str(ROOT / "results/phase26/supercell_table_*.json"))):
    pts += json.load(open(f))["points"]
ks = sorted({round(p["K_over_k0"], 4) for p in pts})
amax = [max(p["alpha_leak_np_per_m"] for p in pts if round(p["K_over_k0"], 4) == k) for k in ks]
R = 0.150
v = np.array([[R * np.cos(np.radians(60 * i)), R * np.sin(np.radians(60 * i))] for i in range(6)])
psi = np.radians(np.arange(0, 360, 1.0))
red = polygon_ray_distance(v, (0.0, 0.0), psi)
b, st = 1.065, 0.5
kk = np.sqrt(b * b + st * st - 2 * b * st * np.cos(psi - np.radians(225)))
a = np.interp(kk, ks, amax)
eta = 1 - np.exp(-2 * a * red)
need = -np.log(0.15) / (2 * red)  # alpha for 85 % on each ray if constant
out = {
    "labels": ["SOLVER_SUPERCELL", "MODEL_ONLY", "GAMMA_0_COLLINEAR_ONLY"],
    "K_over_k0_grid": ks, "alpha_max_np_per_m": amax,
    "K_range_over_k0": [float(kk.min()), float(kk.max())],
    "best_case_conversion_mean": float(eta.mean()), "best_case_conversion_min": float(eta.min()),
    "best_case_conversion_max": float(eta.max()),
    "alpha_needed_for_85pct_mean_np_per_m": float(need.mean()),
    "alpha_available_mean_np_per_m": float(a.mean()),
}
print(json.dumps(out, indent=1))
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
