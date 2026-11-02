"""Phase 30 item 2, stage B (MODEL_ONLY): for a short list of cells from `cell_screen.json`
  (a) arrival-angle sensitivity: Bloch sheet tensor at |kt| = beta0, direction phi = 0..90 deg to the cell's x axis,
      against the local-homogeneous prediction Y_kk(phi) = Y_par cos^2 + Y_perp sin^2 (+ 2 Y_xy cos sin);
      deviation relative to the anisotropy (0 = local tensor, the Phase 26/29 model assumption);
  (b) Floquet harmonic-balance cross-check of the first-order alpha: rotating-axis staircase of N cells,
      `tensor_floquet.floquet_complex_beta_tensor`, same Y_bar and delta, period = N pitch ~ lambda_sw.

    python scripts/research/phase30_cell_stageb.py <screen_json> <out_json> [--workers 20]
Picks are listed in PICKS below (filled in from the ranked table).
"""

import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase30_cell_screen as S  # noqa: E402
from src.basis_functions import MeshSettings, RooftopBasis  # noqa: E402
from src.dispersion import tensor_sheet as ts  # noqa: E402
from src.dispersion.elliptical_cells import EllipticalPatch  # noqa: E402
from src.dispersion.sheet_model import find_beta_from_sheet_admittance  # noqa: E402
from src.dispersion.tensor_floquet import floquet_complex_beta_tensor  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402

PHIS = (0.0, 22.5, 45.0, 67.5, 90.0)


def stage_b(row, mesh_n):
    eps, tand, h, _ = S.SUBSTRATES[row["substrate"]]
    sub = Substrate(row["substrate"].lower(), eps, tand, h)
    pitch = row["pitch_mm"] * 1e-3
    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    geo = EllipticalPatch(row["major_mm"] * 1e-3, row["major_mm"] * 1e-3 * row["ratio"], 0.0,
                          row["slot_L_mm"] * 1e-3, row["slot_W_mm"] * 1e-3, 0.0, "cross" in row["family"], row["exponent"])
    cell = UnitCell(RectangularLattice(pitch, pitch), geo, sub, Metal(), backing="pec")
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    beta0 = row["beta_over_k0"] * S.K0
    out = dict(row)
    yk = {}
    for phi in PHIS:
        p = np.radians(phi)
        y = ts.bloch_sheet_tensor(cell, S.F, mesh, beta0 * np.cos(p), beta0 * np.sin(p), basis=basis)
        kh = np.array([np.cos(p), np.sin(p)])
        yk[phi] = float((kh @ y @ kh).imag)
    y_par, y_perp = yk[0.0], yk[90.0]
    delta = 0.5 * abs(y_par - y_perp)
    dev = {phi: yk[phi] - (y_par * np.cos(np.radians(phi)) ** 2 + y_perp * np.sin(np.radians(phi)) ** 2) for phi in PHIS}
    out["Ykk_mS"] = {str(k): v * 1e3 for k, v in yk.items()}
    out["local_model_dev_over_delta"] = max(abs(v) for v in dev.values()) / max(delta, 1e-12)
    # beta vs arrival angle from Y_kk(phi) alone (scalar resonance along k): spread of the bound constant
    betas = []
    for phi in PHIS:
        try:
            betas.append(find_beta_from_sheet_admittance(1j * yk[phi], sub, S.F).beta_rad_per_m / S.K0)
        except Exception:
            betas.append(float("nan"))
    out["beta_over_k0_vs_phi"] = betas
    out["beta_spread_over_k0"] = float(np.nanmax(betas) - np.nanmin(betas))
    # Floquet cross-check, rotating-axis staircase
    try:
        y_bar = 0.5 * (y_par + y_perp)
        n_cell = max(int(round(2.0 * np.pi / beta0 / pitch)), 4)
        tens = []
        for i in range(n_cell):
            a = np.pi * i / n_cell  # axis angle = psi/2, psi = 2 pi i / N
            c, s = np.cos(a), np.sin(a)
            r = np.array([[c, -s], [s, c]])
            tens.append(1j * (y_bar * np.eye(2) + delta * r @ np.diag([1.0, -1.0]) @ r.T))
        beta_f, _ = floquet_complex_beta_tensor(np.array(tens), n_cell * pitch, sub, S.F, beta0, n_max=4)
        out["N_cells_period"] = n_cell
        out["alpha_floquet_np_m"] = float(abs(beta_f.imag))
        out["alpha_floquet_over_k0"] = float(abs(beta_f.imag)) / S.K0
        out["alpha_over_first_order"] = out["alpha_floquet_np_m"] / row["alpha_np_m"] if row["alpha_np_m"] else None
    except Exception as e:  # noqa: BLE001
        out["floquet_error"] = f"{type(e).__name__}: {e}"
    return out


def _job(a):
    return stage_b(*a)


def main(screen_json, out_json, mesh_n=64, workers=20):
    rows = json.loads(Path(screen_json).read_text())["rows"]
    picks = []
    for want in PICKS:
        for r in rows:
            if r.get("ok") and all(abs(r[k] - v) < 1e-6 if isinstance(v, float) else r[k] == v for k, v in want.items()):
                picks.append(r)
                break
    t0 = time.time()
    with Pool(min(workers, len(picks))) as pool:
        res = pool.map(_job, [(r, mesh_n) for r in picks])
    Path(out_json).write_text(json.dumps({"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "mesh": mesh_n, "rows": res, "seconds": time.time() - t0}, indent=1))
    print("done", len(res), f"{time.time() - t0:.0f}s")


PICKS = []  # list of dicts matching screen rows, set by the driver after ranking

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("screen_json")
    ap.add_argument("out_json")
    ap.add_argument("--picks", help="JSON file: list of dicts selecting screen rows")
    ap.add_argument("--workers", type=int, default=20)
    a = ap.parse_args()
    PICKS = json.loads(Path(a.picks).read_text())
    main(a.screen_json, a.out_json, workers=a.workers)
