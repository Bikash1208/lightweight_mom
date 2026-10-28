"""Phase 29: angle-dependent tensor cell LUT (SOLVER_SUPERCELL-family, MODEL_ONLY): the full 2 x 2 Bloch sheet admittance Y_s of the B3 cell family (slotted circle D 6.0 mm on the 6.4 mm lattice,
RT/duroid 5880 1.575 mm, grounded, 10.54 GHz) vs slot width (the 4 table slots), slot orientation a (0-165 deg, 15 deg steps, relative to the LATTICE) and direction psi of the incident surface wave's
transverse wavenumber kt = 1.065 k0 (0-165 deg, 15 deg steps), from periodic MoM (`tensor_sheet.bloch_sheet_tensor`, Phase 26: (0,0) harmonic, cell-scale near field included, evanescent kt > k0).
The cell is C2-symmetric, so psi and psi + 180 coincide; the square lattice adds a C4 symmetry used as a consistency check. Also stored: the uniform-cell phase constant beta(slot, a, psi) from the tensor
dispersion (`find_beta_tensor` with Y_s at kt = 1.065 k0 taken as kt-independent: an approximation near the root).

    python scripts/research/phase29_angle_lut.py      # results/phase29/angle_lut_v1.npz + angle_lut_v1.json (provenance), ~ minutes on all cores"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import datetime as dt  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
SLOTS = ((4.4, 1.2), (4.6, 1.6), (4.8, 1.6), (4.8, 2.0))
ANG = np.arange(0.0, 180.0, 15.0)
BETA0 = 1.065
MESH = 48


def job(args):
    s, ia = args
    from phase26_slotted_family import K0, SUB, F, make_cell

    from src.basis_functions import MeshSettings, RooftopBasis
    from src.dispersion import tensor_sheet as ts

    mesh = MeshSettings(MESH, MESH, MESH)
    cell = make_cell("circle", 6.0, SLOTS[s][0], SLOTS[s][1], ANG[ia])
    basis = RooftopBasis(cell, mesh.cells_x, mesh.cells_y)
    kt = BETA0 * K0
    y = np.zeros((len(ANG), 2, 2), dtype=complex)
    beta = np.full(len(ANG), np.nan)
    te = np.full(len(ANG), np.nan)
    for ip, psi in enumerate(np.radians(ANG)):
        y[ip] = ts.bloch_sheet_tensor(cell, F, mesh, kt * np.cos(psi), kt * np.sin(psi), basis=basis)
        try:
            b, ch = ts.find_beta_tensor(y[ip], psi, SUB, F)
            beta[ip], te[ip] = b / K0, float(ch) if np.isscalar(ch) else float(np.abs(ch))
        except Exception:  # no root: recorded as NaN, never invented
            pass
    return s, ia, y, beta, te


def main():
    from phase29_pointB_sweep import n_workers

    t0 = time.time()
    jobs = [(s, ia) for s in range(len(SLOTS)) for ia in range(len(ANG))]
    Y = np.zeros((len(SLOTS), len(ANG), len(ANG), 2, 2), dtype=complex)
    B = np.full((len(SLOTS), len(ANG), len(ANG)), np.nan)
    T = np.full_like(B, np.nan)
    with mp.get_context("spawn").Pool(min(len(jobs), n_workers(per_worker_gb=2.0, reserve_gb=8.0, cap=46))) as pool:
        for s, ia, y, b, te in pool.imap_unordered(job, jobs):
            Y[s, ia], B[s, ia], T[s, ia] = y, b, te
    # consistency: C4 of the square lattice, Y(a + 90, psi + 90) = R(90) Y(a, psi) R(90)^T
    R90 = np.array([[0.0, -1.0], [1.0, 0.0]])
    errs = []
    for s in range(len(SLOTS)):
        for ia in range(len(ANG)):
            for ip in range(len(ANG)):
                y2 = Y[s, (ia + 6) % 12, (ip + 6) % 12]
                errs.append(np.linalg.norm(y2 - R90 @ Y[s, ia, ip] @ R90.T) / max(np.linalg.norm(Y[s, ia, ip]), 1e-30))
    sym = {
        "c4_relative_error_median": float(np.median(errs)),
        "c4_relative_error_max": float(np.max(errs)),
        "reciprocity_max_asym": float(
            np.max(np.abs(Y[..., 0, 1] - Y[..., 1, 0]) / np.maximum(np.abs(Y[..., 0, 1]), 1e-30))
        ),
    }
    out = ROOT / "results/phase29"
    np.savez_compressed(
        out / "angle_lut_v1.npz", Y=Y, beta_over_k0=B, te_content=T, slots_mm=np.array(SLOTS), a_deg=ANG, psi_deg=ANG
    )
    head = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    (out / "angle_lut_v1.json").write_text(
        json.dumps(
            {
                "version": "v1",
                "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "BLOCH_SHEET_TENSOR_FROM_PERIODIC_MOM"],
                "generated": f"{dt.datetime.now():%Y-%m-%d %H:%M}",
                "git_head_at_generation": head,
                "source": "scripts/research/phase29_angle_lut.py -> src/dispersion/tensor_sheet.bloch_sheet_tensor (Phase 26)",
                "frequency_ghz": 10.54,
                "kt_over_k0": BETA0,
                "mesh": MESH,
                "cell": "slotted circle D 6.0 mm, 6.4 mm lattice, RT/duroid 5880 1.575 mm grounded",
                "axes": {
                    "Y": "[slot, a_i, psi_j, 2, 2] (siemens, x/y comps; slot along angle a from +x; kt along angle psi from +x)",
                    "slots_mm": [list(x) for x in SLOTS],
                    "a_deg": ANG.tolist(),
                    "psi_deg": ANG.tolist(),
                },
                "nan_beta_count": int(np.isnan(B).sum()),
                "symmetry_checks": sym,
                "seconds": time.time() - t0,
            },
            indent=1,
        )
    )
    print(sym, "nan beta", int(np.isnan(B).sum()), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
