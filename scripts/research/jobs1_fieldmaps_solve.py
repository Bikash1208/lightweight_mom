"""JOB 1 A (model only): near-field maps of the Bloch mode of the island cell (OFF/ON) and the static NOCUT cell at kt = beta (along the slot = x, across = y).
Same MoM path as the scans (mesh 96 island / 64 static, Floquet order 128).  Fields of the solved current sheet (impressed (0,0) field not included) from the layered Green function:
|E_t| and |E_z| at z = +0.1 mm above the top copper, and |E_z| in the ground-plane plane z = -h (E_t = 0 on the PEC; E_z there is what couples into an anti-pad).
E_z from Ampere: E_z = (ky Hx - kx Hy) / (omega eps)  (harmonic exp(-j(kx x + ky y)); eps = eps0 in air, eps0*eps_r in the substrate).
    python jobs1_fieldmaps_solve.py <workers>   -> results/jobs1/fieldmaps/maps_<cell>_<state>_<dir>.npz
MODEL_ONLY / NOT_HFSS_VALIDATED.
"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "5"); os.environ.setdefault("OMP_NUM_THREADS", "5")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
import numpy as np

import diag_full_tensor_lib as L
from src.basis_functions import RooftopBasis
from src.constants import EPS0
from src.frequency import free_space_wavenumber
from src.green_function import LayeredGreenFunction
from src.loads import load_matrix
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.solver import solve

OUT = ROOT / "results" / "jobs1" / "fieldmaps"
NG = 96
ZTOP = 0.1e-3
ORD = 128


def job(a):
    kind, spec, n, state, direction, kt, tag = a
    t0 = time.time()
    k0 = free_space_wavenumber(L.F)
    kx0, ky0 = (kt * k0, 0.0) if direction == "x" else (0.0, kt * k0)
    cell = L.island_cell(spec["cells"], L.P / n, state) if kind == "island" else L.static_cell(spec)
    green = LayeredGreenFunction(cell.substrate, L.F, backing="pec")
    basis = RooftopBasis(cell, n, n)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, ORD)
    if cell.loads:
        z = z + load_matrix(basis, list(cell.loads), green.omega)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    col = 0 if direction == "x" else 1          # E_stack along the propagation direction (TM-like)
    c = cur[:, col]
    P = L.P * 1e-3
    m = np.arange(-ORD, ORD + 1)
    MX, MY = np.meshgrid(m, m, indexing="ij")
    kx = (kx0 + 2 * np.pi * MX / P).ravel(); ky = (ky0 + 2 * np.pi * MY / P).ravel()
    J = np.zeros((kx.size, 2), complex)
    for lo in range(0, kx.size, 3000):
        sl = slice(lo, lo + 3000)
        J[sl] = basis.current_coefficients(c, kx[sl], ky[sl])
    h = cell.substrate.thickness
    w = green.omega
    eps_sub = EPS0 * cell.substrate.eps_complex
    out = {}
    for name, zz, eps in (("top", ZTOP, EPS0), ("gnd", -h, eps_sub)):
        et, ht = green.sheet_fields(kx, ky, J, np.full(kx.size, zz))
        ez = (ky * ht[:, 0] - kx * ht[:, 1]) / (w * eps)
        out[name] = (et, ez)
    xs = (np.arange(NG) - NG // 2) * P / NG
    ax = np.exp(-1j * 2 * np.pi * np.outer(xs, m) / P)
    def real(F):
        return ax @ F.reshape(2 * ORD + 1, 2 * ORD + 1) @ ax.T       # [ix, iy]
    res = dict(Etx_top=real(out["top"][0][:, 0]), Ety_top=real(out["top"][0][:, 1]), Ez_top=real(out["top"][1]), Ez_gnd=real(out["gnd"][1]),
               Etx_gnd=real(out["gnd"][0][:, 0]), Ety_gnd=real(out["gnd"][0][:, 1]))
    fn = OUT / f"maps_{tag}_{state}_{direction}.npz"
    np.savez_compressed(fn, xs_mm=xs * 1e3, kt=kt, mesh=n, seconds=time.time() - t0, **res)
    return fn.name, time.time() - t0


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    j96 = json.load(open(ROOT / "results/sub6_final_board/mesh96/input_island_rows_n96_d789d9a.json"))
    isl = {r["row"]["id"]: r for r in j96["rows"]}["ISL_W2.9_q0.8_Li0.5_c0.6"]
    boff, bon = isl["states"]["OFF_BAR50"]["beta_over_k0"], isl["states"]["ON_BAR50"]["beta_over_k0"]
    rows = {r["row"]["id"]: r for r in (json.loads(l) for l in open(ROOT / "results/sub6_slotdisc/rows_RT5880.jsonl") if l.strip())}
    srow = rows["RT5880_p11.6_g1_q0.8_W3.6"]
    sb = srow["states"]["NOCUT"]["beta_over_k0"]
    (OUT / "meta.json").write_text(json.dumps(dict(island_cells=isl["row"]["cells"], island_real=isl["row"]["real"], beta_off=boff, beta_on=bon, static_row=srow["row"], static_beta=sb), indent=1))
    spec = dict(cells=isl["row"]["cells"])
    jobs = [("island", spec, 96, "OFF_BAR50", d, boff, "island") for d in "xy"] + [("island", spec, 96, "ON_BAR50", d, bon, "island") for d in "xy"] + [("static", srow, 64, "NOCUT", d, sb, "static") for d in "xy"]
    with Pool(nw) as p:
        for r in p.imap_unordered(job, jobs):
            print(r, flush=True)
    print("ALL DONE", flush=True)
