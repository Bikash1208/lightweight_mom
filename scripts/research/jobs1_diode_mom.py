"""JOB E (cell level, MODEL_ONLY / NOT_HFSS_VALIDATED): island cell ISL_W2.9_q0.8_Li0.5_c0.6 (mesh 96, 8-cell diode gap = 0.967 mm, as the board table) with other diodes in the gap.
For every candidate diode (Z_off(f), Z_on(f) in the table below) and both states: the Bloch sheet tensor at kt = 1.017 and 1.097 (two nodes, directions along (x) and across (y) the slot) -> the scan's own self-consistent beta
(sub6_slotdisc_lib.state_summary, linear interpolation / extrapolation over the two kt nodes) -> Delta, alpha (leakage_alpha_first_order_tensor at 25 deg).  Output: alpha_OFF, alpha_ON, ratio (dB), |Delta|, beta.
    python jobs1_diode_mom.py <workers> <out.json>
The 2 x 2 (kt, direction) assemblies are shared by all diodes (one Z assembly each, then one solve per diode state).
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8"); os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("ISL_GAP_CELLS", "8")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
import numpy as np

import diag_full_tensor_lib as L
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.solver import solve

F = S.F
# name: (C_off [F], Rp_off [ohm], Ls [H], R_on [ohm]); iteration 1 = BAR50-02V only (user 2026-11-10), 5.8 GHz values OFF C 0.11 (0.105-0.12) || Rp 3-5k + 0.6 nH, ON 3 ohm + 0.6 nH
DIODES = {
    "BAR50-02V C0.15 Rp5k (scan model)": (0.15e-12, 5e3, 0.6e-9, 3.0),
    "BAR50-02V C0.11 Rp4k (iteration 1 nominal)": (0.11e-12, 4e3, 0.6e-9, 3.0),
    "BAR50-02V C0.12 Rp3k (worst OFF)": (0.12e-12, 3e3, 0.6e-9, 3.0),
    "BAR50-02V C0.105 Rp5k (best OFF)": (0.105e-12, 5e3, 0.6e-9, 3.0),
}


def zload(name, state):
    C, Rp, Ls, Ron = DIODES[name]
    if state == "OFF":
        return lambda f: 1j * 2 * math.pi * f * Ls + 1.0 / (1.0 / Rp + 1j * 2 * math.pi * f * C)
    return lambda f: Ron + 1j * 2 * math.pi * f * Ls


def node(a):
    kt, d = a
    t0 = time.time()
    j96 = json.load(open(ROOT / "results/sub6_final_board/mesh96/input_island_rows_n96_d789d9a.json"))
    cells = {r["row"]["id"]: r for r in j96["rows"]}["ISL_W2.9_q0.8_Li0.5_c0.6"]["row"]["cells"]
    assert L.GAP_CELLS == 8
    mk = lambda: L.island_cell(cells, L.P / 96, "OFF_BAR50")
    cell0 = mk()
    geo = cell0.geometry
    k0 = free_space_wavenumber(F)
    kx, ky = (kt * k0, 0.0) if d == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell0.substrate, F, backing="pec")
    basis = RooftopBasis(cell0, 96, 96)
    z0 = assemble_impedance_matrix(basis, green, kx, ky, L.FORDER)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for name in DIODES:
        for st in ("OFF", "ON"):
            ld = LumpedLoad(0.0, geo.diode_y, length=S.G_C, width=S.DIODE_ALONG, direction="y", z_of_f=zload(name, st))
            z = z0 + load_matrix(basis, [ld], green.omega)
            cur = solve(z, v)
            j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
            out[f"{name}|{st}"] = (j @ np.linalg.inv(np.eye(2) + g00 @ j)).tolist()
    return kt, d, out, time.time() - t0


def main():
    nw, outp = int(sys.argv[1]), sys.argv[2]
    jobs = [(kt, d) for kt in (1.017, 1.097) for d in "xy"]
    nodes = {}
    t0 = time.time()
    with Pool(nw) as p:
        for kt, d, out, sec in p.imap_unordered(node, jobs):
            for k, Y in out.items():
                nodes.setdefault(k, {}).setdefault(kt, {})[d] = np.array(Y)
            print(f"{time.time() - t0:.0f}s node kt {kt} dir {d} done ({sec:.0f}s)", flush=True)
    res = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED island cell with alternative diodes (mesh 96, gap 8 cells)", diodes={k: dict(C_off_pF=v[0] * 1e12, Rp_off_ohm=v[1], Ls_nH=v[2] * 1e9, R_on_ohm=v[3]) for k, v in DIODES.items()}, rows={})
    for name in DIODES:
        per = {}
        for st in ("OFF", "ON"):
            nd = {kt: {d: {st: nodes[f"{name}|{st}"][kt][d]} for d in "xy"} for kt in (1.017, 1.097)}
            r = S.state_summary(nd, st, L.SUB, guess=1.06)
            per[st] = dict(beta_over_k0=r["beta_over_k0"], alpha_np_per_m=r["alpha_np_per_m"], Delta_mS=[complex(r["Delta_mS"]).real, complex(r["Delta_mS"]).imag], abs_Delta_mS=abs(complex(r["Delta_mS"])),
                           Y_par_mS=[complex(r["Y_par_mS"]).real, complex(r["Y_par_mS"]).imag], Y_perp_mS=[complex(r["Y_perp_mS"]).real, complex(r["Y_perp_mS"]).imag], extrapolated=r["extrapolated"])
        a_off, a_on = per["OFF"]["alpha_np_per_m"], per["ON"]["alpha_np_per_m"]
        per["alpha_on_over_off"] = a_on / a_off if a_off > 0 else None
        per["suppression_dB"] = 10 * math.log10(a_off / a_on) if a_on > 0 and a_off > 0 else None
        res["rows"][name] = per
        print(name, "alpha OFF %.3g ON %.3g  ratio %.2f dB" % (a_off, a_on, per["suppression_dB"] or float("nan")), "beta", round(per["OFF"]["beta_over_k0"], 3), round(per["ON"]["beta_over_k0"], 3), flush=True)
    Path(outp).write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
