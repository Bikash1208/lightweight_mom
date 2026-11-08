"""SUB6 slot-disc diode-inductance tolerance (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED).  Mesh 64, F = 2n = 128, RT5880 p11.6.

Rows: the two HFSS cells (g0.73_q0.8_W2.9, g1_q0.9_W2.9) + the 3 next-best rows by alpha_ON (g1_q0.9_W2.4, g0.8_q0.8_W3.6, g0.73_q0.8_W3.6).
Cases: ON_BAR50 Ls = 0.4..0.8 nH (R 3.0), ON R = 4.5 ohm @ 0.6 nH, OFF_BAR50 C = 0.12/0.15/0.20 pF (Rp 5k, L 0.6 nH).
One BRIDGE Z assembly per (row, direction, kt node) shared by all 9 cases (only the load matrix changes); nodes = the row's existing ckpt kt nodes,
extended adaptively (max 5) until every case's beta lies within 0.004 k0 of a node.  alpha: leakage_alpha_first_order_tensor, 25 deg (as the cell scan).

    python sub6_slotdisc_ls_tol.py run [--workers N]      -> results/sub6_slotdisc_ls_tol/{nodes_<id>.json, ls_tol.json, ls_tol.csv}
"""
import argparse
import json
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.basis_functions import MeshSettings, RooftopBasis
from src.frequency import free_space_wavenumber
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.solver import solve

OUT = G.ROOT / "results/sub6_slotdisc_ls_tol"
ROWS = ["RT5880_p11.6_g0.73_q0.8_W2.9", "RT5880_p11.6_g1_q0.9_W2.9", "RT5880_p11.6_g1_q0.9_W2.4", "RT5880_p11.6_g0.8_q0.8_W3.6", "RT5880_p11.6_g0.73_q0.8_W3.6"]
CASES = {}
for ls in (0.4, 0.5, 0.6, 0.7, 0.8):
    CASES[f"ON_Ls{ls:g}"] = dict(R=3.0, L=ls * 1e-9)
CASES["ON_R4.5_Ls0.6"] = dict(R=4.5, L=0.6e-9)
for c in (0.12, 0.15, 0.20):
    CASES[f"OFF_C{c:g}"] = dict(z_of_f=(lambda cc: (lambda f: S.z_off_bar50(f, C=cc * 1e-12)))(c))
GUESS = 1.06
MAXN = 5


def case_task(a):
    """(row dict, direction, kt) -> {case: 2x2 Y}; one BRIDGE assembly."""
    row, d, kt = a
    sub = S.make_sub(row["sub"])
    f = S.F
    k0 = free_space_wavenumber(f)
    cell = S.make_cell(row["p"], row["r"], row["L"], row["W"], sub, "ON_BAR50")
    n, fo = G.N, G.FO
    kx, ky = (kt * k0, 0.0) if d == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(sub, f, backing="pec")
    basis = RooftopBasis(cell, n, n)
    z0 = assemble_impedance_matrix(basis, green, kx, ky, fo)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for name, kw in CASES.items():
        ld = LumpedLoad(cell.geometry.diode_x, 0.0, length=S.G_C, width=S.DIODE_ALONG, direction="y", **kw)
        cur = solve(z0 + load_matrix(basis, [ld], green.omega), v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[name] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(row=row["id"], direction=d, kt=kt, Y=out)


def get_rows():
    rs = {r["id"]: r for r in G.rows(["RT5880"])}
    return [rs[i] for i in ROWS]


def nodes_file(rid):
    return OUT / f"nodes_{rid}.json"


def load(rid):
    p = nodes_file(rid)
    if not p.exists():
        return {}
    return {float(k): {d: {c: G.to_c(y) for c, y in cs.items()} for d, cs in v.items()} for k, v in json.load(open(p)).items()}


def save(rid, nodes):
    OUT.mkdir(parents=True, exist_ok=True)
    nodes_file(rid).write_text(json.dumps({str(k): G.jc(v) for k, v in nodes.items()}))


def summary(nodes, case, sub):
    # reuse the lib's state_summary: it indexes nodes[kt][direction][state]
    return S.state_summary(nodes, case, sub, guess=GUESS)


def main(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    rows = get_rows()
    state = {}
    for r in rows:
        nodes = load(r["id"])
        if not nodes:
            ck = G.load_nodes(G.OUT / "ckpt" / f"{r['id']}.json")
            state[r["id"]] = dict(pending=sorted(ck))  # start from the cell-scan nodes
        else:
            state[r["id"]] = dict(pending=[])
        state[r["id"]]["nodes"] = nodes
    t0 = time.time()
    with Pool(workers) as pool:
        for rnd in range(MAXN + 1):
            tasks = []
            for r in rows:
                st = state[r["id"]]
                for kt in st["pending"]:
                    for d in "xy":
                        tasks.append((r, d, float(kt)))
            if not tasks:
                break
            print(f"[{time.time() - t0:6.0f}s] round {rnd}: {len(tasks)} assemblies", flush=True)
            for o in pool.imap_unordered(case_task, tasks):
                state[o["row"]]["nodes"].setdefault(o["kt"], {})[o["direction"]] = o["Y"]
                if len(state[o["row"]]["nodes"].get(o["kt"], {})) == 2:
                    save(o["row"], state[o["row"]]["nodes"])
            for r in rows:
                st = state[r["id"]]
                st["pending"] = []
                if len(st["nodes"]) >= MAXN:
                    continue
                bs = [summary(st["nodes"], c, r["sub"]) for c in CASES]
                far = [b["beta_over_k0"] for b in bs if b["nearest_node_dist"] > 0.004]
                if far:
                    kt = float(np.median(far))
                    if all(abs(kt - k) > 3e-4 for k in st["nodes"]):
                        st["pending"] = [kt]
    res = []
    for r in rows:
        nodes = state[r["id"]]["nodes"]
        cs = {c: summary(nodes, c, r["sub"]) for c in CASES}
        # d(alpha)/d(Ls) at 0.6 nH: central (0.5, 0.7) and a linear least squares over 0.4-0.8
        ls = np.array([0.4, 0.5, 0.6, 0.7, 0.8])
        al = np.array([cs[f"ON_Ls{x:g}"]["alpha_np_per_m"] for x in ls])
        dcen = (al[3] - al[1]) / 0.2
        dfit = float(np.polyfit(ls, al, 1)[0])
        res.append(dict(id=r["id"], nodes=sorted(nodes), cases=cs, dalpha_dLs_central_np_per_m_per_nH=float(dcen), dalpha_dLs_linfit=dfit,
                        dalpha_dLs_rel_percent_per_0p1nH=float(100 * dcen * 0.1 / al[2])))
        print(r["id"], "dalpha/dLs", dcen, flush=True)
    prov = G.provenance(dict(stage="ls_tol", rows=ROWS, workers=workers, seconds=time.time() - t0))
    (OUT / "ls_tol.json").write_text(json.dumps(G.jc(dict(provenance=prov, results=res)), indent=1, default=str))
    lines = ["row,case,alpha_Np_per_m,beta_over_k0,absDelta_mS,Delta_re_mS,Delta_im_mS,extrapolated,dalpha_dLs_at_0p6nH_Np_per_m_per_nH"]
    for x in res:
        for c, v in x["cases"].items():
            dl = v["Delta_mS"]
            lines.append(f"{x['id']},{c},{v['alpha_np_per_m']:.4g},{v['beta_over_k0']:.4g},{abs(dl):.4g},{dl.real:.4g},{dl.imag:.4g},{v['extrapolated']},{x['dalpha_dLs_central_np_per_m_per_nH']:.4g}")
    (OUT / "ls_tol.csv").write_text("\n".join(lines) + "\n")
    print("DONE", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run"])
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    nw = a.workers or int(min(61, max(4, G.free_gb() * 0.7 / 1.6)))
    main(nw)
