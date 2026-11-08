"""SUB6 island cell: mesh convergence 64 / 96 / 128 for two rows (MODEL_ONLY).  Built on sub6_island_cell.py (compute/sub6-island-cell, unchanged); takes over the unrun script of commit f3f6cd4.

Rows ISL_W2.9_q0.8_Li0.5_c0.6 (board row) and ISL_W2.9_q0.85_Li0.5_c0.5 (refine row), states OFF_BAR50 / ON_BAR50, 5.8 GHz.
Mesh cell dx = p / N, N even (96, 128; 64 = the island-scan baseline) so the x = 0 mirror plane is a cell EDGE: slot, island and diode sheet are mirror-symmetric in x with 2 n cells.  Integer cell counts per mesh
(as in the island scan, re-derived): slot half length / width = floor(L/2/dx + 1/2), island half length = round(Li/2/dx), clearance (-y) = round(c/dx), +y diode gap = round(1.00/dx) (mesh 64 keeps the scan's 5 cells);
the realised dimensions are written for every mesh.  Floquet order F = 2N.
Node set: fixed seed nodes kt = 1.004, 1.054, 1.104 (they bracket beta_ON 1.005 and beta_OFF 1.104 of the 64 scan, which were EXTRAPOLATED there: the scan stopped at nodes 1.054/1.0545/1.06), then adaptive nodes (<= 5) until every state's
beta is within 0.004 k0 of a node.  Lean task: the load matrix enters only through its non-zero entries (no second n x n copy), LU in place.  Resonance proxy at meshes 64 and 96 only (kt = 0, 5.0 / 6.6 GHz).
    python sub6_island_mesh128.py list | run --mesh N [--workers W] | collect      (threads per task: set OPENBLAS_NUM_THREADS before the call)
"""
import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

for _v, _d in (("OMP_NUM_THREADS", "1"), ("MKL_NUM_THREADS", "1")):
    os.environ.setdefault(_v, _d)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import scipy.linalg as sla

import sub6_island_cell as IC
import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice

OUT = IC.ROOT / "results" / "sub6_island_mesh128"
MM, P, R, SUB = 1e-3, IC.P, IC.R_ANCHOR, IC.SUB
MESHES = (64, 96, 128)
ROWS = [dict(id="ISL_W2.9_q0.8_Li0.5_c0.6", W=2.9, q=0.8, li=0.5, c=0.6), dict(id="ISL_W2.9_q0.85_Li0.5_c0.5", W=2.9, q=0.85, li=0.5, c=0.5)]
STATES = IC.ISL_STATES
SEEDS = (1.004, 1.054, 1.104)
RES_F = (5.0e9, 6.6e9)
MAX_NODES, TOL = 5, 0.004


@dataclass(frozen=True)
class Isl(IC.SlotDiscIsland):
    dx_mm: float = 0.0
    gap: int = 5

    @property
    def dxm(self):
        return self.dx_mm * MM

    @property
    def y_top_island(self):
        return (self.n_w - self.gap) * self.dxm


def cells(N, row):
    dx = P / N
    L = 2 * R * row["q"]
    n_l = int(math.floor(L / 2 / dx + 0.5))
    n_w = int(math.floor(row["W"] / 2 / dx + 0.5))
    n_i = max(1, int(round(row["li"] * L / 2 / dx)))
    cb = max(1, int(round(row["c"] / dx)))
    gap = 5 if N == 64 else max(1, int(round(1.0 / dx)))
    return n_l, n_w, n_i, cb, gap


def real(N, row):
    dx = P / N
    n_l, n_w, n_i, cb, gap = cells(N, row)
    ns = sum(1 for k in range(-40, 40) if abs((k + 0.5) * dx) <= S.DIODE_ALONG / MM / 2)
    return dict(mesh=N, dx_mm=dx, slot_L_mm=2 * n_l * dx, slot_W_mm=2 * n_w * dx, island_Li_mm=2 * n_i * dx, island_Wi_mm=(2 * n_w - gap - cb) * dx, clear_minus_y_mm=cb * dx, gap_plus_y_mm=gap * dx,
                clear_short_ends_mm=(n_l - n_i) * dx, sheet_metal_width_mm=ns * dx, cells=[n_l, n_w, n_i, cb, gap],
                nominal=dict(slot_L_mm=2 * R * row["q"], slot_W_mm=row["W"], island_Li_mm=row["li"] * 2 * R * row["q"], clear_mm=row["c"], gap_mm=1.0, sheet_along_mm=S.DIODE_ALONG / MM))


def make_cell(N, c, state):
    n_l, n_w, n_i, cb, gap = c
    geo = Isl(R * MM, n_l, n_w, n_i, cb, True, dx_mm=P / N, gap=gap)
    ld = S.STATES[state][1]
    return UnitCell(RectangularLattice(P * MM, P * MM), geo, S.make_sub(SUB), Metal(), backing="pec",
                    loads=(LumpedLoad(0.0, geo.diode_y, length=S.G_C, width=S.DIODE_ALONG, direction="y", **ld),))


def task(a):
    N, c, direction, kt, f = a
    t0 = time.time()
    cell0 = make_cell(N, c, STATES[0])
    k0 = free_space_wavenumber(f)
    kx, ky = (kt * k0, 0.0) if direction == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell0.substrate, f, backing="pec")
    basis = RooftopBasis(cell0, N, N)
    ent = {}
    for st in STATES:  # load entries only (the dense n x n load matrix is built and dropped before Z exists)
        lm = load_matrix(basis, list(make_cell(N, c, st).loads), green.omega)
        nz = np.nonzero(lm)
        ent[st] = (nz, lm[nz])
        del lm
    z0 = assemble_impedance_matrix(basis, green, kx, ky, 2 * N)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for i, st in enumerate(STATES):
        nz, vals = ent[st]
        if i == 0:
            z = z0.copy()
            z[nz] += vals
        else:
            z0[nz] += vals
            z = z0
        cur = sla.solve(z, v, overwrite_a=True, check_finite=False)
        del z
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[st] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(direction=direction, kt=kt, f=f, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def ck(N, rid):
    return OUT / "ckpt" / f"m{N}_{rid}.json"


def next_nodes(nodes):
    """kt of beta for every state whose beta is further than TOL from a node (one new node per such state), within MAX_NODES."""
    bs = []
    for st in STATES:
        r = S.state_summary(nodes, st, SUB, guess=1.06)
        if r["nearest_node_dist"] > TOL:
            bs.append(r["beta_over_k0"])
    new = []
    for b in sorted(bs):
        if all(abs(b - k) > 3e-4 for k in list(nodes) + new) and len(nodes) + len(new) < MAX_NODES:
            new.append(float(b))
    return new


def run(N, workers):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ckpt").mkdir(exist_ok=True)
    jobs = {r["id"]: cells(N, r) for r in ROWS}
    nodes = {k: G.load_nodes(ck(N, k)) for k in jobs}
    resf = OUT / f"resonance_m{N}.json"
    resv = json.loads(resf.read_text()) if resf.exists() else {}
    t0 = time.time()

    def log(m):
        print(f"[{time.time() - t0:6.0f}s] mesh {N}: {m}", flush=True)

    with Pool(workers) as pool:
        for rnd in range(4):
            todo = []
            for k, c in jobs.items():
                kts = SEEDS if not nodes[k] else next_nodes(nodes[k])
                todo += [(k, (N, c, d, kt, S.F)) for kt in kts for d in "xy"]
            res_t = []
            if rnd == 0 and N in (64, 96):
                res_t = [(k, (N, c, "x", 0.0, fr)) for k, c in jobs.items() for fr in RES_F if f"{k}|{fr:g}" not in resv]
            allt = todo + res_t
            if not allt:
                log("nodes converged")
                break
            log(f"round {rnd}: {len(todo)} node assemblies + {len(res_t)} resonance, {workers} workers")
            for i, (t, res) in enumerate(zip(allt, pool.imap(task, [x[1] for x in allt], chunksize=1))):
                k = t[0]
                if res["kt"] == 0.0:
                    resv[f"{k}|{res['f']:g}"] = {st: G.jc(y) for st, y in res["Y"].items()}
                else:
                    nodes[k].setdefault(res["kt"], {}).setdefault(res["direction"], {}).update(res["Y"])
                log(f"{i + 1}/{len(allt)} {k} {res['direction']} kt={res['kt']:.4f} f={res['f'] / 1e9:.1f} {res['seconds']:.0f}s unknowns {res['unknowns']}")
                for k2 in jobs:
                    G.save_nodes(ck(N, k2), nodes[k2])
                resf.write_text(json.dumps(resv))
    (OUT / f"run_meta_m{N}.json").write_text(json.dumps(dict(mesh=N, workers=workers, seconds=time.time() - t0, realised={r["id"]: real(N, r) for r in ROWS}, nodes={k: sorted(v) for k, v in nodes.items()}), indent=1, default=float))
    print("RUN_DONE", N, flush=True)


def rel(a, b):
    return b / a - 1.0


def collect():
    out = {}
    for r in ROWS:
        d = dict(realised={}, mesh={})
        for N in MESHES:
            if not ck(N, r["id"]).exists():
                continue
            nd = G.load_nodes(ck(N, r["id"]))
            if not nd:
                continue
            d["realised"][N] = real(N, r)
            ms = {st: S.state_summary(nd, st, SUB, guess=1.06) for st in STATES}
            ms["nodes"] = sorted(nd)
            ro = ms["ON_BAR50"]["Delta_mS"] / ms["OFF_BAR50"]["Delta_mS"]
            ms["rho"] = ro
            ms["abs_rho_minus_1"] = float(abs(ro - 1.0))
            ms["gate"] = "OK" if (ms["ON_BAR50"]["alpha_np_per_m"] >= 1.0 or ms["OFF_BAR50"]["alpha_np_per_m"] >= 1.0) and abs(ro - 1.0) >= 1.0 else "NOT_OK"
            rf = OUT / f"resonance_m{N}.json"
            ms["resonance_flip"] = None
            if rf.exists():
                rv = json.loads(rf.read_text())
                a_, b_ = rv.get(f"{r['id']}|{RES_F[0]:g}"), rv.get(f"{r['id']}|{RES_F[1]:g}")
                if a_ and b_:
                    ms["resonance_flip"] = [f"{st}:{'xy'[i]}" for st in STATES for i in (0, 1) if np.sign(G.to_c(a_[st])[i, i].imag) != np.sign(G.to_c(b_[st])[i, i].imag)]
            d["mesh"][N] = ms
        out[r["id"]] = d
    (OUT / "mesh128_check.json").write_text(json.dumps(G.jc(out), indent=1, default=float))
    stop = []
    for rid, d in out.items():
        print(rid)
        for N, ms in d["mesh"].items():
            rl = d["realised"][N]
            print(f"  mesh {N}: dx {rl['dx_mm']:.4f} slot {rl['slot_L_mm']:.3f}x{rl['slot_W_mm']:.3f} Li {rl['island_Li_mm']:.3f} Wi {rl['island_Wi_mm']:.3f} clear {rl['clear_minus_y_mm']:.3f} gap {rl['gap_plus_y_mm']:.3f} sheet {rl['sheet_metal_width_mm']:.3f} nodes {[round(x, 4) for x in ms['nodes']]}")
            for st in STATES:
                s = ms[st]
                print(f"     {st:9s} alpha {s['alpha_np_per_m']:.4f} beta/k0 {s['beta_over_k0']:.4f} |Delta| {abs(s['Delta_mS']):.3f} mS  A_perp {abs(s['A_perp_mS']):.2e} extrap {s['extrapolated']} node dist {s['nearest_node_dist']:.4f}")
            print(f"     rho {ms['rho']:.3f}  |rho-1| {ms['abs_rho_minus_1']:.3f}  gate {ms['gate']}  resonance_flip {ms.get('resonance_flip')}")
        for a, b in ((64, 96), (96, 128), (64, 128)):
            if a in d["mesh"] and b in d["mesh"]:
                for st in STATES:
                    x, y = d["mesh"][a][st], d["mesh"][b][st]
                    ch = dict(alpha=rel(x["alpha_np_per_m"], y["alpha_np_per_m"]), beta=rel(x["beta_over_k0"], y["beta_over_k0"]), Delta=rel(abs(x["Delta_mS"]), abs(y["Delta_mS"])))
                    flag = " <-- > 10 % (RULE 3)" if max(abs(ch["alpha"]), abs(ch["Delta"])) > 0.10 else ""
                    print(f"  {a}->{b} {st:9s} alpha {100 * ch['alpha']:+.1f}%  beta {100 * ch['beta']:+.2f}%  |Delta| {100 * ch['Delta']:+.1f}%{flag}")
                    if flag and (a, b) == (96, 128):
                        stop.append((rid, st, ch))
    print("RULE3_HARD_STOP" if stop else "rule 3 (96 -> 128 > 10 % in alpha or |Delta|): not triggered", stop)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["list", "run", "collect"])
    ap.add_argument("--mesh", type=int, default=128)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.stage == "list":
        for r in ROWS:
            for N in MESHES:
                print(r["id"], json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in real(N, r).items() if k != "nominal"}))
    elif a.stage == "run":
        run(a.mesh, a.workers)
    else:
        collect()
