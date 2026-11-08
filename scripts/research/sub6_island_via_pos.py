"""SUB6 island cell: RF surface-current map |J| from the MoM solution (OFF and ON, 5.8 GHz) and the bias-via position of minimum current (MODEL_ONLY).

MoM: the mesh-64 island cell of sub6_island_cell (unchanged), RT5880 3.175, 5.8 GHz.  For every state (OFF_BAR50, ON_BAR50) and every propagation direction d (x along the slot, y across it) the cell is
driven at the bound-mode wavenumber kt = beta_d(5.8 GHz) (from results/sub6_island_hfss/hfss_pred/pred_table.json, linear interpolation of the two bracketing frequencies) by the unit tangential field along d
(the TM drive), j = Z^-1 v with the BAR50 sheet in; |J(x, y)| = sqrt(|Jx|^2 + |Jy|^2) is evaluated on a 0.05 mm grid with RooftopBasis.evaluate (Bloch current).  Each (state, direction) map is normalised to
its own peak |J| over the metal of the cell (so OFF and ON are compared as fractions of their own peak; the absolute level under the unit drive is also stored).
Via zones (HFSS nominal geometry, mm): island x in [-X_ISL + 0.3, X_ISL - 0.3], y in [-0.55, 0.15] (0.3 mm from the island edges); disc = disc (r 5.435) minus slot, at least 0.3 mm from the rim and from the
slot edges.  Score = max over {OFF, ON} (and over both directions: the HFSS runs are along x and across y; the across-y-only choice is reported too) of the normalised |J|.
    python sub6_island_via_pos.py run | pick
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_island_cell as IC
import sub6_island_hfss_package as PK
import sub6_island_hfss_pred as PR
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.green_function import LayeredGreenFunction
from src.loads import load_matrix
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.solver import solve

OUT = PK.RES / "via_pos"
GRID = 0.05
STATES = ("OFF_BAR50", "ON_BAR50")


def beta_58(state, d):
    t = json.load(open(PK.RES / "hfss_pred" / "pred_table.json"))
    key = f"ISL|{state}|{d}"
    rows = t[key]
    f = [r["f_GHz"] for r in rows]
    b = [r["beta_over_k0"] for r in rows]
    return float(np.interp(5.8, f, b))


def points():
    xs = np.arange(-5.8 + GRID / 2, 5.8, GRID)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    R = np.hypot(X, Y)
    slot = (np.abs(X) <= PK.L / 2) & (np.abs(Y) <= PK.W / 2)
    isl = (np.abs(X) <= PK.X_ISL) & (Y >= PK.Y_ISL0) & (Y <= PK.Y_ISL1)
    disc = (R <= PK.R) & ~slot
    return X, Y, disc, isl


def task(a):
    state, d, kt = a
    row = [r for r in json.load(open(IC.ROOT / "results" / "sub6_island_cell" / "island_rows.json"))["rows"] if r["row"]["id"] == PR.ROW_ID][0]["row"]
    cells = row["cells"]
    f = 5.8e9
    cell = IC.make_island_cell(IC.R_ANCHOR, cells, state, True)
    k0 = free_space_wavenumber(f)
    kx, ky = (kt * k0, 0.0) if d == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(cell.substrate, f, backing="pec")
    basis = RooftopBasis(cell, IC.N_MESH, IC.N_MESH)
    z = assemble_impedance_matrix(basis, green, kx, ky, IC.FORDER)
    if cell.loads:
        z = z + load_matrix(basis, list(cell.loads), green.omega)
    e = np.array([1.0, 0.0]) if d == "x" else np.array([0.0, 1.0])
    v = assemble_excitation_vector(basis, e, kx, ky)
    cur = solve(z, v[:, None])[:, 0] if v.ndim == 1 else solve(z, v)[:, 0]
    X, Y, disc, isl = points()
    m = (disc | isl).ravel()
    xm, ym = X.ravel()[m] * 1e-3, Y.ravel()[m] * 1e-3
    J = np.zeros(m.sum())
    chunk = 3000
    for i in range(0, xm.size, chunk):
        j = basis.evaluate(xm[i:i + chunk], ym[i:i + chunk], cur, kx, ky)
        J[i:i + chunk] = np.hypot(np.abs(j[..., 0]), np.abs(j[..., 1]))
    full = np.full(X.size, np.nan)
    full[m] = J
    return dict(state=state, d=d, kt=kt, J=full.reshape(X.shape), peak=float(np.nanmax(J)))


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    jobs = [(st, d, beta_58(st, d)) for st in STATES for d in ("x", "y")]
    with Pool(len(jobs)) as p:
        res = p.map(task, jobs)
    X, Y, disc, isl = points()
    np.savez_compressed(OUT / "J_maps.npz", X=X, Y=Y, disc=disc, isl=isl, **{f"J_{r['state']}_{r['d']}": r["J"] for r in res}, **{f"peak_{r['state']}_{r['d']}": r["peak"] for r in res}, **{f"kt_{r['state']}_{r['d']}": r["kt"] for r in res})
    print(f"done {time.time() - t0:.0f}s", [(r["state"], r["d"], round(r["kt"], 4), f"{r['peak']:.3g}") for r in res], flush=True)


def pick():
    D = np.load(OUT / "J_maps.npz")
    X, Y, disc, isl = D["X"], D["Y"], D["disc"], D["isl"]
    maps = {(st, d): D[f"J_{st}_{d}"] / D[f"peak_{st}_{d}"] for st in STATES for d in ("x", "y")}
    absmaps = {(st, d): D[f"J_{st}_{d}"] for st in STATES for d in ("x", "y")}
    # via zones
    marg = 0.3
    zone_i = isl & (np.abs(X) <= PK.X_ISL - marg) & (Y >= PK.Y_ISL0 + marg) & (Y <= PK.Y_ISL1 - marg)
    dx = np.maximum(np.abs(X) - PK.L / 2, 0.0)
    dy = np.maximum(np.abs(Y) - PK.W / 2, 0.0)
    d_slot = np.hypot(dx, dy)  # distance to the slot rectangle (0 inside it)
    zone_d = disc & (np.hypot(X, Y) <= PK.R - marg) & (d_slot >= marg)

    def score(sel, dirs):
        s = np.nanmax(np.stack([maps[(st, d)] for st in STATES for d in dirs]), axis=0)
        return np.where(sel, s, np.inf)

    out = {}
    for tag, dirs in (("both_directions", ("x", "y")), ("across_y_only", ("y",)), ("along_x_only", ("x",))):
        res = {}
        for name, zone in (("island", zone_i), ("disc", zone_d)):
            sc = score(zone, dirs)
            i = np.unravel_index(np.argmin(sc), sc.shape)
            res[name] = dict(x_mm=float(X[i]), y_mm=float(Y[i]), score=float(sc[i]),
                             J_norm={f"{st}|{d}": float(maps[(st, d)][i]) for st in STATES for d in ("x", "y")}, J_abs_unit_drive={f"{st}|{d}": float(absmaps[(st, d)][i]) for st in STATES for d in ("x", "y")})
        out[tag] = res
    cur = {}
    for name, (x, y) in (("island", (0.0, PK.Y_VIA_I)), ("disc", (0.0, PK.Y_VIA_D))):
        i = (int(round((x + 5.8 - GRID / 2) / GRID)), int(round((y + 5.8 - GRID / 2) / GRID)))
        cur[name] = dict(x_mm=x, y_mm=y, score_both=float(np.nanmax([maps[k][i] for k in maps])), J_norm={f"{st}|{d}": float(maps[(st, d)][i]) for st in STATES for d in ("x", "y")},
                         J_abs_unit_drive={f"{st}|{d}": float(absmaps[(st, d)][i]) for st in STATES for d in ("x", "y")})
    # regional statistics: max-normalised |J| over each metal piece
    stats = {}
    for name, sel in (("island", isl), ("disc", disc)):
        stats[name] = dict(median_over_piece_max_state_dir=float(np.nanmedian(np.where(sel, np.nanmax(np.stack(list(maps.values())), axis=0), np.nan))))
    json.dump(dict(chosen=out, current=cur, stats=stats, via_zones=dict(margin_mm=marg), norm="each (state, direction) map divided by its own peak |J| over the metal"), open(OUT / "via_positions.json", "w"), indent=1)
    for tag, res in out.items():
        print(tag, {k: (round(v["x_mm"], 2), round(v["y_mm"], 2), round(v["score"], 4)) for k, v in res.items()})
    print("current", {k: (v["x_mm"], v["y_mm"], round(v["score_both"], 4)) for k, v in cur.items()})
    print(stats)


if __name__ == "__main__":
    {"run": run, "pick": pick}[sys.argv[1]]()
