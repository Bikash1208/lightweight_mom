"""SUB6 slot-disc HFSS package: MoM predictions of the DIRECTIONAL bound-mode dispersion beta_x(f), beta_y(f) per cell / state (MODEL_ONLY).

For each cell (CELL1, CELL2), frequency f, class (NOCUT, BRIDGE = OFF_BAR50 + ON_BAR50 sharing one Z assembly) and propagation direction d (x = along the slot,
y = across it) the sheet Y_dd(kt) is computed at adaptive kt nodes (own node set per direction) and the bound-mode condition is solved with the TM sheet admittance
Y_xx (d = x) or Y_yy (d = y) alone: beta_d = find_beta_from_sheet_admittance(j Im Y_dd(kt = beta_d)).  Phase per cell phi = beta_d * p.
No frozen code is changed: the lib functions are imported from sub6_slotdisc_lib.
    python sub6_slotdisc_hfss_pred.py run [--workers N]      (checkpoint: results/sub6_slotdisc/hfss_pred/nodes.json)
    python sub6_slotdisc_hfss_pred.py table                  (-> pred_table.json)
"""
import json
import os
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sub6_slotdisc_lib as lib
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.frequency import free_space_wavenumber

OUT = ROOT / "results" / "sub6_slotdisc" / "hfss_pred"
SUB = ("rt5880", 2.2, 0.0009, 3.175)
P = 11.6
CELLS = {"CELL1": dict(gap=0.73, ratio=0.8, W=2.9), "CELL2": dict(gap=1.0, ratio=0.9, W=2.9)}
FREQS = [5.2, 5.5, 5.8, 6.1, 6.4]  # GHz
N_MESH, FORDER = 64, 128
TOL, MAXNODES, ROUNDS = 0.004, 4, 3


def cell_spec(name):
    c = CELLS[name]
    r = round((P - c["gap"]) / 2, 6)
    return (SUB, P, r, round(c["ratio"] * 2 * r, 6), c["W"])


def key(cell, f, cls, d):
    return f"{cell}|{f}|{cls}|{d}"


def beta_dir(nodes, st, d, f):
    """Self-consistent beta/k0 from Y_dd(kt) of direction d; nodes {kt: {state: 2x2}}."""
    s = lib.make_sub(SUB)
    k0 = free_space_wavenumber(f * 1e9)
    ks = sorted(nodes)
    i0 = 0 if d == "x" else 1

    def y(kt):
        ys = [nodes[k][st][i0, i0] for k in ks]
        if len(ks) == 1:
            return ys[0]
        i = int(np.clip(np.searchsorted(ks, kt) - 1, 0, len(ks) - 2))
        w = (kt - ks[i]) / (ks[i + 1] - ks[i])
        return (1 - w) * ys[i] + w * ys[i + 1]

    kt = ks[0] if len(ks) == 1 else 1.08
    for _ in range(60):
        b = find_beta_from_sheet_admittance(1j * y(kt).imag, s, f * 1e9).beta_rad_per_m / k0
        if abs(b - kt) < 1e-5:
            break
        kt = 0.5 * (kt + b)
    return float(b), complex(y(b)), min(abs(b - k) for k in ks)


def load():
    p = OUT / "nodes.json"
    if not p.exists():
        return {}
    raw = json.loads(p.read_text())
    return {k: {float(kt): {st: np.array(v[0]) + 1j * np.array(v[1]) for st, v in sts.items()} for kt, sts in nd.items()} for k, nd in raw.items()}


def save(nodes):
    raw = {k: {repr(kt): {st: [y.real.tolist(), y.imag.tolist()] for st, y in sts.items()} for kt, sts in nd.items()} for k, nd in nodes.items()}
    (OUT / "nodes.json").write_text(json.dumps(raw))


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    nodes = load()
    keys = [(c, f, cls, d) for c in CELLS for f in FREQS for cls in ("NOCUT", "BRIDGE") for d in ("x", "y")]
    t00 = time.time()
    with Pool(workers) as pool:
        for rnd in range(ROUNDS + 1):
            todo = []
            for (c, f, cls, d) in keys:
                k = key(c, f, cls, d)
                nd = nodes.setdefault(k, {})
                if not nd:
                    new = [1.10]
                else:
                    bs = []
                    for st in lib.CLASS_STATES[cls]:
                        b, _, dist = beta_dir(nd, st, d, f)
                        if dist > TOL:
                            bs.append(b)
                    new = []
                    for b in sorted(bs):
                        if all(abs(b - k2) > 3e-4 for k2 in list(nd) + new) and len(nd) + len(new) < MAXNODES:
                            new.append(b)
                todo += [((cell_spec(c), cls, d, kt, N_MESH, FORDER, f * 1e9), k) for kt in new]
            if not todo:
                print("round", rnd, "no new nodes: converged", flush=True)
                break
            print("round", rnd, "tasks", len(todo), f"t={time.time() - t00:.0f}s", flush=True)
            res = pool.map(lib.node_task, [a for a, _ in todo], chunksize=1)
            for (a, k), r in zip(todo, res):
                nodes[k][a[3]] = r["Y"]
                print("  ", k, f"kt {a[3]:.4f}", f"{r['seconds']:.0f}s", flush=True)
            save(nodes)
    save(nodes)
    print("done", f"{time.time() - t00:.0f}s", flush=True)


def table():
    nodes = load()
    out = {}
    for c in CELLS:
        for st, cls in (("NOCUT", "NOCUT"), ("OFF_BAR50", "BRIDGE"), ("ON_BAR50", "BRIDGE")):
            for d in ("x", "y"):
                rows = []
                for f in FREQS:
                    b, y, dist = beta_dir(nodes[key(c, f, cls, d)], st, d, f)
                    rows.append(dict(f_GHz=f, beta_over_k0=b, phase_deg=float(np.degrees(b * free_space_wavenumber(f * 1e9) * P * 1e-3)),
                                     Y_dd_mS=[y.real * 1e3, y.imag * 1e3], nearest_node_dist=dist, n_nodes=len(nodes[key(c, f, cls, d)])))
                out[f"{c}|{st}|{d}"] = rows
    (OUT / "pred_table.json").write_text(json.dumps(out, indent=1))
    for k, rows in out.items():
        print(k, [(r["f_GHz"], round(r["beta_over_k0"], 3), round(r["phase_deg"], 1), round(r["nearest_node_dist"], 4)) for r in rows])


if __name__ == "__main__":
    if sys.argv[1] == "run":
        w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 40
        run(w)
    else:
        table()
