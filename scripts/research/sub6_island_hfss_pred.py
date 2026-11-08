"""SUB6 island-cell HFSS package: MoM predictions of the DIRECTIONAL bound-mode dispersion beta_x(f), beta_y(f) of the island cell (MODEL_ONLY).

Cell ISL = row ISL_W2.9_q0.8_Li0.5_c0.6 of the island scan (RT5880 3.175, p 11.6, r 5.435, slot 8.696 x 2.9 along x, island Li 4.348, clearance 0.6 nominal, +y gap 1.00 with one
BAR50 sheet).  MoM geometry = the island scan's mesh-64 mask (mesh-snapped dimensions, see sub6_island_cell.py).  States: OFF_BAR50, ON_BAR50 (one assembly) and ISLAND_ABSENT
(REF = closed slotted disc).  Same method as sub6_slotdisc_hfss_pred.py: beta_d from the TM sheet admittance Y_dd(kt = beta_d) alone (d = x along the slot, y across), adaptive kt nodes,
5 frequencies 5.2-6.4 GHz.  Uses sub6_island_cell.node_task (unchanged) and sub6_slotdisc_hfss_pred.beta_dir (unchanged).
    python sub6_island_hfss_pred.py run [--workers N] | table
"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_island_cell as IC
import sub6_slotdisc_hfss_pred as HP
from src.frequency import free_space_wavenumber

OUT = ROOT / "results" / "sub6_island_hfss" / "hfss_pred"
ROW_ID = "ISL_W2.9_q0.8_Li0.5_c0.6"
FREQS = [5.2, 5.5, 5.6, 5.7, 5.8, 6.1, 6.4]  # HP.FREQS plus 5.6 and 5.7: the ON state across y has a resonance between 5.5 and 5.8 GHz
TOL, MAXNODES, ROUNDS = HP.TOL, HP.MAXNODES, HP.ROUNDS
P = IC.P
CLASS_STATES = {"ISL": ["OFF_BAR50", "ON_BAR50"], "REF": ["NOCUT"]}


def row_cells():
    J = json.load(open(ROOT / "results" / "sub6_island_cell" / "island_rows.json"))
    r = [x for x in J["rows"] if x["row"]["id"] == ROW_ID][0]["row"]
    return r, IC.ref_cells(r["W"], r["ratio"])


def key(f, cls, d):
    return f"{f}|{cls}|{d}"


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
    row, rcells = row_cells()
    cells = {"ISL": row["cells"], "REF": rcells}
    nodes = load()
    keys = [(f, cls, d) for f in FREQS for cls in ("ISL", "REF") for d in ("x", "y")]
    t00 = time.time()
    with Pool(workers) as pool:
        for rnd in range(ROUNDS + 1):
            todo = []
            for (f, cls, d) in keys:
                k = key(f, cls, d)
                nd = nodes.setdefault(k, {})
                if not nd:
                    new = [1.10]
                else:
                    bs = []
                    for st in CLASS_STATES[cls]:
                        b, _, dist = HP.beta_dir(nd, st, d, f)
                        if dist > TOL:
                            bs.append(b)
                    new = []
                    for b in sorted(bs):
                        if all(abs(b - k2) > 3e-4 for k2 in list(nd) + new) and len(nd) + len(new) < MAXNODES:
                            new.append(b)
                todo += [((cls, IC.R_ANCHOR, cells[cls], d, kt, f * 1e9), k) for kt in new]
            if not todo:
                print("round", rnd, "no new nodes: converged", flush=True)
                break
            print("round", rnd, "tasks", len(todo), f"t={time.time() - t00:.0f}s", flush=True)
            res = pool.map(IC.node_task, [a for a, _ in todo], chunksize=1)
            for (a, k), r in zip(todo, res):
                nodes[k][a[4]] = r["Y"]
            save(nodes)
    save(nodes)
    print("done", f"{time.time() - t00:.0f}s", flush=True)


def table():
    nodes = load()
    out = {}
    for cls, sts in (("ISL", ("OFF_BAR50", "ON_BAR50")), ("REF", ("NOCUT",))):
        for st in sts:
            for d in ("x", "y"):
                rows = []
                for f in FREQS:
                    b, y, dist = HP.beta_dir(nodes[key(f, cls, d)], st, d, f)
                    rows.append(dict(f_GHz=f, beta_over_k0=b, phase_deg=float(np.degrees(b * free_space_wavenumber(f * 1e9) * P * 1e-3)), Y_dd_mS=[y.real * 1e3, y.imag * 1e3],
                                     nearest_node_dist=dist, n_nodes=len(nodes[key(f, cls, d)])))
                out[f"{cls}|{st}|{d}"] = rows
    (OUT / "pred_table.json").write_text(json.dumps(out, indent=1))
    for k, rows in out.items():
        print(k, [(r["f_GHz"], round(r["beta_over_k0"], 3), round(r["phase_deg"], 1), round(r["nearest_node_dist"], 4)) for r in rows])


if __name__ == "__main__":
    if sys.argv[1] == "run":
        w = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 64
        run(w)
    else:
        table()
