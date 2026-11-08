"""SUB6 slot-disc 2-GAP grid (MODEL_ONLY), job sub6-slotdisc-2gap-grid.  RT5880 3.175 mm, p 11.6, 5.8 GHz, F = 2n.
Cell = SlotDiscNotch (+u notch 1.00 mm, BAR50 sheet) + OPEN second notch g2 at the -u rim (sub6_slotdisc2_lib).  States NOCUT2 (both open), OFF_BAR50, ON_BAR50.
Mesh per g2 (odd n, 1-cell gap): g2 0.12 and 0.36 -> n 97 (1 and 3 cell rows).

    python sub6_slotdisc2_grid.py list  [--g2 0.12,0.36]
    python sub6_slotdisc2_grid.py pilot [--g2 0.12]          anchor-nearest row
    python sub6_slotdisc2_grid.py run   [--g2 0.12,0.36] [--workers N]    resumable; rows -> results/sub6_slotdisc2gap/rows_g2_<g2>.jsonl
"""
import argparse
import json
import math
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc2_lib as T
import sub6_slotdisc_grid as G

OUT = G.ROOT / "results/sub6_slotdisc2gap"
SUB = ("rt5880", 2.2, 0.0009, 3.175)
P, GUESS = 11.6, 1.09
ALL_STATES = [s for v in T.CLASS_STATES.values() for s in v]


def rows(g2s):
    out = []
    for g2 in g2s:
        for gap in G.GAPS:
            r = (P - gap) / 2
            for q in G.RATIOS:
                L = 2 * r * q
                if r - L / 2 < 0.35 - 1e-9:
                    continue
                for w in G.WS:
                    if math.hypot(L / 2, w / 2) >= r - 1e-9:
                        continue
                    out.append(dict(id=f"RT5880_p{P:g}_g2_{g2:g}_g{gap:g}_q{q:g}_W{w:g}", sub_name="RT5880", sub=SUB, p=P, gap=gap, ratio=q, r=r, L=L, W=w, g2=g2, mesh=T.MESH_FOR_G2[g2], g2_modelled_mm=T.gap_cells(P, r, L, w, g2, T.MESH_FOR_G2[g2])[1]))
    out.sort(key=lambda d: (abs(d["gap"] - 0.73) / 0.3) ** 2 + (abs(d["ratio"] - 0.8) / 0.1) ** 2 + (abs(d["W"] - 2.9) / 0.8) ** 2)
    return out


def run_row(row, pool=None):
    t0 = time.time()
    n, fo = row["mesh"], 2 * row["mesh"]
    spec_t = (row["sub"], row["p"], row["r"], row["L"], row["W"], row["g2"])
    ck = OUT / "ckpt" / f"{row['id']}.json"
    nodes = G.load_nodes(ck)
    mapper = pool.map if pool else map
    tasks = lambda kt, f=T.F: [(spec_t, cls, d, kt, n, fo, f) for cls in T.CLASS_STATES for d in "xy"]
    kt = T.next_node(nodes, ALL_STATES, row["sub"]) if nodes else GUESS
    for _ in range(4):
        if kt is None:
            break
        for o in mapper(T.node_task, tasks(kt)):
            nodes.setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
        G.save_nodes(ck, nodes)
        kt = T.next_node(nodes, ALL_STATES, row["sub"])
    res = {}
    for fr in (5.0e9, 6.6e9):
        outs = list(mapper(T.node_task, [(spec_t, cls, "x", 0.0, n, fo, fr) for cls in T.CLASS_STATES]))
        res[fr] = {st: y for o in outs for st, y in o["Y"].items()}
    flip = [f"{st}:{'xy'[i]}" for st in res[5.0e9] for i in (0, 1) if np.sign(res[5.0e9][st][i, i].imag) != np.sign(res[6.6e9][st][i, i].imag)]
    states = {st: T.state_summary(nodes, st, row["sub"], guess=GUESS) for st in ALL_STATES}
    for v in states.values():
        v["alpha_np_per_m"] = float(v["alpha_np_per_m"])
    d = lambda a, b: states[a]["Delta_mS"] / states[b]["Delta_mS"]
    out = dict(row=row, states=states, nodes=sorted(nodes), resonance_flip=flip, resonance_flag=bool(flip),
               rho_nom=d("ON_BAR50", "OFF_BAR50"), rho_open=d("ON_BAR50", "NOCUT2"),
               dA_perp=states["ON_BAR50"]["A_perp_mS"] - states["OFF_BAR50"]["A_perp_mS"],
               dA_par=states["ON_BAR50"]["A_par_mS"] - states["OFF_BAR50"]["A_par_mS"],
               alpha_on=states["ON_BAR50"]["alpha_np_per_m"], alpha_off=states["OFF_BAR50"]["alpha_np_per_m"], seconds=time.time() - t0)
    a_on, a_off = out["alpha_on"], out["alpha_off"]
    out["alpha_ratio_on_over_off"] = a_on / a_off if a_off > 0 else float("nan")
    out["leakage"] = G.leakage_class(a_on, a_off)
    return out


def worker(row):
    return run_row(row, None)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["list", "pilot", "run"])
    ap.add_argument("--g2", default="0.2")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    g2s = [float(x) for x in a.g2.split(",")]
    rs = rows(g2s)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "list":
        print(len(rs), "rows", {g: sum(r["g2"] == g for r in rs) for g in g2s})
    elif a.stage == "pilot":
        row = rs[0]
        t0 = time.time()
        with Pool(4) as pool:
            res = run_row(row, pool)
        print(f"pilot {row['id']} mesh {row['mesh']}: {time.time() - t0:.0f} s wall (4 parallel assemblies), nodes {res['nodes']}")
        (OUT / f"pilot_{row['id']}.json").write_text(json.dumps(G.jc(res), indent=1))
        print(json.dumps({k: G.jc(v) for k, v in res.items() if k in ("rho_nom", "rho_open", "alpha_on", "alpha_off", "alpha_ratio_on_over_off", "leakage", "resonance_flip")}))
    else:
        done = set()
        for g in g2s:
            f = OUT / f"rows_g2_{g:g}.jsonl"
            if f.exists():
                done |= {json.loads(l)["row"]["id"] for l in f.read_text().splitlines() if l.strip()}
        todo = [r for r in rs if r["id"] not in done]
        nw = a.workers or 48
        prov = G.provenance(dict(stage="run", g2=g2s, workers=nw, rows_total=len(rs), rows_todo=len(todo), start=time.strftime("%Y-%m-%d %H:%M:%S")))
        prov["mesh"] = {str(g): T.MESH_FOR_G2[g] for g in g2s}
        (OUT / f"provenance_run_g2_{'_'.join(f'{g:g}' for g in g2s)}.json").write_text(json.dumps(prov, indent=1, default=str))
        print(f"run g2={g2s}: {len(todo)}/{len(rs)} rows, {nw} workers", flush=True)
        t0 = time.time()
        with Pool(min(nw, max(1, len(todo))), maxtasksperchild=1) as pool:
            for k, res in enumerate(pool.imap_unordered(worker, todo)):
                with open(OUT / f"rows_g2_{res['row']['g2']:g}.jsonl", "a") as fh:
                    fh.write(json.dumps(G.jc(res)) + "\n")
                print(f"[{time.time() - t0:7.0f}s] {len(done) + k + 1}/{len(rs)} {res['row']['id']} alpha ON/OFF {res['alpha_on']:.2f}/{res['alpha_off']:.2f} {res['leakage']}", flush=True)
