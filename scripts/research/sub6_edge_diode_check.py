"""SUB6 edge-diode: resonance proxy + mesh 64/96 check (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED).

Rows: theta 90, L/2r 0.7 and 0.8, W 2.4 and 2.9 (4 rows), states OFF_BAR50 / ON_BAR50 (NO_DIODE solved with them: same assembly).
Mesh 64 (F 128) reuses the sub6_edge_diode ckpt; mesh 96 (F 192, h = 0.121 mm: sheet 8.3 x 2.9 cells -> rounded by the load matrix) has its own ckpt and
adaptive kt nodes (as the mesh-64 run).  Resonance proxy: sign of Im Y_uu / Y_vv (sheet-state tensors, kt = 0, slot axis) at 5.0 and 6.6 GHz vs 5.8 GHz.

    python sub6_edge_diode_check.py run96 | res | collect     (WORKERS env, default 24: mesh 96 is ~5 GB per assembly)
"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_edge_diode_cell as E
import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S

OUT = G.ROOT / "results/sub6_edge_diode_check"
WORKERS = int(os.environ.get("WORKERS", "24"))
IDS = [f"EDGE_th90_q{q}_W{w}" for q in (0.7, 0.8) for w in (2.4, 2.9)]
ROWS = [r for r in E.rows() if r["id"] in IDS]
STATES = E.ALL_STATES


def ck96(rid):
    return OUT / "ckpt96" / f"{rid}.json"


def run96():
    (OUT / "ckpt96").mkdir(parents=True, exist_ok=True)
    st = {}
    for r in ROWS:
        nodes = G.load_nodes(ck96(r["id"]))
        nk = S.next_node(nodes, STATES, E.SUB) if nodes else E.GUESS
        st[r["id"]] = dict(nodes=nodes, pending=[] if nk is None else [nk])
    t0 = time.time()
    with Pool(WORKERS, maxtasksperchild=1) as pool:
        for rnd in range(5):
            tasks = [(r, c, d, float(kt), 96, 192, S.F) for r in ROWS for kt in st[r["id"]]["pending"] for c in E.CLASS_STATES for d in "xy"]
            if not tasks:
                break
            print(f"[{time.time() - t0:6.0f}s] round {rnd}: {len(tasks)} assemblies, {WORKERS} workers", flush=True)
            for k, o in enumerate(pool.imap_unordered(E.node_task, tasks)):
                s = st[o["row"]]
                s["nodes"].setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
                G.save_nodes(ck96(o["row"]), s["nodes"])
                print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(tasks)} {o['row']} {o['cls']} {o['direction']} kt={o['kt']:.4f} {o['seconds']:.0f}s", flush=True)
            for r in ROWS:
                nk = S.next_node(st[r["id"]]["nodes"], STATES, E.SUB)
                st[r["id"]]["pending"] = [] if nk is None else [nk]
    print("RUN96_DONE", flush=True)


def res():
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "resonance.json"
    done = json.load(open(p)) if p.exists() else {}
    tasks = [(r, c, "x", 0.0, n, 2 * n, fr) for n in (64, 96) for r in ROWS for c in E.CLASS_STATES for fr in (5.0e9, S.F, 6.6e9) if f"{n}|{r['id']}|{c}|{fr:g}" not in done]
    meshes = {id(t): t[4] for t in tasks}
    with Pool(WORKERS, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(E.node_task, tasks):
            done[f"{o['unknowns'] and (64 if o['unknowns'] < 12000 else 96)}|{o['row']}|{o['cls']}|{o['f']:g}"] = G.jc(o["Y"])
            p.write_text(json.dumps(done))
            print("res", o["row"], o["cls"], o["f"], o["unknowns"], flush=True)
    print("RES_DONE", flush=True)


def collect():
    resj = json.load(open(OUT / "resonance.json")) if (OUT / "resonance.json").exists() else {}
    out = []
    for r in ROWS:
        rec = dict(id=r["id"])
        for n, path in ((64, E.ckpt(r["id"])), (96, ck96(r["id"]))):
            nodes = G.load_nodes(path)
            if not nodes:
                continue
            sm = {s: S.state_summary(nodes, s, E.SUB, guess=E.GUESS) for s in STATES}
            rec[f"mesh{n}"] = {s: dict(alpha=sm[s]["alpha_np_per_m"], beta=sm[s]["beta_over_k0"], Delta_mS=sm[s]["Delta_mS"]) for s in STATES}
            rec[f"mesh{n}"]["LEAKAGE"] = G.leakage_class(sm["ON_BAR50"]["alpha_np_per_m"], sm["OFF_BAR50"]["alpha_np_per_m"])
            rec[f"mesh{n}"]["nodes"] = sorted(nodes)
            fl = {}
            for cls, sts in E.CLASS_STATES.items():
                for fr in ("5e+09", "5.8e+09", "6.6e+09"):
                    k = f"{n}|{r['id']}|{cls}|{fr}"
                    if k in resj:
                        for s, v in resj[k].items():
                            y = G.to_c(v)
                            fl.setdefault(s, {})[fr] = [float(y[0, 0].imag), float(y[1, 1].imag)]
            if fl:
                rec[f"mesh{n}"]["ImY_uu_vv_S"] = fl
                rec[f"mesh{n}"]["resonance_flip"] = [f"{s}:{'uv'[i]}" for s, d in fl.items() if "5e+09" in d and "6.6e+09" in d for i in (0, 1) if np.sign(d["5e+09"][i]) != np.sign(d["6.6e+09"][i])]
        out.append(rec)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "edge_diode_check.json").write_text(json.dumps(G.jc(dict(provenance=G.provenance(dict(stage="edge_diode_check", r_mm=E.R)), results=out)), indent=1, default=str))
    for rec in out:
        print(rec["id"])
        for n in (64, 96):
            m = rec.get(f"mesh{n}")
            if m:
                print(f"  mesh{n}: " + " ".join(f"{s}[a={m[s]['alpha']:.3g} b={m[s]['beta']:.4f} D={abs(m[s]['Delta_mS']):.3g}]" for s in STATES), m["LEAKAGE"], "flip", m.get("resonance_flip"))


if __name__ == "__main__":
    {"run96": run96, "res": res, "collect": collect}[sys.argv[1]]()
