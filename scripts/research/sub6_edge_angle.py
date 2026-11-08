"""SUB6 edge-angle scan (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED).  Job sub6-edge-angle.
Same cell and method as sub6_edge_diode_cell (session 1, commit 34c4a37, imported unchanged): closed slotted disc, r = 5.30 mm, BAR50 sheet 1.00 x 0.35 mm across the 1.00 mm
inter-cell gap to the +x neighbour (periodic line of diodes), mesh 64, F = 128, RT5880 3.175, 5.8 GHz.  Diode-to-slot angle theta = slot angle from the diode axis (+x).
theta 50,60,70,75,80,85,90,95,100,110,120,130 deg; L/2r 0.7, 0.8; W 2.4, 2.9 (48 rows).  States OFF_BAR50, ON_BAR50 (class BRIDGE only, 2 assemblies per node round).
Mirror check: theta and 180 - theta are mirror images about the diode axis (50/130, 60/120, 70/110, 80/100, 85/95): their alpha / beta / |Delta| must agree.

    python sub6_edge_angle.py run [--workers N] | collect
"""
import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_edge_diode_cell as E

G, S = E.G, E.S
OUT = G.ROOT / "results/sub6_edge_angle"
THETAS = (50, 60, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130)
QS, WS = (0.7, 0.8), (2.4, 2.9)
STATES = ["OFF_BAR50", "ON_BAR50"]


def rows():
    out = []
    for q in QS:
        for w in WS:
            for th in THETAS:
                out.append(dict(id=f"ANG_th{th}_q{q:g}_W{w:g}", theta=th, q=q, L=2 * E.R * q, W=w, r=E.R, sub=E.SUB))
    out.sort(key=lambda r: abs(r["theta"] - 90))  # nearest the broadside-like 90 deg first
    return out


def ckpt(rid):
    return OUT / "ckpt" / f"{rid}.json"


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    rs = rows()
    nodes = {r["id"]: G.load_nodes(ckpt(r["id"])) for r in rs}
    pend = {r["id"]: ([E.GUESS] if not nodes[r["id"]] else ([k] if (k := S.next_node(nodes[r["id"]], STATES, E.SUB)) is not None else [])) for r in rs}
    prov = G.provenance(dict(stage="edge_angle", workers=workers, rows=len(rs), start=time.strftime("%Y-%m-%d %H:%M:%S"), base_commit="34c4a37"))
    (OUT / "provenance_run.json").write_text(json.dumps(prov, indent=1, default=str))
    t0 = time.time()
    with Pool(workers, maxtasksperchild=1) as pool:
        for rnd in range(5):
            tasks = [(r, "BRIDGE", d, float(kt), E.N, E.FO, S.F) for r in rs for kt in pend[r["id"]] for d in "xy"]
            if not tasks:
                break
            print(f"[{time.time() - t0:6.0f}s] round {rnd}: {len(tasks)} assemblies, {workers} workers", flush=True)
            for k, o in enumerate(pool.imap_unordered(E.node_task, tasks)):
                nodes[o["row"]].setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
                G.save_nodes(ckpt(o["row"]), nodes[o["row"]])
                print(f"[{time.time() - t0:6.0f}s] {k + 1}/{len(tasks)} {o['row']} {o['direction']} kt={o['kt']:.4f} {o['seconds']:.0f}s", flush=True)
            pend = {r["id"]: ([k] if (k := S.next_node(nodes[r["id"]], STATES, E.SUB)) is not None else []) for r in rs}
    print("RUN_DONE", flush=True)


def collect():
    recs, lines = [], ["id,theta,q,W,alpha_OFF,alpha_ON,beta_OFF,beta_ON,absDelta_OFF_mS,absDelta_ON_mS,A_perp_OFF_mS,A_perp_ON_mS,Yuv_over_Delta_ON,LEAKAGE,no_root"]
    for r in rows():
        nodes = G.load_nodes(ckpt(r["id"]))
        if not nodes:
            continue
        sm = {st: S.state_summary(nodes, st, E.SUB, guess=E.GUESS) for st in STATES}
        nr = [st for st, v in sm.items() if v.get("status") == "NO_ROOT"]
        a_off, a_on = sm["OFF_BAR50"]["alpha_np_per_m"], sm["ON_BAR50"]["alpha_np_per_m"]
        yuv = float("nan")
        if not nr:
            yx, yy = S.interp_y(nodes, "ON_BAR50", sm["ON_BAR50"]["beta_over_k0"])
            yuv = float(max(abs(yx[0, 1]), abs(yy[1, 0])) * 1e3 / abs(sm["ON_BAR50"]["Delta_mS"]))
        cls = G.leakage_class(a_on, a_off)
        rec = dict(row=r, alpha_OFF=a_off, alpha_ON=a_on, states=sm, Yuv_over_Delta_ON=yuv, LEAKAGE=cls, no_root=nr, nodes=sorted(nodes))
        recs.append(rec)
        f = lambda v: f"{abs(v):.4g}"
        lines.append(f"{r['id']},{r['theta']},{r['q']:g},{r['W']:g},{a_off:.4g},{a_on:.4g},{sm['OFF_BAR50']['beta_over_k0']:.4g},{sm['ON_BAR50']['beta_over_k0']:.4g},"
                     f"{f(sm['OFF_BAR50']['Delta_mS'])},{f(sm['ON_BAR50']['Delta_mS'])},{sm['OFF_BAR50']['A_perp_mS']:.4g},{sm['ON_BAR50']['A_perp_mS']:.4g},{yuv:.3g},{cls},{'|'.join(nr) or '-'}")
    (OUT / "edge_angle.json").write_text(json.dumps(G.jc(dict(results=recs)), indent=1, default=str))
    (OUT / "edge_angle.csv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "collect"])
    ap.add_argument("--workers", type=int, default=48)
    a = ap.parse_args()
    run(a.workers) if a.stage == "run" else collect()
