"""PARALLEL TEST item 1 compute: full Bloch tensors at kt = beta*k0 (beta of the table, per state) for theta' = 0, 22.5, 45, 67.5, 90 deg (direction of kt relative to the slot axis x)
for the 10 static NOCUT rows used on the board (mesh 64, as the scan) and the island cell ISL_W2.9_q0.8_Li0.5_c0.6 (mesh 96 as the board table; OFF at beta_OFF, ON at beta_ON; plus a theta' = -45 symmetry check).
    python diag_full_tensor_run.py <workers>      (resumable: results/sub6_final_board/diag_full_tensor/tensors.json)
"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8"); os.environ.setdefault("OMP_NUM_THREADS", "8")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
import numpy as np

import diag_full_tensor_lib as L

OUT = ROOT / "results" / "sub6_final_board" / "diag_full_tensor" / "tensors.json"
ANG = (0.0, 22.5, 45.0, 67.5, 90.0)
STATIC_IDS = ["RT5880_p11.6_g1_q0.8_W3.6", "RT5880_p11.6_g0.6_q0.7_W2.9", "RT5880_p11.6_g0.4_q0.7_W3.6", "RT5880_p11.6_g0.73_q0.7_W3.6", "RT5880_p11.6_g0.4_q0.7_W2.9",
              "RT5880_p11.6_g0.6_q0.8_W1.6", "RT5880_p11.6_g0.6_q0.7_W3.6", "RT5880_p11.6_g1_q0.9_W1.6", "RT5880_p11.6_g1_q0.8_W2.9", "RT5880_p11.6_g0.4_q0.7_W2.4"]
ISL_ID = "ISL_W2.9_q0.8_Li0.5_c0.6"


def jobs():
    rows = {r["row"]["id"]: r for r in (json.loads(l) for l in open(ROOT / "results/sub6_slotdisc/rows_RT5880.jsonl") if l.strip())}
    j = json.load(open(ROOT / "results/sub6_final_board/mesh96/input_island_rows_n96_d789d9a.json"))
    isl = {r["row"]["id"]: r for r in j["rows"]}[ISL_ID]
    out = []
    for i in STATIC_IDS:
        r = rows[i]
        b = r["states"]["NOCUT"]["beta_over_k0"]
        for t in ANG:
            out.append((f"{i}|NOCUT|{t}", ("static", r, 64, t, b, ["NOCUT"])))
    boff, bon = isl["states"]["OFF_BAR50"]["beta_over_k0"], isl["states"]["ON_BAR50"]["beta_over_k0"]
    spec = dict(cells=isl["row"]["cells"])
    for t in ANG + (-45.0,):
        out.append((f"{ISL_ID}|OFF_BAR50|{t}", ("island", spec, 96, t, boff, ["OFF_BAR50"])))
        if t != -45.0:
            out.append((f"{ISL_ID}|ON_BAR50|{t}", ("island", spec, 96, t, bon, ["ON_BAR50"])))
    return out


def work(item):
    key, a = item
    r = L.tensor_task(a)
    st = a[5][0]
    Y = r["Y"][st]
    return key, dict(row=key.split("|")[0], state=st, theta=a[3], kt_over_k0=a[4], mesh=a[2], Y_re=Y.real.tolist(), Y_im=Y.imag.tolist(), seconds=r["seconds"], unknowns=r["unknowns"])


if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    done = json.load(open(OUT)) if OUT.exists() else {}
    todo = [it for it in jobs() if it[0] not in done]
    # long (mesh 96) jobs first
    todo.sort(key=lambda it: -it[1][2])
    print(len(todo), "assemblies,", nw, "workers", flush=True)
    t0 = time.time()
    with Pool(nw) as p:
        for key, res in p.imap_unordered(work, todo):
            done[key] = res
            OUT.write_text(json.dumps(done, indent=1), encoding="utf-8", newline="\n")
            print(f"{time.time() - t0:7.0f}s done {key} {res['seconds']:.0f}s ({len(done)} stored)", flush=True)
    print("ALL DONE", flush=True)
