"""SUB6 slot-disc cell scan grid (MODEL_ONLY).  jobs/sub6-slotdisc-cell-scan.md (+ Amendment 2026-11-07).

    python sub6_slotdisc_grid.py list    [--substrates RT5880,RO3006,TMM6]
    python sub6_slotdisc_grid.py pilot   [--substrates ...]          one row, assemblies in parallel, prints time + projection
    python sub6_slotdisc_grid.py run     [--substrates ...] [--workers N]    resumable; rows -> results/sub6_slotdisc/rows_<SUB>.jsonl
    python sub6_slotdisc_grid.py mesh    [--substrates ...]          mesh 64 / 96 check on the 3 best-leakage rows

Substrates (eps_r, tan d, thickness mm: p list): RT5880 (2.2, 0.0009, 3.175: 11.6); RO3006 (6.15, 0.0020 = repo LAMINATES, 1.28 = thickest
standard single sheet, about half of the 2.5 mm asked: 10, 11.6); TMM6 (6.3, 0.0023, 2.54: 10).  Gap p - 2r in {0.4,0.6,0.8,1.0}, L/(2r) in
{0.5..0.9} with rim r - L/2 >= 0.4 mm, W in {0.4,0.8,1.2,1.6}.  mesh 64, F = 2n.  Single-thread BLAS.
One Z assembly per (geometry class NOCUT / OPEN / BRIDGE, direction x / y, kt node); the four loaded states share the BRIDGE assembly
(load matrix added); nodes are added adaptively until every state's beta is within 0.004 k0 of a node (max 4), Y(kt) linear in between.
Resonance flag: sign of Im Y_xx / Im Y_yy at normal incidence at 5.0 and 6.6 GHz (NOCUT and the BRIDGE states) differs -> RESONANCE_5.0_6.6.
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import platform
import subprocess
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_lib as S

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sub6_slotdisc"
N, FO = 64, 128
SUBS = {
    "RT5880": (("rt5880", 2.2, 0.0009, 3.175), [11.6]),
    "RO3006": (("ro3006", 6.15, 0.0020, 1.28), [10.0, 11.6]),
    "TMM6": (("tmm6", 6.3, 0.0023, 2.54), [10.0]),
}
GUESS = {"RT5880": 1.06, "RO3006": 1.3, "TMM6": 1.3}
GAPS = (0.4, 0.6, 0.8, 1.0)
RATIOS = (0.5, 0.6, 0.7, 0.8, 0.9)
WS = (0.4, 0.8, 1.2, 1.6)
ALL_STATES = [s for v in S.CLASS_STATES.values() for s in v]


def rows(substrates):
    out = []
    for name in substrates:
        sub, ps = SUBS[name]
        for p in ps:
            for gap in GAPS:
                r = (p - gap) / 2
                for q in RATIOS:
                    L = 2 * r * q
                    if r - L / 2 < 0.4 - 1e-9:
                        continue
                    for w in WS:
                        out.append(dict(id=f"{name}_p{p:g}_g{gap:g}_q{q:g}_W{w:g}", sub_name=name, sub=sub, p=p, gap=gap, ratio=q, r=r, L=L, W=w))
    return out


def jc(o):
    if isinstance(o, (complex, np.complexfloating)):
        return [float(o.real), float(o.imag)]
    if isinstance(o, np.ndarray):
        return jc(o.tolist())
    if isinstance(o, dict):
        return {str(k): jc(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jc(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    return o


def to_c(o):
    """inverse of jc for a 2x2 Y: [[ [re,im],[re,im] ],[...]]"""
    a = np.array(o, dtype=float)
    return a[..., 0] + 1j * a[..., 1]


def load_nodes(path):
    if not path.exists():
        return {}
    d = json.load(open(path))
    return {float(kt): {dr: {st: to_c(y) for st, y in sts.items()} for dr, sts in v.items()} for kt, v in d.items()}


def save_nodes(path, nodes):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({str(k): jc(v) for k, v in nodes.items()}))


def tasks_for(spec_t, kt, f=S.F):
    return [(spec_t, cls, d, kt, N, FO, f) for cls in S.CLASS_STATES for d in "xy"]


def run_row(row, pool=None):
    """Full row; returns result dict.  Resumes from results/sub6_slotdisc/ckpt/<id>.json."""
    t0 = time.time()
    spec_t = (row["sub"], row["p"], row["r"], row["L"], row["W"])
    ck = OUT / "ckpt" / f"{row['id']}.json"
    nodes = load_nodes(ck)
    mapper = pool.map if pool else map
    kt = None if nodes else GUESS[row["sub_name"]]
    if nodes:
        kt = S.next_node(nodes, ALL_STATES, row["sub"])
    for rnd in range(4):
        if kt is None:
            break
        outs = list(mapper(S.node_task, tasks_for(spec_t, kt)))
        for o in outs:
            nodes.setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
        save_nodes(ck, nodes)
        kt = S.next_node(nodes, ALL_STATES, row["sub"])
    # resonance proxy
    res = {}
    for fr in (5.0e9, 6.6e9):
        outs = list(mapper(S.node_task, [(spec_t, cls, "x", 0.0, N, FO, fr) for cls in ("NOCUT", "BRIDGE")]))
        res[fr] = {st: y for o in outs for st, y in o["Y"].items()}
    flip = []
    for st in res[5.0e9]:
        for i in (0, 1):
            if np.sign(res[5.0e9][st][i, i].imag) != np.sign(res[6.6e9][st][i, i].imag):
                flip.append(f"{st}:{'xy'[i]}")
    states = {st: S.state_summary(nodes, st, row["sub"], guess=GUESS[row["sub_name"]]) for st in ALL_STATES}
    for st, v in states.items():
        v["alpha_np_per_m"] = float(v["alpha_np_per_m"])
    d = lambda a, b: states[a]["Delta_mS"] / states[b]["Delta_mS"]
    out = dict(row=row, states=states, nodes=sorted(nodes), resonance_flip=flip, resonance_flag=bool(flip),
               rho_nom=d("ON_BAR50", "OFF_BAR50"), rho_max=d("ON_BAR50_MAX", "OFF_BAR50_HI"), rho_open=d("ON_BAR50", "OFF_OPEN"),
               rho_vs_nocut=d("ON_BAR50", "NOCUT"),
               dA_par=states["ON_BAR50"]["A_par_mS"] - states["OFF_BAR50"]["A_par_mS"],
               dA_perp=states["ON_BAR50"]["A_perp_mS"] - states["OFF_BAR50"]["A_perp_mS"],
               alpha_on=states["ON_BAR50"]["alpha_np_per_m"], alpha_off=states["OFF_BAR50"]["alpha_np_per_m"],
               seconds=time.time() - t0)
    a_on, a_off = out["alpha_on"], out["alpha_off"]
    out["alpha_ratio_on_over_off"] = a_on / a_off if a_off > 0 else float("nan")
    out["leakage"] = leakage_class(a_on, a_off)
    return out


def leakage_class(a_on, a_off, thr=1.0):
    """Per-row leakage class at 25 deg (user correction 2026-11-07): OK_BOTH | OK_ON_ONLY | OK_OFF_ONLY | LOW (threshold 1 Np/m)."""
    on, off = a_on >= thr, a_off >= thr
    return "OK_BOTH" if on and off else "OK_ON_ONLY" if on else "OK_OFF_ONLY" if off else "LOW"


def worker(row):
    return run_row(row, None)


def provenance(extra):
    def g(*a):
        try:
            return subprocess.check_output(["git", *a], cwd=ROOT, text=True).strip()
        except Exception:
            return "unknown"
    return dict(machine=platform.node(), os=platform.platform(), git_commit=g("rev-parse", "HEAD"), git_dirty=bool(g("status", "--porcelain", "--untracked-files=no")),
                python=sys.version, numpy=np.__version__, threads={k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
                labels="MODEL_ONLY NOT_HFSS_VALIDATED", mesh=N, floquet_order=FO, **extra)


def free_gb():
    try:
        out = subprocess.check_output(["powershell", "-c", "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB"], text=True)
        return float(out.strip())
    except Exception:
        return 32.0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["list", "pilot", "run", "mesh"])
    ap.add_argument("--substrates", default="RT5880,RO3006,TMM6")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    subs = [s.strip() for s in a.substrates.split(",") if s.strip()]
    bad = [s for s in subs if s not in SUBS]
    if bad:
        sys.exit(f"unknown substrate(s) {bad}; choose from {list(SUBS)}")
    rs = rows(subs)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "list":
        print(len(rs), "rows", {s: sum(r["sub_name"] == s for r in rs) for s in subs})
    elif a.stage == "pilot":
        row = [r for r in rs if abs(r["gap"] - 0.6) < 1e-9 and abs(r["ratio"] - 0.7) < 1e-9 and abs(r["W"] - 0.8) < 1e-9][0]
        t0 = time.time()
        with Pool(6) as pool:
            res = run_row(row, pool)
        dt = time.time() - t0
        n_assm = 6 * len(res["nodes"]) + 4
        print(f"pilot row {row['id']}: {dt:.0f} s wall with 6 parallel assemblies, {len(res['nodes'])} nodes, {n_assm} assemblies; "
              f"single-thread row = ~{n_assm * 330 / 3600:.1f} h-equivalent {n_assm * 330 / 60:.0f} min; projected {len(rs)} rows: {len(rs) * n_assm * 330 / 3600:.0f} core-h")
        (OUT / f"pilot_{row['id']}.json").write_text(json.dumps(jc(res), indent=1))
        print(json.dumps({k: jc(v) for k, v in res.items() if k in ("rho_nom", "rho_max", "rho_open", "alpha_on", "alpha_off", "leakage", "alpha_ratio_on_over_off", "resonance_flip", "nodes")}))
    elif a.stage == "run":
        done = set()
        for s in subs:
            f = OUT / f"rows_{s}.jsonl"
            if f.exists():
                done |= {json.loads(l)["row"]["id"] for l in f.read_text().splitlines() if l.strip()}
        todo = [r for r in rs if r["id"] not in done]
        nw = a.workers or int(min(64, max(4, free_gb() * 0.6 / 1.6)))
        prov = provenance(dict(stage="run", substrates=subs, workers=nw, rows_total=len(rs), rows_todo=len(todo), start=time.strftime("%Y-%m-%d %H:%M:%S")))
        (OUT / f"provenance_run_{'_'.join(subs)}.json").write_text(json.dumps(prov, indent=1, default=str))
        print(f"run {subs}: {len(todo)}/{len(rs)} rows, {nw} workers", flush=True)
        t0 = time.time()
        with Pool(min(nw, max(1, len(todo)))) as pool:
            for k, res in enumerate(pool.imap_unordered(worker, todo)):
                with open(OUT / f"rows_{res['row']['sub_name']}.jsonl", "a") as fh:
                    fh.write(json.dumps(jc(res)) + "\n")
                print(f"[{time.time() - t0:7.0f}s] {len(done) + k + 1}/{len(rs)} {res['row']['id']} alpha ON/OFF {res['alpha_on']:.2f}/{res['alpha_off']:.2f} {res['leakage']}", flush=True)
    else:
        allr = []
        for s in subs:
            f = OUT / f"rows_{s}.jsonl"
            if f.exists():
                allr += [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
        best = sorted(allr, key=lambda r: -min(r["alpha_on"], r["alpha_off"]))[:3]
        out = []
        for b in best:
            row = b["row"]
            spec_t = (tuple(row["sub"]), row["p"], row["r"], row["L"], row["W"])
            nodes = {}
            with Pool(6) as pool:
                for kt in b["nodes"]:
                    for o in pool.map(S.node_task, [(spec_t, cls, d, float(kt), 96, 192, S.F) for cls in S.CLASS_STATES for d in "xy"]):
                        nodes.setdefault(o["kt"], {}).setdefault(o["direction"], {}).update(o["Y"])
            st96 = {st: S.state_summary(nodes, st, row["sub"], guess=GUESS[row["sub_name"]]) for st in ("OFF_BAR50", "ON_BAR50")}
            out.append(dict(id=row["id"], mesh64={st: b["states"][st]["alpha_np_per_m"] for st in st96}, mesh96={st: st96[st]["alpha_np_per_m"] for st in st96}))
            print(out[-1], flush=True)
        (OUT / f"mesh_check_{'_'.join(subs)}.json").write_text(json.dumps(jc(out), indent=1))
