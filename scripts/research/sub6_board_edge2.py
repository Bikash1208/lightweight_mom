"""SUB6 board-edge amendment: EDGE_TO_EDGE variant (MODEL_ONLY).  Imports sub6_board_edge.py (unchanged) and adds the placements
  cx / cy / c : every switched line is a FULL CHORD of the hexagon (a complete lattice row along x, or column along y, from rim to rim, all its cells switched), so both string ends sit at the rim
                (bias by edge traces, no vias).  cx = x chords only, cy = y chords only, c = both orientations with no shared cell (a crossing would be a junction, not a series string).
Chords are the rows / columns whose lattice cells are contiguous from rim to rim; the row j = 0 and the column i = 0 are interrupted by the feed clearance (cell (0,0) removed) and are not chords.
Ring exclusion (excl = 1) here means: no chord within lattice rows / columns |index| <= 2 (the first 2 rings around the feed); the outer-ring exclusion cannot apply (the chord ends ARE the outer ring).
Selection: greedy on the ideal-cell projection objective of the free-segment variant ('b'), stopping when the total cell count is within the first chord that would overshoot the target by > 12 %;
N_sw = the actual cell count of the chosen chord set (the targets 135 / 240 are only targets).  Series string of n cells = n - 1 diodes x 0.95 V at 10 mA.
    python sub6_board_edge2.py run [--workers N]
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_board_edge as E
import sub6_board_island2 as X
import sub6_board_model4 as R4

B, M = E.B, E.M
RES = E.RES
LAUNCH = E.LAUNCH
TARGETS = (135, 240)
SIDES = (0.2, 0.22)


def chords(S, lk, side, excl, axes):
    poly = B.hexagon(side)
    n = lk.n
    out = []
    for ax in axes:
        for fixed in range(-n, n + 1):
            if excl and abs(fixed) <= 2:
                continue
            cells, interrupted = [], False
            for var in range(-n, n + 1):
                i, j = (var, fixed) if ax == "x" else (fixed, var)
                x, y = i * M.PITCH, j * M.PITCH
                inside = all((b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= -1e-9 for a, b in zip(poly, np.roll(poly, -1, axis=0)))
                if not inside:
                    continue
                idx = lk.t[i + n, j + n]
                if idx < 0:
                    interrupted = True
                else:
                    cells.append(int(idx))
            if cells and not interrupted and len(cells) >= 2:
                out.append(dict(axis=ax, fixed=fixed, cells=np.array(cells), L=len(cells)))
    return out


def place_chords(S, lk, target, mode, excl, rho_d, side):
    axes = {"cx": ("x",), "cy": ("y",), "c": ("x", "y")}[mode]
    ch = chords(S, lk, side, excl, axes)
    mi, dAi, dBi = B.march_for(side, 2.0)
    c = -rho_d / (1.0 - rho_d)
    f = target / S.n
    ax_, ay_ = B.unit_fields(dAi)
    bx_, by_ = B.unit_fields(dBi)
    amp = dAi["amp"]
    aA = S.cell_field(B.BEAM_A, amp * ax_ / (1 - f), amp * ay_ / (1 - f))
    aBst = S.cell_field(B.BEAM_B, amp * c * bx_ / (1 - f), amp * c * by_ / (1 - f))
    aBsw = S.cell_field(B.BEAM_B, amp * bx_ / ((1 - rho_d) * f), amp * by_ / ((1 - rho_d) * f))
    sums = np.array([aA.sum(), aBst.sum(), aBst.sum()])
    dvec = np.stack([-aA, aBsw - aBst, rho_d * aBsw - aBst], 1)
    dch = np.array([dvec[k["cells"]].sum(axis=0) for k in ch])
    used = np.zeros(S.n, bool)
    chosen, total = [], 0
    while True:
        ok = np.array([(total + k["L"] <= 1.12 * target) and not used[k["cells"]].any() for k in ch])
        if not ok.any():
            break
        sc = np.where(ok, E.obj_vec(sums[None, :] + dch), -np.inf)
        k = int(np.argmax(sc))
        chosen.append(k)
        used[ch[k]["cells"]] = True
        total += ch[k]["L"]
        sums = sums + dch[k]
        if total >= 0.92 * target:
            break
    mask = np.zeros(S.n, bool)
    lines = []
    for k in chosen:
        mask[ch[k]["cells"]] = True
        lines.append(dict(axis=ch[k]["axis"], cells=ch[k]["cells"], L=ch[k]["L"], fixed=ch[k]["fixed"]))
    return mask, lines


def design_edge2(side, n_sw, place, taper, excl, T):
    if place not in ("cx", "cy", "c"):
        return E.design_edge(side, n_sw, place, taper, excl, T)
    S = B.get_side(side)
    lk = M.Lookup(S)
    rho_d = float(np.mean(T["sw"]["ratio"]))
    mask, lines = place_chords(S, lk, n_sw, place, excl, rho_d, side)
    n_act = int(mask.sum())
    old = E.place_lines
    E.place_lines = lambda *a, **k: (mask, lines)
    try:
        D = E.design_edge(side, n_act, "b", taper, excl, T)
    finally:
        E.place_lines = old
    D["place"], D["n_target"] = place, n_sw
    return D


R4.design = design_edge2
E.DESIGN = design_edge2


def task(a):
    side, n, pl, taper, excl, launch = a
    rows = E.edge_task((side, n, pl, taper, excl, launch))
    for r in rows:
        r["n_target"] = n
    return rows


def run(workers):
    RES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    chord_tasks = [(s, n, m, t, ex, LAUNCH) for s in SIDES for n in TARGETS for m in ("cx", "cy", "c") for t in (0.5, 0.8, 1.0) for ex in (0, 1)]
    with Pool(min(workers, len(chord_tasks))) as pool:
        rows_c = [r for rs in pool.imap_unordered(task, chord_tasks, chunksize=1) for r in rs]
        # free-segment variant at the SAME actual N as each chord set
        combos = sorted({(r["side_mm"] / 1e3, int(r["N_sw_live"]), int(r["excl"])) for r in rows_c})
        free_tasks = [(s, n, pl, t, ex, LAUNCH) for (s, n, ex) in combos for pl in ("a", "b") for t in (0.5, 0.8, 1.0)]
        rows_f = [r for rs in pool.imap_unordered(task, free_tasks, chunksize=1) for r in rs]
    allrows = rows_c + rows_f
    allrows.sort(key=lambda r: (-r["side_mm"], r["N_sw"], r["placement"], r["taper"], r["excl"], r["variant"]))
    B.write_csv(RES / "points_edge_chords.csv", allrows)
    (RES / "provenance_run_chords.json").write_text(json.dumps(B.provenance(dict(stage="run_edge_chords", workers=workers, seconds=time.time() - t0, n_chord_designs=len(chord_tasks), n_free_designs=len(free_tasks),
                                                                      launch_fraction=LAUNCH, v_per_diode=E.V_DIODE)), indent=1, default=str), encoding="utf-8", newline="\n")
    print(f"done {time.time() - t0:.0f}s: chord rows {len(rows_c)}, free rows {len(rows_f)}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run"])
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    run(a.workers)
