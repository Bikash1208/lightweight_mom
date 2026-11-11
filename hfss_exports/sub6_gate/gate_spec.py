"""SUB6 surface-wave gate test: pure-python specification, DC-topology gate and case table.

This module must NEVER import pyaedt and NEVER launch AEDT; the tests import only this file.
Geometry is reused from Option 4 'G-ISL-noline' (opt1_geometry constants + opt4_cell.py build).
"""
from __future__ import annotations

import math

# ---- Option 4 / opt1 constants (mm) ------------------------------------------------------
PERIOD = 11.6
DISC_R = 5.435
SLOT_L = 8.696            # along y
SLOT_W = 2.9              # along x
ISLAND_W = 1.3
DIODE_GAP = 0.3           # sheet length along x (propagation axis)
DIODE_W = 0.6             # sheet width along y
VIA_R = 0.15              # via drill 0.3 mm diameter, solid PEC
PAD_R = 0.3               # top pad 0.6 mm diameter
RT_T = 3.175
CU_T = 0.035
ISLAND_LEFT = -SLOT_W / 2 + DIODE_GAP            # -1.15
ISLAND_CX = ISLAND_LEFT + ISLAND_W / 2           # -0.5  (via x)
ISLAND_LENGTHS = (4.348, 2.5)

# ---- strip model ------------------------------------------------------------------------
SLAB = 1.0                # bare slab between port and first cell edge
N_STATIC = 3
LID_MAIN = 20.0
LID_CHECK = 15.0
GATE_COUNTS = (1, 2, 3)
FREQ_START, FREQ_STOP, FREQ_STEP = 5.5, 6.1, 0.05
SOLVE_FREQ_GHZ = 5.8
N_FREQ = int(round((FREQ_STOP - FREQ_START) / FREQ_STEP)) + 1     # 13
FREQS_GHZ = [round(FREQ_START + i * FREQ_STEP, 6) for i in range(N_FREQ)]
MAX_PASSES = 20
MAX_DELTA_S = 0.01
MIN_PASSES = 2


def bar90_sheet(c_pf=0.15, rp=3000.0, ls_nh=0.4, wl=2.0, f_hz=5.8e9):
    """OFF sheet impedance (C || Rp) + Ls scaled by W/L; same formula as cell_param_pyaedt.bar90_sheet."""
    w = 2 * math.pi * f_hz
    z = 1 / complex(1 / rp, w * c_pf * 1e-12) + 1j * w * ls_nh * 1e-9
    return z.real * wl, z.imag * wl


_c13 = bar90_sheet(c_pf=0.13)
_c18 = bar90_sheet(c_pf=0.18)
STATES = {                                   # name: (Rs, Xs) ohm/sq
    "OFF": (22.23, -335.4),
    "ON": (1.60, 29.2),
    "OFF_C013": (29.56, -390.9),
    "OFF_C018": (15.45, -275.0),
}


# ---- layout -----------------------------------------------------------------------------
def layout(n_gate: int, period: float = PERIOD):
    """Cells along x (kind 'static'|'gate'), strip centred on x = 0; returns (cells, total_length)."""
    kinds = ["static"] * N_STATIC + ["gate"] * n_gate + ["static"] * N_STATIC
    cells_len = len(kinds) * period
    total = cells_len + 2 * SLAB
    x0 = -cells_len / 2
    cells = [{"index": i, "kind": k, "xc": x0 + (i + 0.5) * period} for i, k in enumerate(kinds)]
    return cells, total


def dc_graph(n_gate: int):
    """Wrapped DC graph WITHOUT diodes plus diode edges; nodes: disc_i, island_i (gate cells only), ground.

    y faces are a lattice pair, so every cell wraps onto itself in y (no new nodes); cells are separated in x by
    disc-to-disc air gaps (period - 2*disc_r > 0), so there are no x-edges between cells.
    """
    cells, _ = layout(n_gate)
    plain, diode = {"ground": set()}, []
    for c in cells:
        plain[f"disc_{c['index']}"] = set()
        if c["kind"] == "gate":
            i = f"island_{c['index']}"
            plain[i] = {"ground"}                     # grounded by its via
            plain["ground"].add(i)
            diode.append((f"disc_{c['index']}", i))   # BAR90 sheet bridges disc slot wall -> island wall
    return plain, diode


def _reach(adj, start, goal):
    seen, todo = set(), [start]
    while todo:
        n = todo.pop()
        if n == goal:
            return True
        if n not in seen:
            seen.add(n)
            todo.extend(adj[n])
    return False


def dc_gate(n_gate: int):
    """DC graph gate: island reaches ground via its via; disc reaches ground ONLY through the diode."""
    cells, _ = layout(n_gate)
    plain, diode = dc_graph(n_gate)
    with_diode = {k: set(v) for k, v in plain.items()}
    for a, b in diode:
        with_diode[a].add(b)
        with_diode[b].add(a)
    report, ok = [], True
    for c in cells:
        d = f"disc_{c['index']}"
        entry = {"cell": c["index"], "kind": c["kind"], "disc_reaches_ground_without_diode": _reach(plain, d, "ground"),
                 "disc_reaches_ground_with_diode": _reach(with_diode, d, "ground")}
        if c["kind"] == "gate":
            i = f"island_{c['index']}"
            entry["island_reaches_ground_via_via"] = _reach(plain, i, "ground")
            good = (entry["island_reaches_ground_via_via"] and not entry["disc_reaches_ground_without_diode"]
                    and entry["disc_reaches_ground_with_diode"])
        else:
            good = not entry["disc_reaches_ground_without_diode"] and not entry["disc_reaches_ground_with_diode"]
        entry["ok"] = bool(good)
        ok &= bool(good)
        report.append(entry)
    return {"ok": bool(ok), "cells": report}


def clearances(island_l: float):
    """Analytic copper clearances inside one gate cell (mm)."""
    left, right = ISLAND_LEFT, ISLAND_LEFT + ISLAND_W
    return {
        "island_to_slot_right": SLOT_W / 2 - right,
        "island_to_slot_end_y": (SLOT_L - island_l) / 2,
        "diode_gap_x": left - (-SLOT_W / 2),
        "via_pad_inside_island_x": min(ISLAND_CX - left, right - ISLAND_CX) - PAD_R,
        "via_pad_inside_island_y": island_l / 2 - PAD_R,
        "diode_sheet_inside_island_y": island_l / 2 - DIODE_W / 2,
        "slot_corner_inside_disc": DISC_R - math.hypot(SLOT_L / 2, SLOT_W / 2),
        "disc_to_disc_gap": PERIOD - 2 * DISC_R,
    }


# ---- case table -------------------------------------------------------------------------
def _fmt_l(l):
    return f"{l:.3f}".replace(".", "p")


def cases():
    rows = []

    def add(kind, n, lid, state, island_l):
        if kind == "REF":
            run_id = f"REF_H{int(lid)}"
            rs, xs = STATES["OFF"]          # unused (no gate cells); kept so variables always exist
            island_l = ISLAND_LENGTHS[0]
        else:
            run_id = f"G{_fmt_l(island_l)}_N{n}_{state}_H{int(lid)}"
            rs, xs = STATES[state]
        rows.append({"run_id": run_id, "kind": kind, "n_gate": 0 if kind == "REF" else n, "lid_mm": lid,
                     "state": None if kind == "REF" else state, "island_l": island_l, "Rs": rs, "Xs": xs})

    add("REF", 0, LID_MAIN, None, None)
    for l in ISLAND_LENGTHS:
        for n in GATE_COUNTS:
            for st in ("OFF", "ON"):
                add("GATE", n, LID_MAIN, st, l)
        for st in ("OFF_C013", "OFF_C018"):
            add("GATE", 3, LID_MAIN, st, l)
    add("REF", 0, LID_CHECK, None, None)
    for l in ISLAND_LENGTHS:
        for st in ("OFF", "ON"):
            add("GATE", 3, LID_CHECK, st, l)
    ids = [r["run_id"] for r in rows]
    assert len(ids) == len(set(ids))
    return rows


def row_layout(row):
    """(cells, total_length) for any row: REF/GATE use the gate layout; BENCH rows are n_static static discs."""
    if row["kind"] != "BENCH":
        return layout(row["n_gate"])
    n, slab = row["n_static"], row["slab_mm"]
    cells_len = n * PERIOD
    x0 = -cells_len / 2
    return ([{"index": i, "kind": "static", "xc": x0 + (i + 0.5) * PERIOD} for i in range(n)],
            cells_len + 2 * slab)


def bench_cases():
    """Bench diagnosis D1-D4 (no gate cells). REF_H20 (6 discs, 1 mm slab, H20, 1 mode) and REF_H15 are reused."""
    rows = []

    def add(run_id, n_static, slab, lid, modes):
        rs, xs = STATES["OFF"]
        rows.append({"run_id": run_id, "kind": "BENCH", "n_gate": 0, "lid_mm": float(lid), "state": None,
                     "island_l": ISLAND_LENGTHS[0], "Rs": rs, "Xs": xs, "n_static": n_static, "slab_mm": float(slab),
                     "modes": modes})

    add("BENCH_D1_slab_H20", 0, SLAB, 20, 1)
    for n in (2, 4, 8):
        add(f"BENCH_D2_n{n}_H20", n, SLAB, 20, 1)
    for h in (10, 30, 40):
        add(f"BENCH_D3_n6_H{h}", 6, SLAB, h, 3 if h >= 30 else 1)
    for sl in (5, 15):
        add(f"BENCH_D4_n6_slab{sl}_H20", 6, sl, 20, 1)
    return rows


def all_cases():
    return cases() + bench_cases()


def variables_for(row):
    """Every AEDT design variable (name -> (value, unit)) that is set and read back for one row."""
    return {"period": (PERIOD, "mm"), "disc_r": (DISC_R, "mm"), "slot_l": (SLOT_L, "mm"), "slot_w": (SLOT_W, "mm"),
            "island_l": (row["island_l"], "mm"), "island_w": (ISLAND_W, "mm"), "diode_gap": (DIODE_GAP, "mm"),
            "diode_w": (DIODE_W, "mm"), "via_r": (VIA_R, "mm"), "pad_r": (PAD_R, "mm"), "rt_t": (RT_T, "mm"),
            "cu_t": (CU_T, "mm"), "slab": (row.get("slab_mm", SLAB), "mm"), "lid_h": (row["lid_mm"], "mm"),
            "Rs": (row["Rs"], "ohm"), "Xs": (row["Xs"], "ohm")}
