"""Open-board gate test: pure-python spec (no pyaedt, never launches AEDT).

200 mm circular RT5880 board on a solid PEC ground, uniform lattice of static slotted discs, gate = SD_w2p9_ylines
cells in lattice columns i = 3, 4, 5 (x = 34.8, 46.4, 58.0 mm), centre probe feed as the 5.8 GHz final board.
"""
from __future__ import annotations

import math

import gate_spec as G

R_BOARD = 100.0
GATE_COLS = (3, 4, 5)
PIN_R = 0.635
HOLE_R = 2.05
PIN_ABOVE_GROUND_TOP = 12.1          # final-board probe: pin top 12.1 mm above the ground top, no hat
PAD_LAT, PAD_BOT, PAD_TOP = 15.0, 15.0, 40.0
PROBE_X, PROBE_HALF_Y, PROBE_Z_ABOVE_CU = 75.0, 40.0, 3.0
SOLVE_GHZ = 5.8
MAX_PASSES = 6
MAX_DELTA_S = 0.02
SWEEP = (5.6, 6.1, 0.1)
FREQS = [round(5.6 + 0.1 * k, 6) for k in range(6)]
RIM_MARGIN = 0.1                     # copper must stay this far inside the rim (line clipping)
FEED_CLEARANCE = 0.15
GROUND_BOTTOM = -G.RT_T - G.CU_T


def cells(r_board=R_BOARD, gate_cols=GATE_COLS, gate=True):
    """Lattice cells (square lattice, pitch 11.6, a cell on the origin) whose disc is fully inside the rim and clear of
    the probe opening. Each cell: i, j, x, y, gate (True only when `gate` and i in gate_cols)."""
    out = []
    nmax = int(r_board // G.PERIOD) + 1
    for i in range(-nmax, nmax + 1):
        for j in range(-nmax, nmax + 1):
            x, y = i * G.PERIOD, j * G.PERIOD
            if math.hypot(x, y) + G.DISC_R > r_board + 1e-9:
                continue
            if math.hypot(x, y) - G.DISC_R < HOLE_R + FEED_CLEARANCE:      # clashes with the feed opening
                continue
            out.append({"i": i, "j": j, "x": x, "y": y, "gate": bool(gate and i in gate_cols)})
    return out


def columns(cs):
    """{i: [cells]} for gate cells only, ordered by j."""
    cols = {}
    for c in cs:
        if c["gate"]:
            cols.setdefault(c["i"], []).append(c)
    return {i: sorted(v, key=lambda c: c["j"]) for i, v in cols.items()}


def contiguous(col):
    return all(b["j"] - a["j"] == 1 for a, b in zip(col, col[1:]))


def line_span(col, side, r_board=R_BOARD):
    """(xc_line_centre, y_lo, y_hi) of the continuous y-line of one column side (side -1 = LEFT net, +1 = RIGHT net).

    The line runs from the lower cell face to the upper cell face; each end is clipped so its outer corner stays
    RIM_MARGIN inside the rim."""
    xc = col[0]["x"]
    xl = xc + side * (G.SD_SPLIT_W / 2 + G.LINE_X_OFF)
    xo = abs(xl) + G.LINE_W / 2
    ymax_rim = math.sqrt(max((r_board - RIM_MARGIN) ** 2 - xo ** 2, 0.0))
    hi = min(col[-1]["y"] + G.PERIOD / 2, ymax_rim)
    lo = max(col[0]["y"] - G.PERIOD / 2, -ymax_rim)
    return xl, lo, hi


def expected(r_board=R_BOARD, gate_cols=GATE_COLS):
    ref = cells(r_board, gate_cols, gate=False)
    gt = cells(r_board, gate_cols, gate=True)
    n_gate = sum(c["gate"] for c in gt)
    return {"ref": {"discs": len(ref), "gate_cells": 0, "diode_sheets": 0, "column_nets": 0},
            "gate": {"discs": len(gt) - n_gate, "gate_cells": n_gate, "diode_sheets": n_gate,
                     "column_nets": 2 * len({c["i"] for c in gt if c["gate"]})}}


def dc_gate(gt):
    """Analytic DC graph: per gate column LEFT net and RIGHT net (halves + line); diode per gate cell joins them;
    no node touches ground; static discs isolated."""
    cols = columns(gt)
    plain, diodes = {}, []
    for i, col in cols.items():
        plain[f"L{i}"], plain[f"R{i}"] = set(), set()
        diodes += [(f"L{i}", f"R{i}")] * len(col)
    for c in gt:
        if not c["gate"]:
            plain[f"d{c['i']}_{c['j']}"] = set()
    with_d = {k: set(v) for k, v in plain.items()}
    for a, b in diodes:
        with_d[a].add(b)
        with_d[b].add(a)
    rep, ok = [], True
    for i in cols:
        e = {"column": i, "cells": len(cols[i]), "L_R_connected_without_diode": G._reach(plain, f"L{i}", f"R{i}"),
             "L_R_connected_with_diode": G._reach(with_d, f"L{i}", f"R{i}"), "contiguous": contiguous(cols[i])}
        e["ok"] = (not e["L_R_connected_without_diode"]) and e["L_R_connected_with_diode"] and e["contiguous"]
        ok &= e["ok"]
        rep.append(e)
    return {"ok": bool(ok), "columns": rep, "any_node_grounded": False}


def cases():
    rs0, xs0 = G.STATES["OFF"]
    base = {"kind": "BOARD", "lid_mm": 0.0, "island_l": 0.0, "n_gate": 0, "modes": 1, "renorm": True}
    return [dict(base, run_id="BOARD_REF", state="REF", gate=False, Rs=rs0, Xs=xs0),
            dict(base, run_id="BOARD_OFF", state="OFF", gate=True, Rs=G.STATES["OFF"][0], Xs=G.STATES["OFF"][1]),
            dict(base, run_id="BOARD_ON", state="ON", gate=True, Rs=G.STATES["ON"][0], Xs=G.STATES["ON"][1])]


def variables_for(row):
    r = row.get("r_board", R_BOARD)
    return {"period": (G.PERIOD, "mm"), "disc_r": (G.DISC_R, "mm"), "slot_l": (G.SLOT_L, "mm"), "slot_w": (G.SLOT_W, "mm"),
            "split_w": (G.SD_SPLIT_W, "mm"), "line_w": (G.LINE_W, "mm"), "diode_gap": (G.DIODE_GAP, "mm"),
            "diode_w": (G.DIODE_W, "mm"), "rt_t": (G.RT_T, "mm"), "cu_t": (G.CU_T, "mm"), "r_board": (r, "mm"),
            "pin_r": (PIN_R, "mm"), "hole_r": (HOLE_R, "mm"), "Rs": (row["Rs"], "ohm"), "Xs": (row["Xs"], "ohm")}
