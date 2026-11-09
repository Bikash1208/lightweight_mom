"""SUB6 FINAL BOARD, 13-pad variant: export + gates from the routed result of sub6_final_13pad.py (route13.pkl).  MODEL_ONLY / NOT_HFSS_VALIDATED.
    SUB6_MESH=96 SUB6_VARIANT=13pad python sub6_final_13pad_export.py
Writes fab/{traces.json, cell_table.csv, drill.csv, drill.drl, string_table.csv, layout_summary.json, sub6_final_board.dxf, connectivity.json, drc.json, rf_unchanged.json}.
Gates (independent of the router): CONNECTIVITY_OK/FAIL, DRC_OK/FAIL, RF_UNCHANGED/RF_CHANGED.
"""
import csv
import json
import math
import os
import pickle
import shutil
import sys
from pathlib import Path

os.environ["SUB6_VARIANT"] = "13pad"
os.environ.setdefault("SUB6_MESH", "96")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.spatial import cKDTree

import sub6_final_layout as L
import sub6_final_13pad as M

SRC = M.SRC
OUT = L.OUT
FAB = OUT / "fab"
R3 = lambda v: round(float(v) + 0.0, 3)
PAD_NAMES = {}


def pad_names(assign):
    """pads are named top to bottom: S01.. for the 12 drive pads, RETURN in the middle (k = 6)."""
    names, j = {}, 0
    for k in range(M.N_PAD):
        if k == M.RET_K:
            names[k] = "RETURN"
        else:
            j += 1
            names[k] = f"S{j:02d}"
    return names


def collect_traces(res):
    tr = []
    for n in res["links"]:
        tr.append(dict(kind="link", net=f"L{n['string']}_{n['pos']}", string=int(n["string"]), pos=int(n["pos"]), w=L.W_TR, poly=[tuple(p) for p in n["poly"]]))
    for sid, poly in sorted(res["leads"].items()):
        tr.append(dict(kind="lead", net=f"LEAD{sid}", string=int(sid), pos="start", w=M.W_LEAD, poly=[tuple(p) for p in poly]))
    for i, b in enumerate(res["branches"]):
        tr.append(dict(kind="return", net="RETURN", string=int(b["string"]), pos=f"{b['kind']}{i}", w=float(b["w"]), poly=[tuple(p) for p in b["poly"]]))
    return tr


def seg_dist(P, a, b):
    """distance from points P (n,2) to the segment a-b"""
    ab = np.array(b) - np.array(a)
    t = np.clip(((P - a) @ ab) / max(ab @ ab, 1e-18), 0, 1)
    return np.hypot(*(P - (a + t[:, None] * ab)).T)


def poly_dist(P, poly):
    d = np.full(len(P), 1e9)
    for a, b in zip(poly[:-1], poly[1:]):
        d = np.minimum(d, seg_dist(P, np.array(a), np.array(b)))
    return d


def checks(res, traces, names):
    """independent connectivity + DRC from the final geometry only."""
    vias = res["vias"]
    cell_of = {}
    vpts, vlab = [], []
    for ci, (vi, vd) in vias.items():
        vpts += [vi, vd]
        vlab += [(ci, "I"), (ci, "D")]
    vpts = np.array(vpts)
    pads = res["pads"]
    n_tr = len(traces)
    # union-find over: traces [0..n_tr), via pads [n_tr..n_tr+480), connector pads after
    nV = len(vpts)
    parent = list(range(n_tr + nV + len(pads)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    samples = [np.array(L.resample(t["poly"], 0.02)) for t in traces]
    vt = cKDTree(vpts)
    PADR = L.PAD_R
    touches = []
    for i, (t, s) in enumerate(zip(traces, samples)):
        idx = vt.query_ball_point(s, t["w"] / 2 + PADR + 1e-6)
        for v in {j for lst in idx for j in lst}:
            union(i, n_tr + v)
            touches.append((i, "via", v))
        for k, p in enumerate(pads):
            q = np.abs(s - p) - 0.5
            dd = np.hypot(np.maximum(q[:, 0], 0), np.maximum(q[:, 1], 0))
            if dd.min() <= t["w"] / 2 + 1e-6:
                union(i, n_tr + nV + k)
                touches.append((i, "pad", k))
    allp = np.vstack(samples)
    lab = np.concatenate([np.full(len(s), i) for i, s in enumerate(samples)])
    ttree = cKDTree(allp)
    wv = np.array([t["w"] for t in traces])
    pr = ttree.query_pairs(0.8, output_type="ndarray")
    pr = pr[lab[pr[:, 0]] != lab[pr[:, 1]]]
    dd = np.linalg.norm(allp[pr[:, 0]] - allp[pr[:, 1]], axis=1)
    ii, jj = lab[pr[:, 0]], lab[pr[:, 1]]
    gap = dd - wv[ii] / 2 - wv[jj] / 2
    for a, b in zip(ii[gap <= 1e-6], jj[gap <= 1e-6]):
        union(int(a), int(b))
    comps = {}
    for x in range(len(parent)):
        comps.setdefault(find(x), []).append(x)
    # --- expectations
    cells = res["d"]["cells"]
    strings = res["d"]["strings"]
    fails = []
    pos_of = {}
    for s in strings:
        for k, ci in enumerate(s["cells"]):
            pos_of[ci] = (s["id"], k + 1)

    def comp_of_via(ci, kind):
        j = vlab.index((ci, kind))
        return find(n_tr + j)

    def members(root):
        m = comps[root]
        return dict(traces=[traces[x]["net"] for x in m if x < n_tr], vias=[vlab[x - n_tr] for x in m if n_tr <= x < n_tr + nV], pads=[names[x - n_tr - nV] for x in m if x >= n_tr + nV])

    ret_root = find(n_tr + nV + M.RET_K)
    mem_ret = members(ret_root)
    exp_ret = sorted((s["cells"][-1], "D") for s in strings)
    if sorted(mem_ret["vias"]) != exp_ret:
        fails.append(f"RETURN component vias {sorted(mem_ret['vias'])} != 12 end disc vias")
    if mem_ret["pads"] != ["RETURN"]:
        fails.append(f"RETURN component pads {mem_ret['pads']}")
    # each drive pad -> chain
    seen_cells, chains = set(), {}
    for k in range(M.N_PAD):
        if k == M.RET_K:
            continue
        root = find(n_tr + nV + k)
        m = members(root)
        if m["pads"] != [names[k]]:
            fails.append(f"pad {names[k]} component contains pads {m['pads']}")
        isl = [v for v in m["vias"] if v[1] == "I"]
        dsc = [v for v in m["vias"] if v[1] == "D"]
        if len(m["vias"]) != 1 or len(isl) != 1:
            fails.append(f"pad {names[k]} component vias {m['vias']} (expected exactly one island via)")
            continue
        ci = isl[0][0]
        sid, pos = pos_of[ci]
        chain = []
        while True:
            if pos_of[ci][1] != len(chain) + 1:
                fails.append(f"pad {names[k]}: order broken at cell {ci}")
                break
            if ci in seen_cells:
                fails.append(f"cell {ci} visited twice")
                break
            seen_cells.add(ci)
            chain.append(ci)
            r = comp_of_via(ci, "D")  # diode: island -> disc of the same cell
            if r == ret_root:
                break
            m2 = members(r)
            nxt = [v for v in m2["vias"] if v[1] == "I"]
            exp_traces = [t for t in m2["traces"]]
            if len(m2["vias"]) != 2 or len(nxt) != 1 or m2["pads"] or len(exp_traces) < 1:
                fails.append(f"after cell {ci} (string {sid}) unexpected component {m2}")
                break
            ci = nxt[0][0]
            if pos_of[ci][0] != sid:
                fails.append(f"string mixed at cell {ci}")
                break
        chains[names[k]] = dict(string=int(sid), n_diodes=len(chain), first_pos=pos_of[chain[0]][1], last_pos=pos_of[chain[-1]][1], ends_on_RETURN=bool(comp_of_via(chain[-1], "D") == ret_root))
        if len(chain) != 20 or not chains[names[k]]["ends_on_RETURN"]:
            fails.append(f"pad {names[k]}: chain length {len(chain)}, ends on RETURN {chains[names[k]]['ends_on_RETURN']}")
    all_switched = {ci for ci in vias}
    floating = sorted(all_switched - seen_cells)
    if floating:
        fails.append(f"floating diodes (cells not on a pad-to-RETURN chain): {floating[:10]} ({len(floating)})")
    # vias outside any expected component, shorts: every component with >1 pad
    for root, m in comps.items():
        mm = members(root)
        if len(mm["pads"]) > 1:
            fails.append(f"short between pads {mm['pads']}")
    # every trace belongs to a component with a pad or a link component of exactly 2 vias
    conn = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED", state="CONNECTIVITY_OK" if not fails else "CONNECTIVITY_FAIL", fails=fails[:40], n_fails=len(fails), chains=chains, n_traces=n_tr, n_vias=nV,
                n_diodes_on_chains=len(seen_cells), n_components=len(comps), return_component_traces=len(mem_ret["traces"]), polarity="every diode: anode = island via (string side), cathode = disc via (next link or RETURN); chain traversal island(1)->disc(1)->link->island(2)... pad -> RETURN")
    # ---------------- DRC
    isl_net_ = lambda s_, p_: f"S{s_:02d}" if p_ == 1 else f"N{s_:02d}_{p_ - 1:02d}"
    dsc_net_ = lambda s_, p_: "RETURN" if p_ == 20 else f"N{s_:02d}_{p_:02d}"
    vnet = [(isl_net_(*pos_of[ci]) if kd == "I" else dsc_net_(*pos_of[ci])) for ci, kd in vlab]   # via nets (vias of the trace's own net are exempt from the foreign-via clearance)
    names_net = np.array([t["net"] for t in traces])
    net_id = {n: i for i, n in enumerate(sorted(set(names_net)))}
    nid = np.array([net_id[n] for n in names_net])
    dmin = dict(trace_to_trace_gap_mm=9.0, trace_edge_to_foreign_antipad_edge_mm=9.0, trace_edge_to_foreign_via_pad_edge_mm=9.0, trace_edge_to_foreign_connector_pad_edge_mm=9.0, trace_edge_to_board_outline_mm=9.0, trace_edge_to_probe_opening_edge_mm=9.0)
    drc_fail = []
    # trace-trace (different nets)
    m = nid[ii] != nid[jj]
    if m.any():
        dmin["trace_to_trace_gap_mm"] = float(gap[m].min())
        bad = np.nonzero(m & (gap < M.MIN_CLR - 1e-9 if hasattr(M, "MIN_CLR") else m & (gap < 0.15 - 1e-9)))[0]
        if len(bad):
            drc_fail.append(("trace_to_trace", sorted({(traces[int(ii[b])]["net"], traces[int(jj[b])]["net"]) for b in bad})[:10]))
    # trace vs foreign vias (anti-pad and pad), own = via centre within 0.05 of a trace end point
    for i, (t, s) in enumerate(zip(traces, samples)):
        ends = [np.array(t["poly"][0]), np.array(t["poly"][-1])]
        d, j = vt.query(s, k=3)
        for kk in range(3):
            okm = np.array([(min(np.linalg.norm(vpts[jj_] - e) for e in ends) > 0.05) and (vnet[jj_] != t['net']) for jj_ in j[:, kk]])
            if okm.any():
                dm = float(d[okm, kk].min())
                dmin["trace_edge_to_foreign_antipad_edge_mm"] = min(dmin["trace_edge_to_foreign_antipad_edge_mm"], dm - t["w"] / 2 - L.ANTIPAD_R)
                dmin["trace_edge_to_foreign_via_pad_edge_mm"] = min(dmin["trace_edge_to_foreign_via_pad_edge_mm"], dm - t["w"] / 2 - PADR)
        # connector pads (foreign: pad of another net)
        for k, p in enumerate(pads):
            q = np.abs(s - p) - 0.5
            dp = np.hypot(np.maximum(q[:, 0], 0), np.maximum(q[:, 1], 0))
            if min(np.linalg.norm(p - e) for e in ends) < 1e-6 or names.get(k) == t["net"]:
                continue
            dmin["trace_edge_to_foreign_connector_pad_edge_mm"] = min(dmin["trace_edge_to_foreign_connector_pad_edge_mm"], float(dp.min()) - t["w"] / 2)
        out = res["d"]["outline"]
        o = np.array(out)
        for a, b in zip(o, np.roll(o, -1, axis=0)):
            ln = math.hypot(*(b - a))
            cr = ((b[0] - a[0]) * (s[:, 1] - a[1]) - (b[1] - a[1]) * (s[:, 0] - a[0])) / ln
            dmin["trace_edge_to_board_outline_mm"] = min(dmin["trace_edge_to_board_outline_mm"], float(cr.min()) - t["w"] / 2)
        dmin["trace_edge_to_probe_opening_edge_mm"] = min(dmin["trace_edge_to_probe_opening_edge_mm"], float(np.hypot(s[:, 0], s[:, 1]).min()) - t["w"] / 2 - 2.05)
    # trace crossing an anti-pad = edge clearance to a foreign antipad < 0 (checked by the same quantity)
    req = dict(trace_to_trace_gap_mm=0.15, trace_edge_to_foreign_antipad_edge_mm=M.ANTIPAD_MIN, trace_edge_to_foreign_via_pad_edge_mm=0.15, trace_edge_to_foreign_connector_pad_edge_mm=0.15, trace_edge_to_board_outline_mm=0.15, trace_edge_to_probe_opening_edge_mm=0.15)
    viol = {k: round(v, 4) for k, v in dmin.items() if v < req[k] - 1e-9}
    ret_w = [t["w"] for t in traces if t["kind"] == "return"]
    drc = dict(label="MODEL_ONLY / NOT_HFSS_VALIDATED", state="DRC_OK" if not viol and not drc_fail else "DRC_FAIL", violations=viol, trace_pair_violations=drc_fail,
               min_clearances_mm={k: round(v, 4) for k, v in dmin.items()}, requirement_mm=req, min_trace_width_mm=round(min(t["w"] for t in traces), 3), min_RETURN_width_mm=round(min(ret_w), 3),
               RETURN_widths_mm={str(w): int(sum(1 for x in ret_w if abs(x - w) < 1e-9)) for w in sorted(set(ret_w))}, RETURN_length_by_width_mm={str(w): round(sum(L.poly_len(t["poly"]) for t in traces if t["kind"] == "return" and abs(t["w"] - w) < 1e-9), 1) for w in sorted(set(ret_w))},
               antipad_crossings=int(dmin["trace_edge_to_foreign_antipad_edge_mm"] < 0), n_traces=n_tr)
    return conn, drc


def rect_pts(cx, cy, hx, hy, th, ox=0.0, oy=0.0):
    return [(R3(cx + L.rot(np.array([ox + sx * hx, oy + sy * hy]), th)[0]), R3(cy + L.rot(np.array([ox + sx * hx, oy + sy * hy]), th)[1])) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]


def write_dxf(cells, traces, pads, names, hexp, vias):
    import ezdxf
    doc = ezdxf.new("R2010", setup=True)
    doc.units = 4
    msp = doc.modelspace()
    for name, col in (("TOP_DISC", 1), ("TOP_SLOT", 2), ("TOP_ISLAND", 3), ("DIODE_SHEET", 6), ("VIA_DRILL", 5), ("VIA_PAD", 4), ("ANTIPAD", 8), ("BOT_TRACE", 30), ("BOT_RETURN", 31), ("BOT_PAD", 40), ("OUTLINE", 7), ("SMA_HOLE", 7), ("CONNECTOR", 50)):
        doc.layers.add(name, color=col)
    for c in cells:
        th = math.radians(float(c["slot_axis_deg"]))
        x, y = float(c["x_mm"]), float(c["y_mm"])
        msp.add_circle((x, y), float(c["disc_r_mm"]), dxfattribs={"layer": "TOP_DISC"})
        msp.add_lwpolyline(rect_pts(x, y, float(c["slot_L_mm"]) / 2, float(c["slot_W_mm"]) / 2, th), close=True, dxfattribs={"layer": "TOP_SLOT"})
        if c["class"] == "SWITCHED":
            msp.add_lwpolyline(rect_pts(x, y, 2.174, 0.65, th, 0.0, -0.20), close=True, dxfattribs={"layer": "TOP_ISLAND"})
            msp.add_lwpolyline(rect_pts(x, y, 0.175, 0.5, th, 0.0, 0.95), close=True, dxfattribs={"layer": "DIODE_SHEET"})
            for p in ((float(c["island_via_x"]), float(c["island_via_y"])), (float(c["disc_via_x"]), float(c["disc_via_y"]))):
                for lay, r in (("VIA_DRILL", 0.15), ("VIA_PAD", 0.30), ("ANTIPAD", 0.50), ("BOT_PAD", 0.30)):
                    msp.add_circle(p, r, dxfattribs={"layer": lay})
    for t in traces:
        msp.add_lwpolyline([(R3(a), R3(b)) for a, b in t["poly"]], dxfattribs={"layer": "BOT_RETURN" if t["kind"] == "return" else "BOT_TRACE", "const_width": t["w"]})
    for k, p in enumerate(pads):
        msp.add_lwpolyline(rect_pts(R3(p[0]), R3(p[1]), 0.5, 0.5, 0.0), close=True, dxfattribs={"layer": "CONNECTOR"})
        msp.add_text(names[k], dxfattribs={"layer": "CONNECTOR", "height": 0.8, "insert": (R3(p[0] + 1.2), R3(p[1] - 0.4))})
    msp.add_lwpolyline(hexp, close=True, dxfattribs={"layer": "OUTLINE"})
    msp.add_circle((0.0, 0.0), 0.635, dxfattribs={"layer": "SMA_HOLE"})
    doc.saveas(FAB / "sub6_final_board.dxf")


def sha(p):
    import hashlib
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    FAB.mkdir(parents=True, exist_ok=True)
    res = pickle.load(open(OUT / "route13.pkl", "rb"))
    if res["fails"]:
        raise SystemExit(f"UNROUTED leads {res['fails']}")
    pads = res["pads"]
    names = pad_names(res["assign"])
    traces = collect_traces(res)
    d = res["d"]
    conn, drc = checks(res, traces, names)
    print(conn["state"], conn["fails"][:5], flush=True)
    print(drc["state"], drc["min_clearances_mm"], "min RETURN width", drc["min_RETURN_width_mm"], flush=True)
    (FAB / "connectivity.json").write_text(json.dumps(conn, indent=1), encoding="utf-8", newline="\n")
    (FAB / "drc.json").write_text(json.dumps(drc, indent=1), encoding="utf-8", newline="\n")
    (FAB / "traces.json").write_text(json.dumps([dict(kind=t["kind"], net=t["net"], string=t["string"], pos=t["pos"], w=t["w"], poly=[(R3(a), R3(b)) for a, b in t["poly"]]) for t in traces]), encoding="utf-8", newline="\n")
    # RF files: copied byte-for-byte from the 24-pad layout
    for fn in ("cell_table.csv", "drill.csv", "drill.drl"):
        shutil.copyfile(SRC / "fab" / fn, FAB / fn)
    cells = list(csv.DictReader(open(FAB / "cell_table.csv")))
    hexp = [tuple(map(float, p)) for p in d["outline"]]
    write_dxf(cells, traces, pads, names, [(R3(a), R3(b)) for a, b in hexp], res["vias"])
    # string table
    srows = []
    pad_of_string = {sid: k for sid, k in res["assign"].items()}
    for s in d["strings"]:
        sid = s["id"]
        Ll = sum(L.poly_len(t["poly"]) for t in traces if t["kind"] == "link" and t["string"] == sid)
        Ld = L.poly_len(res["leads"][sid])
        Lb = sum(L.poly_len(t["poly"]) for t in traces if t["kind"] == "return" and t["string"] == sid)
        k = pad_of_string[sid]
        srows.append(dict(string=sid, n_diodes=len(s["cells"]), link_trace_mm=R3(Ll), lead_trace_mm=R3(Ld), return_branch_mm=R3(Lb), V_forward_at_10mA=round(len(s["cells"]) * L.V_DIODE, 2), I_mA=L.I_DIODE * 1e3, drive_pad=names[k], drive_pad_index_from_top=k + 1, return_pad="RETURN", return_pad_index_from_top=M.RET_K + 1))
    with open(FAB / "string_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0]))
        w.writeheader()
        w.writerows(srows)
    with open(FAB / "connector_pinout.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pad_index_from_top", "pad", "net", "string", "x_mm", "y_mm", "voltage_at_10mA_V", "current_mA"])
        inv = {k: sid for sid, k in pad_of_string.items()}
        for k, p in enumerate(pads):
            if k == M.RET_K:
                w.writerow([k + 1, "RETURN", "RETURN (common, 12 string cathode ends)", "all", R3(p[0]), R3(p[1]), "0 (reference)", 120])
            else:
                w.writerow([k + 1, names[k], f"drive string {inv[k]}", inv[k], R3(p[0]), R3(p[1]), 19.0, 10])
    # RF unchanged gate: byte compare with the 24-pad files (cell_table/drill from dfa0fc22 lineage, gerbers from 526f98de)
    rf = rf_unchanged()
    (FAB / "rf_unchanged.json").write_text(json.dumps(rf, indent=1), encoding="utf-8", newline="\n")
    print(rf["state"], flush=True)
    meta = dict(label="FINAL_CANDIDATE_13PAD / MODEL_ONLY / NOT_HFSS_VALIDATED (bias re-route; RF unchanged)", gates=dict(connectivity=conn["state"], drc=drc["state"], rf=rf["state"]), verification=dict(min_clearances_mm=drc["min_clearances_mm"], requirement_mm=0.15, antipad_requirement_mm=M.ANTIPAD_MIN, all_ge_0p15=bool(drc["state"] == "DRC_OK"), antipad_crossings=drc["antipad_crossings"],
                min_trace_width_mm=drc["min_trace_width_mm"], min_RETURN_width_mm=drc["min_RETURN_width_mm"]),
                strings=srows, n_cells=len(cells), n_switched=sum(c["class"] == "SWITCHED" for c in cells), n_vias=480, outline_hex_circumradius_mm=R3(d["r_out"]), connector_x_mm=R3(d["conn"][0]), n_pads=M.N_PAD, pad_pitch_mm=M.PITCH, return_pad_index_from_top=M.RET_K + 1,
                pad_names={str(k): names[k] for k in range(M.N_PAD)}, pad_y_mm=[R3(p[1]) for p in pads], n_nets=len(traces),
                note="13-pad bias layout: 12 drive pads S01..S12 (top to bottom, RETURN in the middle), common RETURN net on B_Cu. Cell local frame unchanged: diode anode on island, cathode on disc.")
    (FAB / "layout_summary.json").write_text(json.dumps(meta, indent=1), encoding="utf-8", newline="\n")
    print("EXPORT_13PAD_DONE", flush=True)


def rf_unchanged():
    """RF files of the 13-pad layout must equal the 24-pad ones (cell_table.csv, drill.csv, drill.drl, F_Cu, In1_Cu, PTH) - computed from the same inputs by the same script."""
    out = {}
    ok = True
    pairs = [("fab/cell_table.csv", "fab/cell_table.csv"), ("fab/drill.csv", "fab/drill.csv"), ("fab/drill.drl", "fab/drill.drl")]
    for a, b in pairs:
        same = sha(OUT / a) == sha(SRC / b)
        out[a] = "identical" if same else "DIFFERENT"
        ok &= same
    g = FAB / "gerber"
    for fn in ("sub6_final-F_Cu.gtl", "sub6_final-In1_Cu.g2", "sub6_final-PTH.drl"):
        if (g / fn).exists():
            same = sha(g / fn) == sha(SRC / "fab" / "gerber" / fn)
            out["gerber/" + fn] = "identical" if same else "DIFFERENT"
            ok &= same
        else:
            out["gerber/" + fn] = "not yet generated"
    out["state"] = "RF_UNCHANGED" if ok and all(v != "not yet generated" for v in out.values()) else ("RF_CHANGED" if not ok else "RF_PENDING_GERBER")
    out["reference"] = "24-pad layout results/sub6_final_board/mesh96 (cell table and drill CSV of dfa0fc22/d47822d lineage; Gerbers of 526f98de)"
    return out


if __name__ == "__main__":
    main()
