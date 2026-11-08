"""SUB6 slot-disc DC-block by a second, OPEN printed gap (MODEL_ONLY / IDEAL_BIAS / NOT_HFSS_VALIDATED).  Job sub6-slotdisc-dcblock.

Cell = SlotDiscNotch (diode notch 1.00 mm at the +u slot end, BAR50 sheet) + a second notch of width g2 through the opposite rim (-u slot end),
left OPEN (no component), so the two disc halves are DC-isolated.  g2 in {0.10, 0.15, 0.20, 0.30, 0.50, 1.00} mm, plus g2 = 0 (single notch) as the
same-mesh reference.  RT5880 3.175 mm, p 11.6, 5.8 GHz, rows = the 5 ls_tol rows (CELL1 = g0.73_q0.8_W2.9 first).

MESH: the mask is sampled at cell centres, so a gap is resolved only as an integer number of cell rows.  Mesh 64 (h = 0.181 mm) cannot resolve
0.10-0.15 mm at all, so the run uses an ODD mesh n = 97 (F = 2n = 194, h = 0.1196 mm): the centre row y = 0 is removed for g2 < 0.239 mm (modelled
gap 1 row = 0.120 mm), 3 rows (0.359) for 0.30, 5 rows (0.598) for 0.50, 9 rows (1.076) for 1.00.  The modelled width is reported with every row.
States: ON_BAR50 Ls = 0.5 / 0.6 / 0.7 nH (R 3.0), OFF_BAR50 C = 0.15 pF (5 kohm || C, + 0.6 nH).  One BRIDGE assembly per (row, g2, direction, kt) shared by the
4 loads.  Nodes start at kt = 1.05 (OFF) and the ls_tol ON beta of the row, added adaptively (max 4) until every case beta is within 0.004 k0.

    python sub6_slotdisc_dcblock.py selftest | pilot | run [--workers N] | res | collect
"""
import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_slotdisc_grid as G
import sub6_slotdisc_lib as S
from src.basis_functions import RooftopBasis
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.green_function import LayeredGreenFunction
from src.loads import LumpedLoad, load_matrix
from src.materials import Metal
from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
from src.periodic_bc import RectangularLattice
from src.solver import solve

MM = 1e-3
OUT = G.ROOT / "results/sub6_slotdisc_dcblock"
N = 97
FO = 2 * N
ROWS = ["RT5880_p11.6_g0.73_q0.8_W2.9", "RT5880_p11.6_g1_q0.9_W2.9", "RT5880_p11.6_g1_q0.9_W2.4", "RT5880_p11.6_g0.8_q0.8_W3.6", "RT5880_p11.6_g0.73_q0.8_W3.6"]
G2S = (0.0, 0.12, 0.36, 0.60, 1.08)  # 0 = single-notch reference (same mesh); design values = mesh-realisable widths at n = 97 (1, 3, 5, 9 rows; 2 rows (0.24) is not symmetric about y = 0 on an odd mesh)
CASES = {
    "ON_Ls0.5": dict(R=3.0, L=0.5e-9),
    "ON_Ls0.6": dict(R=3.0, L=0.6e-9),
    "ON_Ls0.7": dict(R=3.0, L=0.7e-9),
    "OFF_C0.15": dict(z_of_f=lambda f: S.z_off_bar50(f, C=0.15e-12)),
}
MAXN = 4


@dataclass(frozen=True)
class SlotDiscNotch2(S.SlotDiscNotch):
    notch2_width: float = 0.0  # open notch from the -x slot end to the disc edge, |y| <= notch2_width / 2

    def in_feature(self, x, y):
        inside = super().in_feature(x, y)
        if self.notch2_width > 0:
            x, y = np.asarray(x, float), np.asarray(y, float)
            inside = inside & ~((x <= -self.slot_length / 2) & (np.abs(y) <= self.notch2_width / 2))
        return inside


def make_cell(row, g2, state="ON_BAR50", bridge_state=True):
    sub = S.make_sub(row["sub"])
    geo = SlotDiscNotch2(row["r"] * MM, row["L"] * MM, row["W"] * MM, S.G_C, bridge=True, notch2_width=g2 * MM)
    loads = (LumpedLoad(geo.diode_x, 0.0, length=S.G_C, width=S.DIODE_ALONG, direction="y", R=3.0, L=0.6e-9),)
    return UnitCell(RectangularLattice(row["p"] * MM, row["p"] * MM), geo, sub, Metal(), backing="pec", loads=loads)


def modelled_gap_mm(row, g2, n=N):
    h = row["p"] / n
    yc = -row["p"] / 2 + (np.arange(n) + 0.5) * h
    return float(h * np.sum(np.abs(yc) <= g2 / 2 + 1e-12))


def case_task(a):
    row, g2, d, kt, n, fo = a
    sub = S.make_sub(row["sub"])
    f = S.F
    k0 = free_space_wavenumber(f)
    cell = make_cell(row, g2)
    kx, ky = (kt * k0, 0.0) if d == "x" else (0.0, kt * k0)
    green = LayeredGreenFunction(sub, f, backing="pec")
    basis = RooftopBasis(cell, n, n)
    t0 = time.time()
    z0 = assemble_impedance_matrix(basis, green, kx, ky, fo)
    v = np.stack([assemble_excitation_vector(basis, e, kx, ky) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(kx, ky), dtype=complex)
    out = {}
    for name, kw in CASES.items():
        ld = LumpedLoad(cell.geometry.diode_x, 0.0, length=S.G_C, width=S.DIODE_ALONG, direction="y", **kw)
        cur = solve(z0 + load_matrix(basis, [ld], green.omega), v)
        j = np.stack([basis.current_coefficients(cur[:, k], kx, ky) for k in range(2)], axis=1)
        out[name] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(row=row["id"], g2=g2, direction=d, kt=kt, Y=out, seconds=time.time() - t0, unknowns=int(basis.n_unknowns))


def res_task(a):
    """Normal-incidence Y at 5.0 and 6.6 GHz for ON (Ls 0.6) and OFF (C 0.15) -> sign flips of Im Yxx / Im Yyy (the cell-scan resonance proxy)."""
    row, g2, fr, n, fo = a
    sub = S.make_sub(row["sub"])
    cell = make_cell(row, g2)
    green = LayeredGreenFunction(sub, fr, backing="pec")
    basis = RooftopBasis(cell, n, n)
    z0 = assemble_impedance_matrix(basis, green, 0.0, 0.0, fo)
    v = np.stack([assemble_excitation_vector(basis, e, 0.0, 0.0) for e in np.eye(2)], axis=1)
    g00 = np.asarray(green.spectral_dyad(0.0, 0.0), dtype=complex)
    out = {}
    for name, kw in CASES.items():
        if name not in ("ON_Ls0.6", "OFF_C0.15"):
            continue
        ld = LumpedLoad(cell.geometry.diode_x, 0.0, length=S.G_C, width=S.DIODE_ALONG, direction="y", **kw)
        cur = solve(z0 + load_matrix(basis, [ld], green.omega), v)
        j = np.stack([basis.current_coefficients(cur[:, k], 0.0, 0.0) for k in range(2)], axis=1)
        out[name] = j @ np.linalg.inv(np.eye(2) + g00 @ j)
    return dict(row=row["id"], g2=g2, f=fr, Y=out)


def get_rows():
    rs = {r["id"]: r for r in G.rows(["RT5880"])}
    return [rs[i] for i in ROWS]


def nf(rid, g2):
    return OUT / f"nodes_{rid}_g2_{g2:g}.json"


def load(rid, g2):
    p = nf(rid, g2)
    if not p.exists():
        return {}
    return {float(k): {d: {c: G.to_c(y) for c, y in cs.items()} for d, cs in v.items()} for k, v in json.load(open(p)).items()}


def save(rid, g2, nodes):
    OUT.mkdir(parents=True, exist_ok=True)
    nf(rid, g2).write_text(json.dumps({str(k): G.jc(v) for k, v in nodes.items()}))


def free_gb():
    return G.free_gb()


def start_nodes(rid):
    ls = json.load(open(G.ROOT / "results/sub6_slotdisc_ls_tol/ls_tol.json"))["results"]
    b = [x for x in ls if x["id"] == rid][0]["cases"]["ON_Ls0.6"]["beta_over_k0"]
    return [1.05, round(float(b), 4)]


def run(workers):
    OUT.mkdir(parents=True, exist_ok=True)
    rows = get_rows()
    st = {}
    for r in rows:
        for g2 in G2S:
            nodes = load(r["id"], g2)
            st[(r["id"], g2)] = dict(nodes=nodes, pending=[] if nodes else start_nodes(r["id"]))
    t0 = time.time()
    mem_gb = 12.0  # n = 97: 2*97^2 ~ 17k unknowns, Z + LU + temporaries (measured in pilot)
    with Pool(workers, maxtasksperchild=1) as pool:
        for rnd in range(MAXN + 1):
            tasks = [(r, g2, d, float(kt), N, FO) for r in rows for g2 in G2S for kt in st[(r["id"], g2)]["pending"] for d in "xy"]
            if not tasks:
                break
            print(f"[{time.time() - t0:7.0f}s] round {rnd}: {len(tasks)} assemblies, {workers} workers", flush=True)
            k = 0
            for o in pool.imap_unordered(case_task, tasks):
                k += 1
                s = st[(o["row"], o["g2"])]
                s["nodes"].setdefault(o["kt"], {})[o["direction"]] = o["Y"]
                if len(s["nodes"][o["kt"]]) == 2:
                    save(o["row"], o["g2"], s["nodes"])
                print(f"[{time.time() - t0:7.0f}s] {k}/{len(tasks)} {o['row'][12:]} g2={o['g2']:g} {o['direction']} kt={o['kt']:.4f} {o['seconds']:.0f}s {o['unknowns']} unk", flush=True)
            for r in rows:
                for g2 in G2S:
                    s = st[(r["id"], g2)]
                    s["pending"] = []
                    if len(s["nodes"]) >= MAXN:
                        continue
                    far = [x["beta_over_k0"] for x in (S.state_summary(s["nodes"], c, r["sub"], guess=1.06) for c in CASES) if x["nearest_node_dist"] > 0.004]
                    if far:
                        kt = round(float(np.median(far)), 4)
                        if all(abs(kt - k_) > 3e-4 for k_ in s["nodes"]):
                            s["pending"] = [kt]
    print("RUN_DONE", flush=True)


def res(workers):
    rows = get_rows()
    p = OUT / "resonance.json"
    done = json.load(open(p)) if p.exists() else {}
    tasks = [(r, g2, fr, N, FO) for r in rows for g2 in G2S for fr in (5.0e9, 6.6e9) if f"{r['id']}|{g2:g}|{fr:g}" not in done]
    with Pool(workers, maxtasksperchild=1) as pool:
        for o in pool.imap_unordered(res_task, tasks):
            done[f"{o['row']}|{o['g2']:g}|{o['f']:g}"] = G.jc(o["Y"])
            p.write_text(json.dumps(done))
            print("res", o["row"][12:], o["g2"], o["f"], flush=True)
    print("RES_DONE", flush=True)


def collect():
    rows = get_rows()
    resj = json.load(open(OUT / "resonance.json")) if (OUT / "resonance.json").exists() else {}
    ref = {x["id"]: x for x in json.load(open(G.ROOT / "results/sub6_slotdisc_ls_tol/ls_tol.json"))["results"]}
    out, lines = [], ["row,g2_mm,modelled_gap_mm,alpha_OFF,alpha_ON_Ls0.6,beta_OFF,beta_ON,absDelta_OFF_mS,absDelta_ON_mS,A_perp_OFF_mS,A_perp_ON_mS,dalpha_dLs_pct_per_0p1nH,resonance_flip,gate"]
    for r in rows:
        for g2 in G2S:
            nodes = load(r["id"], g2)
            if not nodes:
                continue
            cs = {c: S.state_summary(nodes, c, r["sub"], guess=1.06) for c in CASES}
            a = {c: cs[c]["alpha_np_per_m"] for c in CASES}
            dl = (a["ON_Ls0.7"] - a["ON_Ls0.5"]) / 0.2
            rel = abs(100 * dl * 0.1 / a["ON_Ls0.6"]) if a["ON_Ls0.6"] > 0 else float("nan")
            flip = []
            for c in ("ON_Ls0.6", "OFF_C0.15"):
                k5, k6 = f"{r['id']}|{g2:g}|5e+09", f"{r['id']}|{g2:g}|6.6e+09"
                if k5 in resj and k6 in resj:
                    y5, y6 = G.to_c(resj[k5][c]), G.to_c(resj[k6][c])
                    flip += [f"{c}:{'xy'[i]}" for i in (0, 1) if np.sign(y5[i, i].imag) != np.sign(y6[i, i].imag)]
            have_res = f"{r['id']}|{g2:g}|5e+09" in resj
            gate = "REFERENCE" if g2 == 0 else ("DCBLOCK_EQUIVALENT" if (a["ON_Ls0.6"] >= 1 and a["OFF_C0.15"] < 0.1 and rel <= 5.0) else "NOT_EQUIVALENT")
            rec = dict(id=r["id"], g2_mm=g2, modelled_gap_mm=modelled_gap_mm(r, g2) if g2 else 0.0, alpha_OFF=a["OFF_C0.15"], alpha_ON=a["ON_Ls0.6"], alpha_ON_Ls=[a["ON_Ls0.5"], a["ON_Ls0.6"], a["ON_Ls0.7"]],
                       beta_OFF=cs["OFF_C0.15"]["beta_over_k0"], beta_ON=cs["ON_Ls0.6"]["beta_over_k0"], absDelta_OFF=abs(cs["OFF_C0.15"]["Delta_mS"]), absDelta_ON=abs(cs["ON_Ls0.6"]["Delta_mS"]),
                       A_perp_OFF=cs["OFF_C0.15"]["A_perp_mS"], A_perp_ON=cs["ON_Ls0.6"]["A_perp_mS"], dalpha_dLs_pct_per_0p1nH=rel, resonance_flip=flip if have_res else None, nodes=sorted(nodes), gate=gate,
                       extrapolated={c: cs[c]["extrapolated"] for c in CASES})
            out.append(rec)
            lines.append(f"{r['id']},{g2:g},{rec['modelled_gap_mm']:.3f},{rec['alpha_OFF']:.4g},{rec['alpha_ON']:.4g},{rec['beta_OFF']:.4g},{rec['beta_ON']:.4g},{rec['absDelta_OFF']:.4g},{rec['absDelta_ON']:.4g},{rec['A_perp_OFF']:.4g},{rec['A_perp_ON']:.4g},{rel:.3g},{'|'.join(flip) if have_res else 'NA'},{gate}")
    prov = G.provenance(dict(stage="dcblock", mesh=N, floquet_order=FO, rows=ROWS, g2=G2S))
    (OUT / "dcblock.json").write_text(json.dumps(G.jc(dict(provenance=prov, ls_tol_ref_mesh64={k: v["cases"]["ON_Ls0.6"]["alpha_np_per_m"] for k, v in ref.items()}, results=out)), indent=1, default=str))
    (OUT / "dcblock.csv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def selftest():
    rows = get_rows()
    r = rows[0]
    h = r["p"] / N
    yc = -r["p"] / 2 + (np.arange(N) + 0.5) * h
    print("centre row y =", yc[N // 2])
    for g2 in G2S:
        print("g2", g2, "modelled gap rows", int(np.sum(np.abs(yc) <= g2 / 2 + 1e-12)), "=", f"{modelled_gap_mm(r, g2):.3f} mm")
    c0 = make_cell(r, 0.0)
    c1 = make_cell(r, 0.0)
    s = S.make_cell(r["p"], r["r"], r["L"], r["W"], S.make_sub(r["sub"]), "ON_BAR50")
    xs = np.linspace(-r["p"] / 2, r["p"] / 2, 241) * MM
    X, Y = np.meshgrid(xs, xs)
    assert np.array_equal(c0.geometry.in_feature(X, Y), s.geometry.in_feature(X, Y)), "g2 = 0 must equal SlotDiscNotch"
    c2 = make_cell(r, 0.5)
    d = c0.geometry.in_feature(X, Y) & ~c2.geometry.in_feature(X, Y)
    assert d.any() and (X[d] < 0).all() and (np.abs(Y[d]) <= 0.25 * MM + 1e-12).all(), "notch2 must only remove -x material in |y| <= g2/2"
    # DC isolation: with the diode region NOT metal, the disc must split into 2 components for every g2 > 0 (and stay 1 for g2 = 0)
    from scipy.ndimage import label
    hx = r["p"] / N
    cx = -r["p"] / 2 + (np.arange(N) + 0.5) * hx
    GX, GY = np.meshgrid(cx * MM, cx * MM, indexing="ij")
    for g2 in G2S:
        geo = SlotDiscNotch2(r["r"] * MM, r["L"] * MM, r["W"] * MM, S.G_C, bridge=False, notch2_width=g2 * MM)
        lab, nl = label(geo.in_feature(GX, GY))
        print(f"g2={g2:g}: metal components without the diode = {nl}")
        assert nl == (1 if g2 == 0 else 2), "DC isolation check"
    print("selftest ok")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["selftest", "pilot", "run", "res", "collect"])
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--g2", default="", help="comma list of g2 values to run (default all); decision 2026-11-08: 0,0.12 only")
    a = ap.parse_args()
    if a.g2:
        G2S = tuple(float(x) for x in a.g2.split(","))
    if a.stage == "selftest":
        selftest()
    elif a.stage == "pilot":
        rows = get_rows()
        t = case_task((rows[0], 0.10, "x", 1.05, N, FO))
        print("pilot", t["seconds"], t["unknowns"], {k: np.round(v[0, 0], 5) for k, v in t["Y"].items()})
    elif a.stage == "collect":
        collect()
    else:
        nw = a.workers or int(min(61, max(4, free_gb() * 0.7 / 12.0)))
        (res if a.stage == "res" else run)(nw)
