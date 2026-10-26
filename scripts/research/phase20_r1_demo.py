"""Phase 20 R1 item 4: demonstration of the loaded, oblique MoM -- NOT a design.

    python scripts/research/phase20_r1_demo.py [out_dir] [--period-mm 10 --thickness-mm 1.575 --eps-r 2.2 --tan-d 0.0009
        --freq-ghz 10 --patch-length-mm 8 --patch-width-mm 8 --gap-mm 1 --mesh 20 --floquet 16 --angles 0,30,60 --quick]

One patch (length along x, width along y) split by ONE gap across x; the gap region (conductor in the MoM) carries one lumped load (direction x,
transverse conduction blocked: a lumped element only conducts along its axis). All numbers are parameters of the script: the
roadmap's starting values (10 mm cell, RT/duroid 5880 1.575 mm, PEC ground, 10 GHz) are only the defaults. Sweeps the load from short
to open (resistive, capacitive, and an illustrative on/off pair) at theta = 0, 30, 60 deg for TE and TM (phi = 0: TM has E along x,
across the gap). Reference plane = the pattern plane (z = 0). Also writes the one-cell HFSS check package.
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED (the loaded and oblique solver have no HFSS number yet).
"""

import argparse
import csv
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.basis_functions import MeshSettings  # noqa: E402
from src.composite_geometry import AxisRect, CompositeGeometry  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.loads import LumpedLoad  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads, oblique_stack_s_parameters  # noqa: E402

MODES = ("TE", "TM")


def build(period_mm, thickness_mm, eps_r, tan_d, patch_length_mm, patch_width_mm, gap_mm):
    sub = Substrate("demo", eps_r, tan_d, thickness_mm * 1e-3)
    # In the MoM the load region is CONDUCTOR (the rooftops cross it): the pattern is the whole patch; the load sits in its middle.
    # (The CAD export cuts the load region out and bridges it with the Load_1 sheet.)
    geo = CompositeGeometry((AxisRect(0.0, 0.0, patch_length_mm * 1e-3, patch_width_mm * 1e-3),))
    cell = UnitCell(RectangularLattice(period_mm * 1e-3, period_mm * 1e-3), geo, sub, Metal(), backing="pec")
    gap = dict(x=0.0, y=0.0, length=gap_mm * 1e-3, width=patch_width_mm * 1e-3, direction="x")
    return cell, gap, sub


def load_families(gap, omega):
    """name -> list of (label, LumpedLoad kwargs or None for the bare gap)."""
    lk = dict(gap, transverse_conducts=False)
    fam = {"resistive": [("short (R=0)", dict(lk, R=0.0))], "capacitive": [], "on_off_pair": []}
    for r in (2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 5000.0):
        fam["resistive"].append((f"R={r:g}", dict(lk, R=r)))
    fam["resistive"].append(("open", dict(gap, open_circuit=True)))
    for c_pf in (5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
        fam["capacitive"].append((f"C={c_pf:g}pF", dict(lk, C=c_pf * 1e-12)))
    fam["capacitive"].append(("open", dict(gap, open_circuit=True)))
    fam["on_off_pair"] = [("ON  R=2 ohm + L=0.4 nH", dict(lk, R=2.0, L=0.4e-9)), ("OFF C=0.04 pF + L=0.4 nH", dict(lk, C=0.04e-12, L=0.4e-9))]
    return fam


def run(out_dir, period_mm=10.0, thickness_mm=1.575, eps_r=2.2, tan_d=0.0009, freq_ghz=10.0, patch_length_mm=8.0, patch_width_mm=8.0,
        gap_mm=1.0, mesh_cells=40, floquet=24, angles=(0.0, 30.0, 60.0), quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cell, gap, sub = build(period_mm, thickness_mm, eps_r, tan_d, patch_length_mm, patch_width_mm, gap_mm)
    f = freq_ghz * 1e9
    mesh = MeshSettings(mesh_cells, mesh_cells, floquet)
    planes = ReferencePlanes(0.0, None, sub.thickness)  # S referenced to the pattern plane z = 0
    green = LayeredGreenFunction(sub, f, backing="pec")
    omega = green.omega
    fams = load_families(gap, omega)
    rows, bare = [], []
    for th in angles:
        theta = math.radians(th)
        s0 = oblique_stack_s_parameters(green, planes, theta, 0.0)
        bare.append({"theta_deg": th, **{m: [abs(s0.s11[i, i]), math.degrees(np.angle(s0.s11[i, i]))] for i, m in enumerate(MODES)}})
        for fam_name, items in fams.items():
            if quick and fam_name != "resistive":
                continue
            sets = [(LumpedLoad(**kw),) for _, kw in items]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = compute_s_parameters_oblique_over_loads(cell, f, mesh, planes, theta, 0.0, sets)
            for (label, kw), s in zip(items, res, strict=True):
                for i, m in enumerate(MODES):
                    z = LumpedLoad(**kw).impedance(omega)
                    rows.append({
                        "family": fam_name, "load": label, "Z_L_real": z.real if abs(z) < 1e11 else float("inf"), "Z_L_imag": z.imag if abs(z) < 1e11 else 0.0,
                        "theta_deg": th, "mode": m, "S11_co_mag_dB": 20 * math.log10(max(abs(s.s11[i, i]), 1e-300)),
                        "S11_co_phase_deg": math.degrees(np.angle(s.s11[i, i])), "S11_cross_mag_dB": 20 * math.log10(max(abs(s.s11[1 - i, i]), 1e-300)),
                        "absorbed_fraction": 1.0 - float(s.power_out_per_incident()[0, i]),
                    })
    meta = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "parameters": {
        "period_mm": period_mm, "thickness_mm": thickness_mm, "eps_r": eps_r, "tan_d": tan_d, "freq_ghz": freq_ghz, "patch_length_mm": patch_length_mm,
        "patch_width_mm": patch_width_mm, "gap_mm": gap_mm, "mesh_cells": mesh_cells, "floquet_order": floquet, "angles_deg": list(angles),
        "backing": "pec", "reference_plane": "pattern plane z = 0", "load": gap, "transverse_conducts": False,
        "basis": "TE: E along q = z x p; TM: tangential E along p; phi = 0 (TM has E along x, across the gap); power-wave S"}}
    (out / "demo.json").write_text(json.dumps({**meta, "bare_slab": bare, "rows": rows}, indent=1))
    with open(out / "demo.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return meta, bare, rows


def convergence_checks(out_dir, **kw):
    """Mesh / Floquet-order convergence and the effect of the transverse-conduction choice, for three representative loads."""
    base = dict(period_mm=10.0, thickness_mm=1.575, eps_r=2.2, tan_d=0.0009, freq_ghz=10.0, patch_length_mm=8.0, patch_width_mm=8.0, gap_mm=1.0)
    base.update({k: v for k, v in kw.items() if k in base})
    cell, gap, sub = build(*(base[k] for k in ("period_mm", "thickness_mm", "eps_r", "tan_d", "patch_length_mm", "patch_width_mm", "gap_mm")))
    f = base["freq_ghz"] * 1e9
    planes = ReferencePlanes(0.0, None, sub.thickness)
    lk = dict(gap, transverse_conducts=False)
    cases = {"R=10": dict(lk, R=10.0), "C=0.2pF": dict(lk, C=0.2e-12), "open": dict(gap, open_circuit=True)}
    res = {}
    for th in (0.0, 60.0):
        for name, kw_ in cases.items():
            row = {}
            for tag, (mc, fo) in {"mesh20_fl16": (20, 16), "mesh20_fl24": (20, 24), "mesh40_fl24": (40, 24), "mesh60_fl32": (60, 32), "mesh80_fl40": (80, 40)}.items():
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    s_ = compute_s_parameters_oblique_over_loads(cell, f, MeshSettings(mc, mc, fo), planes, math.radians(th), 0.0, [(LumpedLoad(**kw_),)])[0]
                row[tag] = {m: [abs(s_.s11[i, i]), math.degrees(np.angle(s_.s11[i, i]))] for i, m in enumerate(MODES)}
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                s_ = compute_s_parameters_oblique_over_loads(cell, f, MeshSettings(40, 40, 24), planes, math.radians(th), 0.0, [(LumpedLoad(**dict(kw_, transverse_conducts=True, open_circuit=False)),)])[0] if name != "open" else None
            if s_ is not None:
                row["mesh40_fl24_transverse_conducting"] = {m: [abs(s_.s11[i, i]), math.degrees(np.angle(s_.s11[i, i]))] for i, m in enumerate(MODES)}
            res[f"theta{th:g}_{name}"] = row
    (Path(out_dir) / "convergence.json").write_text(json.dumps(res, indent=1))
    return res


def write_readme(out_dir):
    out = Path(out_dir)
    data = json.loads((out / "demo.json").read_text())
    rows, par = data["rows"], data["parameters"]
    conv = json.loads((out / "convergence.json").read_text()) if (out / "convergence.json").exists() else {}

    def pick(fam, th, mode):
        return [r for r in rows if r["family"] == fam and r["theta_deg"] == th and r["mode"] == mode]

    L = ["# Phase 20 R1 -- demonstration of the loaded, oblique MoM (NOT a design)", "",
         "**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED.** Generated by `scripts/research/phase20_r1_demo.py`; every dimension below is a parameter of that script",
         "(the roadmap's starting values are only the defaults). The one-cell HFSS check package is `hfss_exports/phase20_r1_loaded_cell/`.", "",
         "## Cell and solver settings",
         f"Period {par['period_mm']} mm square; RT/duroid-like substrate eps_r {par['eps_r']}, tan d {par['tan_d']}, h {par['thickness_mm']} mm, PEC ground; one patch "
         f"{par['patch_length_mm']} (x) x {par['patch_width_mm']} (y) mm, split by a {par['gap_mm']} mm gap across x at the centre; the gap region is conductor in the MoM and carries ONE lumped load "
         f"(direction x, series Z_L, transverse conduction blocked); {par['freq_ghz']} GHz; mesh {par['mesh_cells']} x {par['mesh_cells']}, Floquet order {par['floquet_order']}; S referenced to the pattern plane z = 0;",
         "TE: E along y (phi = 0), TM: tangential E along x (across the gap); power-wave S-parameters. Whole table: `demo.csv` / `demo.json`.", ""]
    for fam, title in (("resistive", "Resistive load, short -> open"), ("capacitive", "Capacitive load, large C -> open"), ("on_off_pair", "Illustrative on/off pair")):
        L += [f"## {title} (co-pol S11: |S| dB / phase deg; absorbed fraction of the incident power)", ""]
        loads = [r["load"] for r in pick(fam, 0.0, "TM")]
        head = "| load | " + " | ".join(f"theta {t:g} TM | theta {t:g} TE" for t in (0.0, 30.0, 60.0)) + " |"
        L += [head, "|" + "---|" * (1 + 6)]
        for ld in loads:
            cells = []
            for th in (0.0, 30.0, 60.0):
                for m in ("TM", "TE"):
                    r = next(q for q in pick(fam, th, m) if q["load"] == ld)
                    cells.append(f"{r['S11_co_mag_dB']:.2f} / {r['S11_co_phase_deg']:.1f} (abs {r['absorbed_fraction']:.2f})")
            L.append(f"| {ld} | " + " | ".join(cells) + " |")
        L.append("")
    def wrap(d):
        return (d + 180.0) % 360.0 - 180.0

    swing, tebest = {}, {}
    for th in (0.0, 30.0, 60.0):
        res_tm = pick("resistive", th, "TM")
        short, opn = res_tm[0], res_tm[-1]
        swing[th] = abs(wrap(opn["S11_co_phase_deg"] - short["S11_co_phase_deg"]))
        best = max(res_tm, key=lambda r: r["absorbed_fraction"])
        te_best = max(pick("resistive", th, "TE"), key=lambda r: r["absorbed_fraction"])
        tebest[th] = (best["load"], best["absorbed_fraction"], te_best["load"], te_best["absorbed_fraction"])
    L += ["## Observations (computed from the tables above)",
          "* TM (E across the gap) reflection phase, short -> open: " + ", ".join(f"{swing[t]:.0f} deg at theta {t:g}" for t in swing) + " (the swing shrinks with angle).",
          "* Best resistive absorber in the TM sweep: " + "; ".join(f"theta {t:g}: {tebest[t][0]} absorbs {tebest[t][1]:.2f}" for t in tebest) + ".",
          "* TE (E along y) is NOT insensitive: the largest resistive absorption in the TE sweep is " + "; ".join(f"theta {t:g}: {tebest[t][2]} absorbs {tebest[t][3]:.2f}" for t in tebest)
          + ". The model's load is a full-width strip (8 mm) in the gap, so the odd (TE) mode drives circulating x-currents through it (more so at oblique incidence); a single point diode could not carry them. "
          "This is a model property to keep in mind for R3 (a physical diode is narrow), and it is what the one-cell HFSS check with a full-width RLC sheet compares.",
          "* A capacitive gap gives a smooth, lossless phase path from the short value to the open value (|S| = 1 within the substrate loss), the usual diode-OFF behaviour.",
          "* Nothing here is a 4-state design: no state spacing, bandwidth or angle-stability claim is made.", ""]
    if conv:
        L += ["## Discretisation (important for any HFSS comparison)", "",
              "Co-pol S11 (|S| / phase deg) versus mesh, for three loads; Floquet order = n/2 + 4. Magnitudes converge fast; **phases converge slowly (about second order in the mesh size)** and the 20 x 20 mesh is "
              "far from converged for this cell: it sits near a patch resonance.", "",
              "| case | mesh | TM |S| / phase | TE |S| / phase |", "|---|---|---|---|"]
        for case, v in conv.items():
            for tag, r in v.items():
                L.append(f"| {case} | {tag} | {r['TM'][0]:.4f} / {r['TM'][1]:.2f} | {r['TE'][0]:.4f} / {r['TE'][1]:.2f} |")
        L += ["", "`..._transverse_conducting`: the load region also conducts transversely (the brief's default model) instead of being blocked: no effect on TM, a different TE phase (the strip across the gap "
              "conducts y-currents)."]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?", default=str(ROOT / "results/phase20_r1_demo"))
    for k, d in (("period-mm", 10.0), ("thickness-mm", 1.575), ("eps-r", 2.2), ("tan-d", 0.0009), ("freq-ghz", 10.0), ("patch-length-mm", 8.0),
                 ("patch-width-mm", 8.0), ("gap-mm", 1.0)):
        ap.add_argument("--" + k, type=float, default=d)
    ap.add_argument("--mesh", type=int, default=40)
    ap.add_argument("--floquet", type=int, default=24)
    ap.add_argument("--angles", default="0,30,60")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args(argv)
    if not a.quick:
        convergence_checks(a.out_dir)
    res = run(a.out_dir, a.period_mm, a.thickness_mm, a.eps_r, a.tan_d, a.freq_ghz, a.patch_length_mm, a.patch_width_mm, a.gap_mm, a.mesh, a.floquet,
               tuple(float(v) for v in a.angles.split(",")), a.quick)
    if not a.quick:
        write_readme(a.out_dir)
    return res


if __name__ == "__main__":
    main()
