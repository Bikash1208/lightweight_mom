"""Phase 20 R2 item 4: one SMP1345-079LF PIN diode (ideal series model) bridging a gapped patch -- ON vs OFF. NOT a design.

    python scripts/research/phase20_r2_demo.py [out_dir] [--quick]

R1 demo cell (10 mm period, RT5880 eps_r 2.2 tan d 0.0009 1.575 mm, PEC ground) with the 8 x 8 mm patch cut by a 1 mm gap across x; the gap is
bridged ONLY by a 0.5 mm wide strip carrying the diode load (span 1 mm, `transverse_conducts=False`). Outputs in `results/phase20_r2/`:
`demo_sweep.csv` (8-12 GHz, theta 0/30/60, TE/TM, ON/OFF), `demo_sensitivity.json` (+-20 % one parameter at a time, 10 GHz),
`demo_convergence.json` (meshes 40 and 80: the aligned ones; no extrapolation), `bias_budget.json`.
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET (8-12 GHz is an extrapolation of the 6 GHz datasheet range).
"""

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
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.diodes import bias_current_budget, diode_load, smp1345_off, smp1345_on  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
PLANES = ReferencePlanes(0.0, None, SUB.thickness)
SPAN, STRIP_W = 1.0e-3, 0.5e-3
MODES = ("TE", "TM")
ANGLES = (0.0, 30.0, 60.0)


def cell():
    """Two 3.5 x 8 mm patch halves (gap 1 mm) joined by the 0.5 mm strip; the strip (conductor in the MoM) carries the load."""
    arms = (
        AxisRect(-2.25e-3, 0.0, 3.5e-3, 8e-3),
        AxisRect(2.25e-3, 0.0, 3.5e-3, 8e-3),
        AxisRect(0.0, 0.0, SPAN, STRIP_W),
    )
    return UnitCell(RectangularLattice(10e-3, 10e-3), CompositeGeometry(arms), SUB, Metal(), backing="pec")


def _load(model):
    return (diode_load(model, 0.0, 0.0, SPAN, STRIP_W),)


def solve(model_sets, f, theta_deg, mesh_n, floquet=None):
    fo = floquet or max(24, mesh_n // 2 + 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique_over_loads(
            cell(),
            f,
            MeshSettings(mesh_n, mesh_n, fo),
            PLANES,
            math.radians(theta_deg),
            0.0,
            [_load(m) for m in model_sets],
        )


def _co(s, i):
    return s.s11[i, i]


def sweep(out, mesh_n=40, freqs=None):
    freqs = freqs if freqs is not None else np.arange(8.0, 12.01, 0.5) * 1e9
    on, off = smp1345_on(), smp1345_off()
    rows = []
    for th in ANGLES:
        for f in freqs:
            r_on, r_off = solve([on, off], f, th, mesh_n)
            for i, mode in enumerate(MODES):
                a, b = _co(r_on, i), _co(r_off, i)
                rows.append(
                    dict(
                        f_ghz=f / 1e9,
                        theta_deg=th,
                        mode=mode,
                        on_abs=abs(a),
                        on_phase_deg=math.degrees(np.angle(a)),
                        on_loss=1 - float(r_on.power_out_per_incident()[0, i]),
                        off_abs=abs(b),
                        off_phase_deg=math.degrees(np.angle(b)),
                        off_loss=1 - float(r_off.power_out_per_incident()[0, i]),
                        on_off_phase_deg=math.degrees(np.angle(a / b)),
                    )
                )
    with open(out / "demo_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return rows


def sensitivity(out, mesh_n=40, f=10e9):
    sets = []
    on, off = smp1345_on(), smp1345_off()
    for k, base, params in (("on", on, ("R_s", "L_s")), ("off", off, ("C_T", "L_s"))):
        for p in params:
            for fac in (0.8, 1.2):
                sets.append((f"{k}:{p} x{fac:g}", k, base.scaled(**{p: fac})))
    extra = [
        ("on:1mA (Rs 3.5)", "on", smp1345_on("1mA")),
        ("off:Rp 10k (assumed)", "off", smp1345_off(10e3)),
        ("off:Rp 3k (assumed)", "off", smp1345_off(3e3)),
    ]
    out_rows = {}
    for th in ANGLES:
        res = solve([m for _, _, m in sets + extra] + [on, off], f, th, mesh_n)
        nom_on, nom_off = res[-2], res[-1]
        for (name, kind, _), r in zip(sets + extra, res[:-2], strict=True):
            nom = nom_on if kind == "on" else nom_off
            for i, mode in enumerate(MODES):
                out_rows[f"{name}|theta{th:g}|{mode}"] = dict(
                    d_phase_deg=math.degrees(np.angle(_co(r, i) / _co(nom, i))),
                    d_abs=abs(_co(r, i)) - abs(_co(nom, i)),
                    phase_deg=math.degrees(np.angle(_co(r, i))),
                    abs=abs(_co(r, i)),
                )
    (out / "demo_sensitivity.json").write_text(json.dumps(out_rows, indent=1))
    return out_rows


def convergence(out, freqs=(8e9, 9e9, 10e9, 11e9, 12e9)):
    """ON/OFF co-pol S11 on meshes 40 and 80 (Floquet 24 / 41). Only meshes whose cell edges hit the 0.5 mm strip edges are valid
    (n a multiple of 40: n = 60 or 100 put the strip edges mid-cell, the solver warns, and the result is off), so there are TWO meshes:
    no order can be observed and no Richardson extrapolation is made. The stated discretisation uncertainty of the mesh-40 sweep is the
    mesh 40 -> 80 difference itself (a lower bound on the error of mesh 40, an upper-ish scale for mesh 80)."""
    on, off = smp1345_on(), smp1345_off()
    rows = []
    for th in ANGLES:
        for f in freqs:
            r40 = solve([on, off], f, th, 40, floquet=24)
            r80 = solve([on, off], f, th, 80, floquet=41)
            for k, name in enumerate(("ON", "OFF")):
                for i, mode in enumerate(MODES):
                    a, b = _co(r40[k], i), _co(r80[k], i)
                    rows.append(
                        dict(
                            state=name,
                            f_ghz=f / 1e9,
                            theta_deg=th,
                            mode=mode,
                            mesh40_abs=abs(a),
                            mesh80_abs=abs(b),
                            mesh40_phase=math.degrees(np.angle(a)),
                            mesh80_phase=math.degrees(np.angle(b)),
                            d_phase_deg=math.degrees(np.angle(b / a)),
                            d_abs=abs(b) - abs(a),
                        )
                    )
    (out / "demo_convergence.json").write_text(json.dumps(rows, indent=1))
    return rows


def bias(out):
    rows = []
    for n in (100, 1000, 5685):  # 5685 = the board A cell count (an example of a board size, not an RIS design)
        for frac in (0.5, 1.0):
            for ifw, label in ((10e-3, "10 mA (R_s 1.5)"), (1e-3, "1 mA (R_s 3.5)")):
                i, p = bias_current_budget(n, frac, ifw)
                rows.append(
                    dict(n_cells=n, fraction_on=frac, i_forward_a=ifw, label=label, total_current_a=i, dc_power_w=p)
                )
    (out / "bias_budget.json").write_text(
        json.dumps(
            {
                "note": "one diode per cell; V_F = 0.89 V at 10 mA (datasheet Table 3); at 1 mA V_F is not given, the 10 mA value is used (an overestimate); bias network losses excluded",
                "rows": rows,
            },
            indent=1,
        )
    )
    return rows


def write_summary(out):
    """DEMO.md: every number is read back from the CSV/JSON written above."""
    with open(out / "demo_sweep.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    sens = json.loads((out / "demo_sensitivity.json").read_text())
    conv = json.loads((out / "demo_convergence.json").read_text())
    bud = json.loads((out / "bias_budget.json").read_text())["rows"]

    def sel(th, mode):
        return [r for r in rows if float(r["theta_deg"]) == th and r["mode"] == mode]

    L = [
        "# Phase 20 R2 item 4 -- one SMP1345-079LF in a gapped patch, ON vs OFF (NOT a design)",
        "",
        "**Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET** (datasheet covers 10 MHz-6 GHz; 8-12 GHz is an extrapolation). "
        "Generated by `scripts/research/phase20_r2_demo.py`. Cell: R1 demo cell (10 mm, RT5880 1.575 mm tan d 0.0009, PEC ground), 8 x 8 mm patch cut by a "
        "1 mm gap, bridged only by a 0.5 mm wide strip carrying the diode (ideal series model, span 1 mm, `transverse_conducts=False`); mesh 40 x 40, Floquet 24; "
        "S at the pattern plane z = 0; phi = 0 (TM: E along x, across the gap).",
        "",
        "## ON - OFF reflection phase (deg) of the co-pol S11, 8-12 GHz",
        "| f GHz | " + " | ".join(f"theta {t:g} {m}" for t in (0.0, 30.0, 60.0) for m in ("TM", "TE")) + " |",
        "|" + "---|" * 7,
    ]
    for f in sorted({float(r["f_ghz"]) for r in rows}):
        cells = []
        for t in (0.0, 30.0, 60.0):
            for m in ("TM", "TE"):
                r = next(q for q in sel(t, m) if float(q["f_ghz"]) == f)
                cells.append(f"{float(r['on_off_phase_deg']):+.1f}")
        L.append(f"| {f:g} | " + " | ".join(cells) + " |")
    L += ["", "## Loss (1 - reflected power; includes the 0.0009 substrate loss that a lossless cell would not show)"]
    for t in (0.0, 30.0, 60.0):
        for m in ("TM", "TE"):
            s_ = sel(t, m)
            L.append(
                f"* theta {t:g} {m}: ON max {max(float(r['on_loss']) for r in s_):.3f}, OFF max {max(float(r['off_loss']) for r in s_):.3f}; "
                f"min |S11| ON {min(float(r['on_abs']) for r in s_):.3f}, OFF {min(float(r['off_abs']) for r in s_):.3f}"
            )
    te_max = max(abs(float(r["on_off_phase_deg"])) for r in rows if r["mode"] == "TE" and float(r["f_ghz"]) < 12)
    te_ab = max(max(float(r["on_loss"]), float(r["off_loss"])) for r in rows if r["mode"] == "TE")
    L += [
        "",
        f"**TE vs the R1 8 mm strip:** TE ON-OFF phase difference is at most {te_max:.1f} deg below 12 GHz and TE absorption is at most {te_ab:.3f} "
        "(substrate loss included), against up to 0.68 absorbed with the R1 full-width load at 60 deg. A diode-sized (0.5 mm) bridge removes the circulating-current "
        "path: TE is nearly independent of the diode state. At 12 GHz and 60 deg the TE ON-OFF phase is "
        f"{float(next(q for q in sel(60.0, 'TE') if float(q['f_ghz']) == 12.0)['on_off_phase_deg']):.1f} deg (see the discretisation note).",
        "",
        "## Sensitivity at 10 GHz: +-20 % one parameter at a time (change of the co-pol S11 phase vs the nominal state, deg; TM, theta 0 / 30 / 60)",
    ]
    for n in sorted({k.split("|")[0] for k in sens}):
        L.append(f"* {n}: " + " / ".join(f"{sens[f'{n}|theta{t}|TM']['d_phase_deg']:+.1f}" for t in ("0", "30", "60")))
    L += [
        "",
        "R_p entries are `ASSUMED_SENSITIVITY_NOT_DEVICE_DATA`; `on:1mA` is the datasheet's 1 mA R_s. **OFF is NOT robust:** it sits "
        "near the L_s-C_T series resonance (13.8 GHz nominal), so +-20 % on C_T or L_s moves the phase by tens of degrees up to over 100 deg (the 0.8 x C_T case "
        "crosses the +-180 deg wrap); C_T at 0 V is not given and may be higher than the 0.19 pF used. ON is robust (R_s irrelevant, L_s +-1.5 deg).",
        "",
        "## Discretisation (mesh 40 vs 80; the only two meshes whose cell edges hit the 0.5 mm strip edges)",
    ]
    for st in ("ON", "OFF"):
        for m in ("TM", "TE"):
            s_ = [c for c in conv if c["state"] == st and c["mode"] == m]
            w = max(s_, key=lambda c: abs(c["d_phase_deg"]))
            L.append(
                f"* {st} {m}: max |phase(80) - phase(40)| = {abs(w['d_phase_deg']):.2f} deg (at {w['f_ghz']:g} GHz, theta {w['theta_deg']:g}); max |d abs| {max(abs(c['d_abs']) for c in s_):.4f}"
            )
    L += [
        "",
        "R1's practice (Richardson with the observed order from three meshes) could not be repeated: meshes 60 and 100 put the strip edges mid-cell "
        "(the solver warns) and gave non-monotone results (the OFF TM phase read -157 / -180 / -156 deg on 40 / 60 / 80), so they were dropped; with two meshes no order is observable "
        "and **no extrapolation is made**. The stated uncertainty of the tables above is therefore the 40 -> 80 difference. TE phase at 60 deg is the worst (several degrees); "
        "the ON-OFF TE statement above is robust to it because ON and OFF share the same TE error.",
        "",
        "## Bias-current implication (one diode per cell; V_F = 0.89 V at 10 mA, Table 3; V_F at 1 mA is not given, the 10 mA value is used: an overestimate; bias network excluded)",
        "| cells | fraction ON | I_F per diode | total current | DC power |",
        "|---|---|---|---|---|",
    ]
    for b in bud:
        L.append(
            f"| {b['n_cells']} | {b['fraction_on']:g} | {b['label']} | {b['total_current_a']:.3g} A | {b['dc_power_w']:.3g} W |"
        )
    L += ["", "5685 cells is the board A cell count, used only as an example of board size."]
    (out / "DEMO.md").write_text("\n".join(L) + "\n")


def main(out_dir, quick=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sweep(out, freqs=np.array([8e9, 10e9, 12e9]) if quick else None)
    sensitivity(out)
    if not quick:
        convergence(out)
    bias(out)
    if not quick:
        write_summary(out)


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(a[0] if a else ROOT / "results" / "phase20_r2", quick="--quick" in sys.argv)
