"""Writes results/sub6_island_cell/rows.csv and ISLAND_CELL.md from island_rows.json (MODEL_ONLY)."""
import csv
import json
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_island_cell"
J = json.load(open(D / "island_rows.json"))
rows, refs = J["rows"], J["refs"]
c = lambda v: complex(v[0], v[1])  # noqa: E731
hdr = ["id", "W", "L_over_2r", "Li_over_L", "clear_nom_mm", "real_slotW_mm", "real_Li_mm", "real_clear_mm", "real_end_clear_mm", "alpha_OFF", "alpha_ON", "beta_OFF", "beta_ON",
       "Delta_OFF_re", "Delta_OFF_im", "Delta_ON_re", "Delta_ON_im", "Aperp_OFF_abs", "Aperp_ON_abs", "rho_re", "rho_im", "abs_rho_minus_1", "leakage", "resonance", "gate",
       "ref_alpha", "ref_beta", "ref_Delta_abs"]
out = []
for r in rows:
    s, row = r["states"], r["row"]
    ref = refs[f"REF_W{row['W']:g}_q{row['ratio']:g}"]["states"]["NOCUT"]
    rho = c(r["rho"])
    do, dn = c(s["OFF_BAR50"]["Delta_mS"]), c(s["ON_BAR50"]["Delta_mS"])
    out.append([row["id"], row["W"], row["ratio"], row["li_over_l"], row["clear"], f"{row['real']['slot_W_mm']:.3f}", f"{row['real']['island_Li_mm']:.3f}", f"{row['real']['clear_minus_y_mm']:.3f}",
                f"{row['real']['clear_short_ends_mm']:.3f}", f"{r['alpha_off']:.4g}", f"{r['alpha_on']:.4g}", f"{s['OFF_BAR50']['beta_over_k0']:.4f}", f"{s['ON_BAR50']['beta_over_k0']:.4f}",
                f"{do.real:.4g}", f"{do.imag:.4g}", f"{dn.real:.4g}", f"{dn.imag:.4g}", f"{abs(c(s['OFF_BAR50']['A_perp_mS'])):.4g}", f"{abs(c(s['ON_BAR50']['A_perp_mS'])):.4g}",
                f"{rho.real:.4g}", f"{rho.imag:.4g}", f"{r['abs_rho_minus_1']:.4g}", r["leakage"], "FLIP:" + ",".join(r["resonance_flip"]) if r["resonance_flip"] else "none", r["gate"],
                f"{ref['alpha_np_per_m']:.4g}", f"{ref['beta_over_k0']:.4f}", f"{abs(c(ref['Delta_mS'])):.4g}"])
with open(D / "rows.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    w.writerows(out)
n = len(out)
ok = [o for o in out if o[hdr.index("gate")] == "OK"]
cls = {}
for o in out:
    cls[o[hdr.index("leakage")]] = cls.get(o[hdr.index("leakage")], 0) + 1
amax = sorted(out, key=lambda o: -max(float(o[hdr.index("alpha_ON")]), float(o[hdr.index("alpha_OFF")])))[:12]
res = sum(o[hdr.index("resonance")] != "none" for o in out)
tab = "| " + " | ".join(["id", "real slot W / island Li / clear -y (mm)", "alpha OFF / ON", "beta OFF / ON", "Delta OFF / ON (mS)", "|A_perp| OFF / ON", "rho", "|rho-1|", "class", "resonance", "gate"]) + " |\n|" + "---|" * 11 + "\n"
for o in amax:
    g = lambda k: o[hdr.index(k)]  # noqa: E731
    tab += (f"| {g('id')} | {g('real_slotW_mm')} / {g('real_Li_mm')} / {g('real_clear_mm')} | {float(g('alpha_OFF')):.3g} / {float(g('alpha_ON')):.3g} | {g('beta_OFF')} / {g('beta_ON')} | "
            f"{g('Delta_OFF_re')}{float(g('Delta_OFF_im')):+.3g}j / {g('Delta_ON_re')}{float(g('Delta_ON_im')):+.3g}j | {float(g('Aperp_OFF_abs')):.3g} / {float(g('Aperp_ON_abs')):.3g} | "
            f"{float(g('rho_re')):.3g}{float(g('rho_im')):+.3g}j | {float(g('abs_rho_minus_1')):.3g} | {g('leakage')} | {g('resonance')} | {g('gate')} |\n")
md = f"""# SUB6 island cell scan (stage: sub6-island-cell)

**MODEL_ONLY / NOT_HFSS_VALIDATED.** Slotted disc (anchor r 5.435, p 11.6, RT5880 3.175, 5.8 GHz) + floating metal island in the slot + ONE BAR50 sheet (1.00 x 0.35 mm) across the 1.00 mm clearance between the island's +y long edge and the slot's +y edge. Code `scripts/research/sub6_island_cell.py` (+ `sub6_slotdisc_lib.py`, `sub6_slotdisc_grid.py` unchanged); data `rows.csv` ({n} rows), `island_rows.json` (full states, references), `ckpt/` (MoM nodes), `provenance_run.json`.

## Geometry conventions (mine, flagged)
* **Mesh 64 on p = 11.6 mm: one mesh cell = 0.18125 mm** (the mask is sampled at cell centres, so every edge is snapped to a multiple of it). Slot half-length and half-width = floor(dim/2/0.18125 + 1/2) cells; island half-length = round(Li/2/0.18125); clearance 0.3 mm -> 2 cells = **0.3625 mm**, 0.6 mm -> 3 cells = **0.544 mm**; the 1.00 mm +y gap -> 5 cells = **0.906 mm** (same convention as the 1.00 mm rim notch of the cell scan, sheet length 1.00 mm). The realised dimensions are in `rows.csv`. A 0.3 mm gap is only 1.7 mesh cells: resolution-limited, treat the clearance dependence as qualitative.
* The swept clearance applies to the island's -y long edge (the "other side"); the clearance to the two short slot ends follows from Li/L ((L - Li)/2, 1.1-2.7 mm), not from the 0.3/0.6 sweep.
* ISLAND_ABSENT = the closed slotted disc (NOCUT), one reference per (W, L/2r) ({len(refs)} references); W = 4.4 with L/2r = 0.9 skipped (slot corner leaves 0.07 mm of rim, below the 0.2 mm rule of the cell scan).
* States: OFF_BAR50 ((0.15 pF || 5 k) + 0.6 nH), ON_BAR50 (3.0 ohm + 0.6 nH), as in the cell scan. rho = Delta_ON / Delta_OFF (complex). Gate OK: alpha_ON or alpha_OFF >= 1 Np/m AND |rho - 1| >= 1.

## Result
Gate OK: **{len(ok)} of {n}** rows. Leakage classes: {cls}. Resonance proxy flagged: {res} rows.

Top 12 by max(alpha_OFF, alpha_ON):
{tab}
Alpha is the first-order tensor leakage at 25 deg (as the cell scan); beta is the self-consistent TM constant of the mean sheet; the A_perp of the island cell is the antisymmetric part of Y_yy at kt = beta (mS).

## Limits
Mesh 64 (not converged against 96 here), MoM only, no HFSS check; a 0.3 mm clearance is 2 mesh cells; the island is floating metal in a hole and its coupling to the BAR50 model (lumped sheet 1.00 x 0.35 mm) is the only DC/RF connection to the disc; no bias network is modelled. Stopped here: no board model, no HFSS package.
"""
(D / "ISLAND_CELL.md").write_text(md, encoding="utf-8", newline="\n")
print(n, len(ok), cls, res)
