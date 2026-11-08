"""SUB6 island-cell HFSS package, part 2: predictions.json, states.csv, results_template.csv, HFSS_SETUP.md (reads hfss_pred/pred_table.json).  MODEL_ONLY.
    python sub6_island_hfss_docs.py
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import sub6_island_hfss_package as pk
import sub6_slotdisc_hfss_docs as sd
import sub6_slotdisc_hfss_package as sp

OUT, RES = pk.OUT, pk.RES
PHASES_SCAN = range(60, 111, 5)
STATE_LABEL = {"NOCUT": "ABSENT", "OFF_BAR50": "OFF", "ON_BAR50": "ON"}
BOARD_USE = "ISL_W2.9_q0.8_Li0.5_c0.6 is the island row used by 4 of the 6 best board designs of the follow-up (200 mm N 135 b; 220 mm N 120 a, N 135 a, N 240 a; 200 mm N 120 a uses Li0.7_c0.6, 200 mm N 240 a with no ring exclusion uses Li0.5_c0.3)"


def main():
    tab = json.loads((RES / "hfss_pred" / "pred_table.json").read_text())
    eq = {"OFF": sp.parallel_equiv(sp.z_off(pk.F0), pk.F0), "ON": sp.parallel_equiv(sp.z_on(pk.F0), pk.F0)}
    eq_f = {f"{f:.1f}": {"OFF": sp.parallel_equiv(sp.z_off(f * 1e9), f * 1e9), "ON": sp.parallel_equiv(sp.z_on(f * 1e9), f * 1e9)} for f in (5.2, 5.5, 5.8, 6.1, 6.4)}
    J = json.load(open(pk.ROOT / "results" / "sub6_island_cell" / "island_rows.json"))
    grid = [x for x in J["rows"] if x["row"]["id"] == "ISL_W2.9_q0.8_Li0.5_c0.6"][0]
    cell = dict(row_id=grid["row"]["id"], real_mesh_dims_mm=grid["row"]["real"], mean_sheet_5p8GHz_mesh64={st: dict(beta_over_k0=grid["states"][st]["beta_over_k0"], alpha_np_per_m=grid["states"][st]["alpha_np_per_m"]) for st in ("OFF_BAR50", "ON_BAR50")},
                gate=grid["gate"], leakage=grid["leakage"], rho=grid["rho"], resonance_flip=grid["resonance_flip"])
    pred = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "PREDICTED_BEFORE_HFSS_RUN"], "mesh_phase_vs_f": 64, "floquet_order": 128, "p_mm": pk.P,
            "series_models": {"OFF": "(0.15 pF || 5 kohm) + 0.6 nH", "ON": "3.0 ohm + 0.6 nH"}, "parallel_equivalents_5p8GHz": eq, "parallel_equivalents_vs_f": eq_f,
            "cell": cell, "board_use": BOARD_USE, "states": {}, "runs": [],
            "method": "directional bound-mode condition: beta_d from the TM sheet admittance Y_dd(kt = beta_d) alone (d = x along the slot, y across), linear in kt between adaptive nodes; phase = beta_d * p; open air, no lid; f(phase) quadratic fit of 5 frequencies 5.2-6.4 GHz",
            "model_vs_hfss_geometry": "the MoM geometry is the mesh-64 mask of the island scan (clearance 0.544, +y gap 0.906, island width 1.450, length 4.350, slot 8.700 x 2.900 mm); the HFSS geometry uses the nominal dimensions (0.600, 1.000, 1.300, 4.348, 8.696 x 2.900). Differences of up to 0.15 mm in the island width and the gaps are expected to shift the HFSS frequencies from the predictions; the variant-2 predictions are the NO-VIA values (the MoM has no vias).",
            "variant2_note": "variant 2 runs use the same phases as variant 1; f_MoM is the predicted NO-VIA frequency; the difference HFSS(vias) - HFSS(no vias) is the via/anti-pad/stub effect"}
    fits = {}
    for cls, sts in (("ISL", ("OFF_BAR50", "ON_BAR50")), ("REF", ("NOCUT",))):
        for st in sts:
            for d in ("x", "y"):
                tr_all = tab[f"{cls}|{st}|{d}"]
                tr = [r for r in tr_all if r["nearest_node_dist"] < 0.02]
                dropped = [r["f_GHz"] for r in tr_all if r["nearest_node_dist"] >= 0.02]
                f = np.array([r["f_GHz"] for r in tr])
                ph = np.array([r["phase_deg"] for r in tr])
                # branches: a drop of the phase with rising f by more than 5 deg marks a resonance (ON across y: pole near 5.65 GHz); each monotone segment is fitted separately
                cuts = [0] + [i for i in range(1, len(f)) if ph[i] < ph[i - 1] - 5.0] + [len(f)]
                branches = []
                for i0, i1 in zip(cuts[:-1], cuts[1:]):
                    ff, pp = f[i0:i1], ph[i0:i1]
                    fit = np.polyfit(pp, ff, min(2, len(ff) - 1)) if len(ff) > 1 else None
                    branches.append(dict(f_range_GHz=[float(ff.min()), float(ff.max())], phase_range_deg=[float(pp.min()), float(pp.max())], fit=fit, f=ff, ph=pp))
                main = next((br for br in branches if br["f_range_GHz"][0] - 1e-9 <= 5.8 <= br["f_range_GHz"][1] + 1e-9), branches[0])
                fits[(st, d)] = branches
                phi58 = float(np.interp(5.8, main["f"], main["ph"]))
                phic = int(5 * round(phi58 / 5))
                lo, hi = main["phase_range_deg"][0] - 2.0, main["phase_range_deg"][1] + 2.0
                rec = [q for q in (phic - 5, phic, phic + 5) if lo <= q <= hi]
                for extra in (phic + 10, phic - 10):
                    if len(rec) < 3 and lo <= extra <= hi:
                        rec.append(extra)
                rec = sorted(rec)
                pred["states"][f"{STATE_LABEL[st]}|{d}"] = dict(
                    beta_over_k0_5p8GHz=float(np.interp(5.8, f, [r["beta_over_k0"] for r in tr])), phase_deg_5p8GHz=phi58, n_branches=len(branches),
                    branches=[dict(f_range_GHz=br["f_range_GHz"], phase_range_deg=br["phase_range_deg"], contains_5p8GHz=br is main) for br in branches],
                    phase_vs_f=[dict(f_GHz=float(a), beta_over_k0=r["beta_over_k0"], phase_deg=float(b)) for a, b, r in zip(f, ph, tr)],
                    f_at_phase_GHz={str(q): {("main" if br is main else "other"): float(np.polyval(br["fit"], q)) for br in branches if br["fit"] is not None and br["phase_range_deg"][0] - 3 <= q <= br["phase_range_deg"][1] + 3} for q in PHASES_SCAN},
                    recommended_phases_deg=rec, tier1_phase_deg=phic, max_nearest_node_dist=max(r["nearest_node_dist"] for r in tr), dropped_unconverged_f_GHz=dropped,
                    note=("RESONANCE: two branches; the 5.8 GHz branch is the upper one; the lower branch (island resonance, beta/k0 up to 1.39) gives a LOWER mode at the same phase: record modes 1-3" if len(branches) > 1 else ""))
    bare = {str(q): dict(open_air=sd.bare_f(q, None), lid_150mm=sd.bare_f(q, pk.LID), lid_100mm=sd.bare_f(q, 100.0)) for q in (80, 90, 100)}
    pred["sanity_bare_substrate_GHz"] = bare
    runs = []
    # priority: 1 = the OFF/ON difference across y (variant 1) + reference; 2 = the same with vias; 3 = along x variant 1; 4 = along x variant 2; 5 = tier 2 (+-5 deg) of everything
    plan = [("novia", "ISL_novia.step", "ABSENT", "ISL_absent.step"), ("novia", "ISL_novia.step", "OFF", "ISL_novia.step"), ("novia", "ISL_novia.step", "ON", "ISL_novia.step"),
            ("vias", "ISL_vias.step", "OFF", "ISL_vias.step"), ("vias", "ISL_vias.step", "ON", "ISL_vias.step")]
    for d in ("y", "x"):
        for variant, _, st, fnm in plan:
            e = pred["states"][f"{st}|{d}"]
            for i, q in enumerate(e["recommended_phases_deg"]):
                tier = 1 if q == e["tier1_phase_deg"] else 2
                if tier == 2:
                    pri = 5
                elif d == "y":
                    pri = 1 if variant == "novia" else 2
                else:
                    pri = 3 if variant == "novia" else 4
                stl = {"ABSENT": "NOCUT", "OFF": "OFF_BAR50", "ON": "ON_BAR50"}[st]
                brs = fits[(stl, d)]
                main = next((br for br in brs if br["f_range_GHz"][0] - 1e-9 <= 5.8 <= br["f_range_GHz"][1] + 1e-9), brs[0])
                oth = [br for br in brs if br is not main and br["fit"] is not None and br["phase_range_deg"][0] - 3 <= q <= br["phase_range_deg"][1] + 3]
                runs.append(dict(variant=variant, state=st, file=fnm, direction="across_y" if d == "y" else "along_x", ph_x_deg=q if d == "x" else 0, ph_y_deg=q if d == "y" else 0, tier=tier, priority=pri,
                                 f_MoM_GHz=float(np.polyval(main["fit"], q)), f_MoM_other_branch_GHz=float(np.polyval(oth[0]["fit"], q)) if oth else None))
    runs.sort(key=lambda r: (r["priority"], r["variant"] != "novia", {"OFF": 0, "ON": 1, "ABSENT": 2}[r["state"]], r["direction"], r["ph_x_deg"] + r["ph_y_deg"]))
    for n, r in enumerate(runs, 1):
        r["run_id"] = f"I{n:02d}"
    pred["runs"] = runs
    pred["sanity_runs"] = [dict(run_id="S1", what="bare substrate (ground + substrate only, no metal on top), $ph_y = 90, $ph_x = 0", f_model_GHz=bare["90"]["open_air"], priority=1,
                                note="model 6.291 GHz open air / 6.291 GHz with the 150 mm lid; run once; a different value means a setup error (units, lid, periodic pairs)")]
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    sh0, sh1 = pk.sheet_current_line()
    with open(OUT / "states.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["variant", "state", "step_file", "DIODE_SHEET", "R_ohm", "C_pF", "L_nH", "series_model", "note"])
        w.writerow(["reference", "ABSENT", "ISL_absent.step", "none (no object)", "", "", "", "", "closed slotted disc, no island, no sheet"])
        for variant, fnm in (("novia", "ISL_novia.step"), ("vias", "ISL_vias.step")):
            w.writerow([variant, "OFF", fnm, "Lumped RLC parallel", f"{eq['OFF']['R_ohm']:.4g}", f"{eq['OFF']['C_pF']:.4g}", "unticked", "(0.15 pF || 5 kohm) + 0.6 nH", f"5.8 GHz parallel equivalent; current line {sh0} -> {sh1}"])
            w.writerow([variant, "ON", fnm, "Lumped RLC parallel", f"{eq['ON']['R_ohm']:.4g}", "unticked", f"{eq['ON']['L_nH']:.4g}", "3.0 ohm + 0.6 nH", "5.8 GHz parallel equivalent; same geometry as OFF"])
    with open(OUT / "results_template.csv", "w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["run_id", "priority", "variant", "state", "step_file", "direction", "ph_x_deg", "ph_y_deg", "tier", "f_MoM_noVia_GHz", "f_MoM_other_branch_GHz", "f_HFSS_GHz", "Q", "Ez_decays_from_patch_YN", "mode_no", "note"])
        w.writerow(["S1", 1, "bare", "none", "(ground + substrate only)", "across_y", 0, 90, 1, f"{bare['90']['open_air']:.3f}", "", "", "", "", "", "sanity: bare substrate"])
        for r in runs:
            w.writerow([r["run_id"], r["priority"], r["variant"], r["state"], r["file"], r["direction"], r["ph_x_deg"], r["ph_y_deg"], r["tier"], f"{r['f_MoM_GHz']:.3f}", "" if r["f_MoM_other_branch_GHz"] is None else f"{r['f_MoM_other_branch_GHz']:.3f}", "", "", "", "",
                        "ON across y: two MoM branches, the lower one is the island resonance: record modes 1-3" if r["f_MoM_other_branch_GHz"] is not None else ""])
    setup(pred, eq, bare)
    print("written; runs", len(runs))


def setup(pred, eq, bare):
    light90 = sd.C0 * math.radians(90) / (pk.P * 1e-3) / (2 * math.pi) / 1e9
    sh0, sh1 = pk.sheet_current_line()
    s = pred["states"]
    dy = {k: s[f"{k}|y"]["tier1_phase_deg"] for k in ("ABSENT", "OFF", "ON")}
    dx = {k: s[f"{k}|x"]["tier1_phase_deg"] for k in ("ABSENT", "OFF", "ON")}

    def fmain(k):
        br = [b_ for b_ in s[f"{k}|y"]["branches"] if b_["contains_5p8GHz"]][0]
        pts = [(r["f_GHz"], r["phase_deg"]) for r in s[f"{k}|y"]["phase_vs_f"] if br["f_range_GHz"][0] - 1e-9 <= r["f_GHz"] <= br["f_range_GHz"][1] + 1e-9]
        fa, pa = zip(*pts)
        return float(np.interp(dy[k], pa, fa))

    fy = {k: fmain(k) for k in ("ABSENT", "OFF", "ON")}
    on_note = s["ON|y"]["note"]
    txt = f"""# HFSS eigenmode: SUB6 island cell (RT5880 3.175 mm, p 11.6, 5.8 GHz), variants without / with bias vias. EXPORT STANDARD (2026-11-05). Predictions MODEL_ONLY (`predictions.json`), NOT_HFSS_VALIDATED.

**Cell.** Row `ISL_W2.9_q0.8_Li0.5_c0.6`: disc r 5.435, slot 8.696 x 2.900 along x, floating metal island 4.348 x 1.300 mm centred in x, y from -0.850 to +0.450 (0.60 clearance to the slot -y edge, 1.00 mm gap to the slot +y edge), ONE `DIODE_SHEET` bridging the 1.00 mm gap. {BOARD_USE}; the 0.6 mm clearance row was preferred as asked.
**1. Assign material `pec` to `METAL_PEC` (all copper: disc patch, island, ground; variant 2 includes the via barrels, pads and stubs). STEP imports as vacuum; without this nothing works.** (If HFSS splits `METAL_PEC` into several parts, assign all of them: 3 solids, the island is separate and floating.)
**2. Air box `AIR`: one vacuum box, ONE face per side wall.** Variant 1 and the reference: z = -3.210 ... +{pk.LID:g}. Variant 2: z = {pk.Z_AIR_V2} ... +{pk.LID:g} (1.0 mm below the bottom copper; the unassigned bottom face is PEC). It overlaps all solids on purpose; keep `AIR` first. **Periodic pairs go on the four side faces of `AIR` only** (x = -5.8 / +5.8, y = -5.8 / +5.8).
**3.** `SUBSTRATE`: RT5880, eps_r 2.2, tan d 0.0009 (z -3.175 ... 0). Variant 2 only: `SUBSTRATE_FR4` (z -3.710 ... -3.210): **FR-4, eps_r {pk.FR4_EPS}, tan d {pk.FR4_TAND} (typical values ASSUMED: replace with your laminate)**. **4. PEC lid:** Perfect E on the top face of `AIR` (z = +{pk.LID:g}). Bare substrate (MoM, 90 deg): open air {bare['90']['open_air']:.3f} GHz, lid 150 mm {bare['90']['lid_150mm']:.3f} GHz.
**5.** Project Variables `$ph_x`, `$ph_y` (deg). Lattice Pair 1: x faces of `AIR`, `$ph_x`; Pair 2: y faces, `$ph_y`. **ALONG the slot (x): `$ph_x = phi`, `$ph_y = 0`; ACROSS the slot (y): `$ph_x = 0`, `$ph_y = phi`.**
**6. Files / states** (`states.csv`):
| state | file | `DIODE_SHEET` (Lumped RLC, PARALLEL) |
|---|---|---|
| ABSENT (reference) | `ISL_absent.step` | none (closed slotted disc) |
| OFF | `ISL_novia.step` (variant 1) / `ISL_vias.step` (variant 2) | **R = {eq['OFF']['R_ohm'] / 1e3:.3f} kohm, C = {eq['OFF']['C_pF']:.4f} pF, L unticked** |
| ON | same files | **R = {eq['ON']['R_ohm']:.1f} ohm, L = {eq['ON']['L_nH']:.4f} nH, C unticked** |

Parallel equivalents are exact at 5.8 GHz only (OFF series model (0.15 pF || 5 kohm) + 0.6 nH = {eq['OFF']['Z_series'][0]:.2f} {eq['OFF']['Z_series'][1]:+.2f}j ohm; ON 3.0 ohm + 0.6 nH = {eq['ON']['Z_series'][0]:.2f} {eq['ON']['Z_series'][1]:+.2f}j ohm); other frequencies in `predictions.json` (`parallel_equivalents_vs_f`).
**Sheet `DIODE_SHEET`: 1.00 (across, y: +0.450 ... +1.450) x 0.35 (along, x: -0.175 ... +0.175) mm at z = +0.0175, x = 0; its long edges... no: its SHORT ends sit on the island +y edge (y = 0.450) and on the slot +y edge (y = 1.450).** Current line across the gap: ({sh0[0]:g}, {sh0[1]:g}, {sh0[2]:g}) -> ({sh1[0]:g}, {sh1[1]:g}, {sh1[2]:g}).
**Variant 2 vias** (drill 0.3, pad 0.6): island via at (0, {pk.Y_VIA_I}), disc via at (0, +{pk.Y_VIA_D}); both go from the top metal down through the substrate, through 1.0 mm anti-pads in the ground, through the FR-4, to a 0.6 mm pad on the bottom copper (35 um, z -3.745 ... -3.710) with a 2.0 x 0.3 mm stub (island stub towards +x, disc stub towards -x); the two stubs are not connected to each other (2.55 mm apart); top pads lie inside the metal (omitted). Via barrels are PEC cylinders, fused with the island / the disc patch.
**7.** Eigenmode: min frequency 4.5 GHz, 6 modes, delta f 0.1 %, max 15 passes; initial mesh 1 mm on `SUBSTRATE`, 0.5 mm on the patch faces, 0.1 mm on `DIODE_SHEET`, the gaps and (variant 2) the via barrels and anti-pads; `SUBSTRATE_FR4` 1 mm. **Clean Up Solutions before each run.** Patch 35 um (z 0 ... 0.035), ground 35 um (z -3.210 ... -3.175).
**Mode: mode 1 (lowest, TM-like bound mode below the light line; light line at 90 deg {light90:.3f} GHz).** Record f and Q; check Ez over the patch, decaying upward (a higher mode with Ez ~ 0 is TE-like: do not use).
**Run order = `priority` in `results_template.csv`** (1 first). Tier 1 = the phase where the model puts the mode at 5.8 GHz (rounded to 5 deg); tier 2 = +-5 deg (priority 5), to bracket 5.8 GHz if tier 1 is more than about 0.15 GHz off. **Priority 1: S1 bare substrate, then ACROSS y: variant 1 OFF ({dy['OFF']} deg), ON ({dy['ON']} deg) and the ABSENT reference ({dy['ABSENT']} deg);** priority 2: OFF and ON across y in variant 2 (vias) at the same phases; priority 3-4: ALONG x ({dx['OFF']} deg) variants 1 and 2.
**What the model expects (MoM, no vias, mesh 64, mesh-snapped dimensions).** ACROSS y at the tier-1 phases (ABSENT {dy['ABSENT']} deg, OFF {dy['OFF']} deg, ON {dy['ON']} deg, from linear interpolation of the 5.8 GHz branch): ABSENT {fy['ABSENT']:.3f} GHz, OFF {fy['OFF']:.3f} GHz, ON {fy['ON']:.3f} GHz (the `f_MoM` of the runs are quadratic fits; same phase, different numbers: compare at equal phase in `results_template.csv`) (`predictions.json` `states`, `runs`). **If HFSS shows OFF = ON across y, the sheet or the RLC is not active.** **ON across y has a RESONANCE in the model (island + closed diode, pole near 5.65 GHz, beta/k0 up to 1.39 at 108 deg): at the same phase there are two modes, the lower one at about 5.2-5.6 GHz (island resonance, phases 80-109 deg) and the TM mode near the light line above 5.7 GHz (beta/k0 1.00-1.02); the 5.8 GHz branch is the UPPER one, so for the ON state HFSS mode 1 will probably be the LOWER mode: record modes 1-3 with Ez and compare with `f_MoM_noVia_GHz` (upper) and `f_MoM_other_branch_GHz` (lower).** The variant-2 predictions are the same numbers (no vias in the MoM): the difference HFSS(vias) - HFSS(no vias) at the same phase is the via/anti-pad/stub effect; report it per run.
**Known model-vs-HFSS differences.** MoM geometry is the mesh-64 mask (clearance 0.544, gap 0.906, island 4.350 x 1.450, slot 8.700 x 2.900) while HFSS uses nominal 0.600 / 1.000 / 4.348 x 1.300 / 8.696 x 2.900; the MoM lumped sheet is 1.00 mm long over a 0.906 mm gap. Frequencies can differ by more than the notch-cell package (0.1-0.3 GHz is possible); the OFF-ON ordering and sign across y are the robust check. FR-4 values are assumed. Not run in AEDT by the author.
"""
    txt = txt.replace("its long edges... no: its SHORT ends sit", "its SHORT ends sit")
    (OUT / "HFSS_SETUP.md").write_text(txt, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
