# HFSS: mushroom gate cell, four ON variants (+ two OFF twins) — eigenmode, phases 90 and 180 only. MODEL_ONLY predictions.

> **EXPORT STANDARD (2026-11-05, `hfss_exports/EXPORT_STANDARD.md`) - applies to this sheet:** (1) **assign material `pec` to every metal object and the ground: STEP imports as vacuum, without this nothing works**; (2) periodic pairs on **single side faces of ONE air box that spans from the ground bottom to the lid** (substrate side walls inside the box; pairs on the air faces above the substrate alone leave the substrate walls PEC); (3) Clean Up Solutions before each run. The STEP files of this package predate the standard (separate `Air` body starting at the substrate top): draw the air box (x, y = +-3.2 mm, z from the ground bottom to the lid) in HFSS and put the pairs on its faces.


Same stack, boundaries and diode representation as `hfss_exports/sidequest_mushroom_cell/HFSS_SETUP.md` (RT5880 1.575 mm, pitch 6.4 mm, PEC `Ground`/`Patch`/`Pad`/`Via`, PEC lid = top face of `Air` at z = +30 mm,
master/slave periodic, y phase 0, eigenmode, min 2 GHz, 8 modes, lid modes removed as before). Sheet `D_ON`: Lumped RLC **R = 2.44 ohm**; `D_OFF`: **C = 28.16 fF**; current along +x across the gap.
All variants: pad 0.8 mm, via top at the pad, diode between pad and patch. Table `variants.csv`.

| file | patch w [mm] | patch-to-patch gap [mm] | diode gap [mm] | via r [mm] |
|---|---|---|---|---|
| `V1_ON.step` (+ `V1_OFF.step`) | 2.4 | 4.0 | 0.4 | 0.2 |
| `V2_ON.step` (+ `V2_OFF.step`) | 2.8 | 3.6 | 0.4 | 0.2 |
| `V3_ON.step` | 3.2 | 3.2 | 0.4 | 0.2 |
| `V4_ON.step` | 2.8 | 3.6 | 0.8 | 0.5 |

**Runs:** ON, phases **90 and 180 deg**, for V1-V4 = **8 runs**. The two OFF twins (V1, V2) at the same two phases are 4 optional extra runs (they discriminate the two models by beta/k0 at 10.54 GHz: M2 1.01 vs M1 1.06).
Report for each run: the lowest TM-like band (Ez at the via/pad, below the light line 11.71 GHz at 90 deg, 23.42 GHz at 180 deg) and, if present, the next band.

**Calibration.** The baseline run (patch 5.6) showed ON f(180) = 4.92 GHz, flat 4.9-5.1 GHz from 90 to 180 deg, versus 10.98 GHz in the first model. Two recalibrated models bracket the answer:
M2 = wire-medium via with L_eff scaled x4.77 (10.4 nH), keeps the flat local-resonance shape (baseline: 4.52 / 4.72 / 4.84 / 4.92 GHz at 90-180 deg; HFSS 4.9-5.1); M1 = plain slab + series L_br = 6.0 nH, reproduces the ON second band (7.66 / 11.2 vs HFSS 7.74 / 11.14 GHz) but disperses.

| variant | M2 f(90) | M2 f(180) | M1 f(90) | M1 f(180) | OFF beta/k0 at 10.54 (M2 / M1) |
|---|---|---|---|---|---|
| V1 | 7.15 | **9.94** | 3.68 | 6.66 | 1.01 / 1.06 |
| V2 | 6.98 | 9.41 | 3.67 | 6.57 | 1.01 / 1.06 |
| V3 | 6.76 | 8.81 | 3.64 | 6.44 | 1.01 / 1.07 |
| V4 | 7.00 | 9.45 (9.01 if L scales with gap and via) | 3.67 | 6.58 | 1.01 / 1.06 |
**Uncertainty:** M2 error on the baseline: ON f(90) -8 %, OFF f(120) +11 %, OFF beta -0.04; allow +-10 % on M2 values; M1 is the low-side alternative (-30 % at 90 deg). Second ON band (90 deg), M1: 7.9 / 11.5 GHz.
**Honest limit:** no variant reaches 10.5-11 GHz at 180 deg in either model; V1 (smallest patch) is the closest (M2 9.9 GHz). The remaining lever is the series inductance (about 4x too high versus the first model), which patch size does not change.
**Reading:** V2 vs V4 differ only in diode gap (0.4 -> 0.8) and via radius (0.2 -> 0.5): equal f means the inductance is not via/gap dominated; a drop (-5 %) means it is. V1-V3 tell whether f scales as 1/sqrt(C_sheet) (M2) or hardly at all (M1).
Falsifiers: ON f(180) of V1 below 6.5 GHz (M1 side) or above 11 GHz; V1-V3 spread (f(180) V1 - V3) below 0.5 GHz (patch size irrelevant) or above 3 GHz.
