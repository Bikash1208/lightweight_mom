# phase20_r2_diode_cell -- one-cell HFSS check of the SMP1345-079LF ON and OFF cell (Phase 20 R2)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY), diode values PROVISIONAL_DATASHEET.** Nothing here has an HFSS number yet; this package is meant
to produce one. It is a **demonstration cell, not a design**. Diode values: Skyworks SMP1345 Series datasheet 200046U (2018-06-19): R_s 1.5 ohm
(10 mA, 100 MHz, Table 3), C_T 0.19 pF (1 V, 1 MHz, Table 3; C_T at 0 V is not given), L_s 0.7 nH (SC-79, Table 1); R_p not available (OFF is lossless).
The datasheet covers 10 MHz-6 GHz: these values are **extrapolated** to 8-12 GHz. Details and sources: `results/phase20_r2/README.md`.

## 0. First: the one-line RLC-boundary check (do this before anything else)
HFSS Lumped RLC boundaries are documented here as R, L, C in **parallel**; I could not confirm that for your HFSS version (the one source found does not say).
The package avoids depending on it (every boundary carries ONE element type), but verify the semantic anyway, because it decides whether an R_p could ever be added:
terminate a 50 ohm lumped port with a 1 mm x 0.5 mm sheet carrying one Lumped RLC boundary with **R = 50 ohm and C = 0.1 pF both checked**, and read S11 at 10 GHz:
* parallel R||C: |Z| = 47.7 ohm at -17.4 deg  ->  |S11| = 0.155
* series R+C:   |Z| = 166.8 ohm at -72.6 deg  ->  |S11| = 0.847
Report which one you get. Expected for the package itself: nothing changes either way (single-type boundaries).

## The cell
Period 10 x 10 mm; RT5880 eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground under the substrate; a metal patch 8 x 8 mm (x, y) with a **2.0 mm slot across x** at the centre,
bridged by a **0.5 mm wide strip** (the two 0.5 mm **pads** at the slot edges, the 1.0 mm diode span in the middle). Pads are required: with the boundaries touching the wide arms
the 2-D current crowding changes the OFF phase by tens of degrees against the ideal series model (`results/phase20_r2/README.md`). `cell.step` parts (see `step_parts.json`):
`Substrate` (z -1.575..0), `Metal` (35 um slab z 0..0.035: the two arms plus the pads and the 0.5 mm **island** between the boundaries, three solids), `Ground`,
`Load_1` and `Load_2` (single sheets at z = 0.0175 mm, 0.25 mm long in x, 0.5 mm wide in y, centred at x = -0.375 and +0.375 mm, bridging the metal faces).

## HFSS setup
1. Import `cell.step`. **Perfect E** on `Metal` and `Ground` (the MoM metal has zero thickness; state whether you replaced the 35 um slab by sheets at z = 0).
2. Master/slave on the four side walls, **scan angle theta = 0 and 30 deg**, phi = 0 (plane of incidence xz); air box above (>= lambda/4 at 8 GHz) with ONE **Floquet port**
   on its top face, modes 1 and 2 (TE and TM w.r.t. the plane of incidence: at phi = 0 TE has E along y, TM has E along x, i.e. across the slot). **De-embed** to the pattern plane z = 0.
3. Diode boundaries: **Lumped RLC** on each load sheet, current line along +x across the sheet. One element type each:

| state | Load_1 (x = -0.375 mm) | Load_2 (x = +0.375 mm) |
|---|---|---|
| ON  | R = 1.5 ohm (L, C unchecked) | L = 0.7 nH (R, C unchecked) |
| OFF | C = 0.19 pF (R, L unchecked) | L = 0.7 nH (R, C unchecked) |

4. Frequency sweep 8-12 GHz (0.5 GHz steps or interpolating); adaptive solve at 10 GHz, delta S <= 0.01. Export Floquet-mode S11 (magnitude dB and phase), modes 1 and 2, theta = 0, 30 deg, ON and OFF.
   That is 8 solves-worth of data: 2 states x 2 angles, TE and TM from the same port.

## What to compare (`mom_reference.csv`)
Rows `model = hfss_chain` are the MoM of exactly this geometry (two loads and the island); `ideal_series` is the single-load series model of the same diode (the model used in `results/phase20_r2/DEMO.md`).
`mom_S11_*` are **extrapolated** from the 40 x 40 and 80 x 80 meshes (both aligned with the 0.25 mm strip edges; Floquet order 80 / 160 = 2n, the converged cutoff; the first version used the minimum order and is kept as `mom_reference_mincutoff.csv`) assuming S(n) = S_inf + c/n^p with p = 1.5, the median order
observed at 10 GHz on meshes 40/80/120 (1.2-1.65 per row, `results/phase20_r2/DEMO.md`); `uncertainty_*` is the spread between the p = 1 and p = 2 extrapolations (up to 3.0 deg TE,
3.1 deg TM); the raw 40 and 80 mesh values are listed too. TE is nearly state-independent for this narrow bridge (an HFSS TE ON/OFF difference of more than a few degrees would be news),
but the TE phase itself is the least converged quantity (several degrees between meshes).
1. **Convention first:** run the bare grounded slab (no metal, no loads; `bare_grounded_slab_reference.csv`, analytic) at theta = 0, 30 deg. If HFSS's TM phase differs from the listed value by 180 deg (sign
   convention of the TM mode), apply the same shift to every TM entry before comparing; TE should agree directly. Report what you saw.
2. Then ON and OFF. TM ON-OFF phase difference of the extrapolated MoM prediction (hfss_chain): theta 0 deg: -16 deg at 10 GHz, largest +91 deg at 8 GHz; theta 30 deg: -14 deg at 10 GHz, largest -96 deg at 8.5 GHz. (The OFF series resonance of L_s and C_T is 13.8 GHz nominal.) This cell's contrast is NOT that of the item-4 demo cell (`results/phase20_r2/DEMO.md`, direct 0.5 mm bridge, no pads): the 2 mm slot and the pads change it.
3. **OFF is a sensitive state**: +-20 % on C_T or L_s moves its phase by tens of degrees (up to >100 deg). A difference of that size between HFSS and the MoM is most likely the diode values (C_T at the OFF bias, the
   lead inductance of the pads and island), not a solver error: do the HFSS OFF run with C_T 0.19 pF and also 0.23 pF (+20 %) so the sensitivity can be read directly.
4. Honest limits: zero-thickness metal in the MoM; the island and pads are conductor strips, the loads are series sheets (`transverse_conducts=False`); the chain and the ideal series model differ by up to 113.7 deg
   (OFF, TM) in the MoM itself; the MoM mesh is aligned only at n = 40, 80, 120 (mesh 160 does not fit in 36 GiB)
