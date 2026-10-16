# phase20_r1_loaded_cell -- one-cell HFSS check of the loaded, oblique MoM (Phase 20 R1, item 4)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY).** The loaded-cell and the oblique-incidence MoM have no HFSS number yet; this package
lets you produce one. It is a **demonstration cell, not a design**: all dimensions are parameters of `scripts/research/phase20_r1_demo.py`.

## The cell
Period 10.0 x 10.0 mm; substrate eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground under the substrate; one metal patch
8.0 (x) x 8.0 (y) mm centred in the cell, split by ONE gap of 1.0 mm across x at the centre (two arms), bridged by one lumped
load; f = 10.0 GHz. `cell.step` parts: `Substrate` (z -1.575..0), `Metal` (the two arms, 35 um slab z 0..0.035), `Ground`, `Load_1` (a single
sheet, z = 0.0175 mm, x -0.5..+0.5, y -4..+4 mm, bridging the arms' side faces).

## HFSS setup
1. Import `cell.step`. Assign **Perfect E** to `Metal` and `Ground` (the MoM metal has zero thickness: if you prefer, replace the 35 um slab by sheets at z = 0; the slab
   shifts phases by a fraction of a degree to a few degrees -- state which you used).
2. Boundaries: master/slave on the four side walls of the 10 x 10 mm cell with the **scan angle** theta = 0, 30, 60 deg, phi = 0 (plane of incidence = xz); air box on top
   (>= lambda/4 above the metal) with a **Floquet port** on its top face, modes 1 and 2 (these are TE and TM with respect to the plane of incidence: at phi = 0 TE has E along y and
   TM has E along x, i.e. across the gap). **De-embed** the port to the pattern plane z = 0 (the S11 phases in `mom_reference.csv` are referenced to z = 0).
3. **Load_1**: assign a **Lumped RLC** boundary with the current line along +x across the sheet (from x = -0.5 to +0.5 mm, y = 0). **HFSS RLC boundaries combine R, L and C in
   PARALLEL, so the reference lists only pure-R and pure-C loads**; the MoM models the series impedance R + jwL + 1/(jwC) (identical for a single element). For `short` assign Perfect E to
   `Load_1` (or fill the gap); for `open` delete the sheet (empty gap).
4. Solve at 10 GHz (adaptive, delta S <= 0.01); export S11 for Floquet modes 1 and 2 (magnitude dB and phase), at each angle and load.

## What to compare (`mom_reference.csv`)
`mom_S11_mag_dB`, `mom_S11_phase_deg`: **Richardson-extrapolated** from the 80 x 80 and 100 x 100 meshes (Floquet order 40 / 50; assumed second-order convergence, S(n) = S_inf + c/n^2; the observed order of the
60/80/100 sequence is listed per row in `observed_order_p`; median 2.6, range 1.3-3.1: where p < 2 the p = 2 extrapolation under-corrects and the stated uncertainty is understated); `mom_uncertainty_phase_deg` / `mom_uncertainty_mag_dB` = the difference to the same extrapolation from the 60 and 80 meshes; the raw 100 x 100 and the 20 x 20
values are listed too. MoM S are power-wave normalised, like HFSS Floquet-port S-parameters.
1. **Convention first:** run the bare grounded slab (no metal, no load; `bare_grounded_slab_reference.csv`, analytic) at each angle and mode. If HFSS's TM phase differs from the
   listed value by a constant 180 deg (a sign convention of the TM mode), apply the same shift to every TM entry before comparing; TE should agree directly. Report what you saw.
2. Then the loads. Expected qualitative result (see `results/phase20_r1_demo/README.md`): the TM (E across the gap) phase swings by ~160 deg at theta = 0 (less at larger angles) between short and open,
   and a resistive load absorbs most of the power near R ~ 10-20 ohm at theta = 0. TE is NOT load-independent: the load is a full-width (8 mm) strip, so the odd (TE) mode drives circulating currents
   through it, which absorb up to ~0.7 of the TE power at 60 deg. A narrow real diode would not carry those currents; this check compares the full-width RLC sheet with the full-width MoM strip.
3. Honest limits of the model for this comparison: zero-thickness metal; load region modelled as a conducting strip carrying Z_L along x only; no vias; rooftop MoM. This cell sits near a patch
   resonance and the rooftop discretisation error in PHASE is large: the 20 x 20 mesh is off by up to 25 deg (TE) and 23 deg (TM) from the extrapolated value, the extrapolation itself is
   uncertain by up to 0.7 deg (TE) and 0.4 deg (TM) (listed per row), magnitudes by up to 0.03 dB. Treat an HFSS-vs-MoM phase difference of that order as discretisation, not as a load-model error.
