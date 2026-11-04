# Switchable via-to-ground (mushroom-type) gate cell — estimate (MODEL_ONLY)

Scripts: `scripts/research/sidequest_mushroom_cell.py` (theory + MoM patch sheet), `sidequest_mushroom_package.py` (STEP + predictions); tests `tests/test_sidequest_mushroom.py`;
package `hfss_exports/sidequest_mushroom_cell/`. The dog-bone board `hfss_exports/sidequest_gate_test/` is archived (`ARCHIVED.md`), not run.

## Concept and LC theory
Patch (5.6 mm) on RT5880 1.575 mm, 6.4 mm lattice, a via (r 0.2 mm) from a 0.8 mm pad to the ground; the MA4AGP907 bridges the 0.4 mm gap between pad and patch (in plane, current along x).
ON: forward diode (2.44 ohm), via connected = high-impedance surface; OFF: 28 fF in series with the via, via open, a capacitive patch sheet.
Standard LC model (bibliographic-only, NOT verified from this session: D. Sievenpiper et al., IEEE Trans. Microwave Theory Tech. 47(11), 2059-2074, 1999, doi 10.1109/22.798001):
L = mu0 h = 1.98 nH, C = w eps0 (1 + eps_r)/pi acosh((2w+g)/g), f0 = 1/(2 pi sqrt(LC)), BW = sqrt(L/C)/eta0. For w = 5.6, g = 0.8: C = 0.172 pF, f0 = 8.6 GHz, BW 28 %.
The formula overestimates the patch capacitance by 1.76x against the MoM sheet (C = b/omega = 0.0975 pF at 10.54 GHz), so the design uses the MoM sheet admittance Y_sheet(f) of the patch+pad cell with a
wire-medium (dense-via) slab, Z_slab = j eta_s tan(k_s h) = j 144 ohm (L_eff = 2.18 nH), in series with the measured diode; OFF: plain slab in parallel with the capacitively loaded via (approximation).

## Dimensions and result
| patch w [mm] | ON via-patch resonance f0 [GHz] |
|---|---|
| 5.2 | 12.36 |
| **5.6** | **10.98** |
| 6.0 | 9.45 |
w = 5.6 mm: **ON TM stopband 10.15 (phase-180 TM edge) to 10.98 GHz = 8 % wide, covers 10.54 +- 0.3 GHz (10.24-10.84)**, margin 0.09 GHz below, 0.14 GHz above. **OFF beta/k0 = 1.080 at 10.54 GHz** (target ~1.06-1.065: reflection of a 3-row strip against the
hologram ~ -19 to -25 dB without transition rows; the OFF sheet is capacitive, b = 6.5 mS patch + 2.5 mS via branch). The gap is narrow: the TM band flattens against the resonance; the stopband width is set by the sheet slope, not by BW from the LC formula.
Limits: lumped model (wire-medium assumption a/lambda_slab = 0.33 is borderline), mesh-32 MoM sheet, no TE bands, no attenuation inside the gap (needs a driven strip run), loss of the ON diode only through R = 2.44 ohm; +-5 % in frequency.
The earlier dog-bone failure came from a missing physical blocking mechanism (all roots, `GATE_CELL.md` correction); here the block is the resonance of the via-patch circuit, which a lumped model describes correctly, but only HFSS shows the real band edges.

## Solver plan (planning only, no code): thin-wire vertical-via extension of the MoM (additive)
1. New `src/vias.py`: `ViaLoad(x, y, radius, z_top, z_bottom=ground, Z_series(f))`, default none => every existing result bit-identical; an optional pad/junction basis on the top metal.
2. Unknowns: one vertical current I_z (piecewise-sinusoidal/rooftop, 2-4 sub-segments) per via plus a junction function joining it to the patch rooftops (attachment mode on the pad).
3. Needs the vertical-source / vertical-field spectral Green's functions of the grounded slab (Jz -> Ex, Ey, Ez; Jx, Jy -> Ez): closed forms (TL analogy) but absent from `green_function.py`: new module.
4. Z-matrix blocks: patch-via, via-via (incl. Floquet periodic images with the Bloch phase), via self term (thin-wire radius) + series load Z_diode on the junction (reuse `LumpedLoad.z_of_f`).
5. Reuse the existing Floquet sum with cutoff F = 2n; verify convergence of the slowly decaying vertical terms (Ewald-type acceleration if needed).
6. Effort: about 6-9 working days (GF 2, basis/junction 2, assembly + loads 2, tests/validation 2-3).
7. Tests: no-via bit-identical (`tests/_platform_tol.py`); reciprocity and energy balance (lossless via, ON diode R -> absorbed power = 1 - sum|S|^2); zero-height / open-via limits reduce to the existing solver;
   wire-medium limit vs the analytic mushroom reflection phase (Sievenpiper LC) and vs the HFSS cell of this package (band edges, OFF beta); mesh/Floquet/segment convergence; rotation/mirror symmetry.
8. Risks: Ez / vertical-current coupling convergence; junction singularity and the mesh-dependent via inductance (thin-wire reduced kernel, radius << cell); spatial dispersion of the via lattice;
   the 15 % normal-incidence-tensor error already seen across slots; cost (extra unknowns small, GF evaluation larger); no HFSS evidence yet, so the extension should wait for the eigenmode result of this package.
