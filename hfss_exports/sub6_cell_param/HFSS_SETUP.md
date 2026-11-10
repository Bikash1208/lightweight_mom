> **STATUS 2026-11-10:** this setup note describes the BAR50 calibration cell and its original draft. Calibration is now reproduced natively
> (VERIFY passes), with the closed PEC box floor at z = -33.21 mm (not -34.045 as written below). Current BAR90 setup: see `RUN.md` (BAR90 section)
> and `REMOTE_RUN.md`.

# Nominal BASELINE cell — setup reference, BLOCKED

**Assign material `pec` to METAL_PEC after STEP import.** This is the identical checked
BASELINE_vias_fr4_stubs STEP from the preceding package, not a newly solved cell.
SUBSTRATE: RT5880 eps_r 2.2, tanD 0.0009; SUBSTRATE_FR4: eps_r 4.4, tanD 0.02; AIR: vacuum.

Period 11.6 mm; disc r 5.435; slot along y, 8.696 x 2.900 mm; island 4.348 y x 1.300 x,
centre (0.200,0); via centres (0.530,-1.780), (3.630,-3.630), drill 0.300, bottom pads/anti-pads
0.600 mm. FR4 0.800 mm; retain both 2.000 x 0.300 mm +y bottom stubs at z -4.045..-4.010.

AIR intentionally overlaps all solids. Four entire side faces form lattice pairs, x phase $ph_x,
y phase $ph_y. Top z150 Perfect E. Calibration bottom z-34.045 is unassigned, hence default PEC.
VERIFY reproduces this closed box. The separate artefact test compares PEC floors 30/60 mm below
the lowest copper and a 377-ohm/square impedance floor if AEDT accepts it in Eigenmode.

ONE impedance boundary on DIODE_SHEET, 1.00 x 0.35 mm at z0.0175. OFF Rs/Xs 5.424/-79.320;
ON 1.05/+7.65 ohm/sq, constant per solve; series 0.6 nH already included. No W/L rescaling.
Physical gap remains 1.00 mm. ACROSS = phase along x; ALONG = phase along y.

After resolution only: eigenmode min4.5 GHz, 8 modes, 0.1% delta, 15 maximum passes.
VERIFY phase90, target f/Q OFF along5.771/147, OFF across6.202/180,
ON along4.835/115, ON across4.905/79. Frequency tolerance1%, Q30%, four cases must pass.
All target values are USER_REPORTED_HFSS, not results computed by this package.

Criteria: flat low-Q (Q<1000) modes excluded from4.8–6.8 GHz; surface-wave Q>=1000;
OFF |beta_across-beta_along|/k0>=0.09; ON anisotropy<=0.2 OFF; fab radial annular ring/clearances>=0.15 mm.
Flatness analysis convention proposed in draft: f-span<=0.20 GHz over85–95 deg; requires field-identified branches.
Phase90 alone cannot qualify these criteria. No beta extrapolation or MoM ON-state result is accepted.
