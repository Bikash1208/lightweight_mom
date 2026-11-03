# HFSS_TESTBOARD_SETUP - Phase 30 single-feed test board, eps_r 11.2, 10.05 GHz (manual, from STEP; MODEL_ONLY predictions, written before any HFSS run)

Purpose: measure the edge power of the modulated board and separate leakage from monopole direct radiation and edge diffraction by comparing ACTIVE (572 slotted cells, M-tapered) with CONTROL (same board, every cell the unslotted squircle). Run after the library is accepted.

**Files** (mm, origin = board centre on the ground plane): `testboard_active.step`, `testboard_control.step` (~38 MB, 572 patch sheets each), `ports_and_excitations.csv`, `layout.csv` (per cell: x, y, slot L, W, angle, model alpha, beta/k0).
**Bodies**: SUB (solid 100 x 100 x 0.635, pin hole), GND (sheet z = 0, coax opening r 2.05), PATCH (572 sheets at z = 0.635; 90 mm circle of 3.2 mm cells, squircle 3.05 mm, gap 0.15 mm), PIN (solid r 0.635, z 0..6.635), PORT (annulus sheet r 0.635..2.05, z = 0).
**Import**: Modeler > Import (STEP), mm, heal OFF; new project per file; Driven Modal.

**Materials**: SUB = new `RO3010` eps_r **11.2**, tan delta **0.0022**. Region/air = vacuum. GND, PATCH, PIN = **Perfect E** (lossless).
**Region and boundary**: Draw > Region, padding **15 mm on all six sides** (about lambda_0/2 = 14.9 mm; lambda_0 = 29.8 mm at 10.05 GHz) > **Radiation** boundary on the region (or PML; report which).
**Port**: PORT sheet > **Lumped Port**, 50 ohm, renormalise, no de-embedding; integration line radial from the pin (r 0.635) to the opening edge (r 2.05); excite 1 W, 0 deg.
**Solution**: single frequency **10.05 GHz**, max passes **15**, max delta S **0.02**, first-order basis; optional interpolating sweep 9.5 - 10.6 GHz (S11 only). Far field: infinite sphere, theta 0..180 / phi 0..355 step 1 deg (3-5 deg is enough for the first look).
**Edge power (the number we need)**: Draw a cylinder r **43.5 mm**, height 20 mm, base at z = 0 (centred on the board axis; use its lateral face as a sheet, name `EDGE`). Fields Calculator: Poynting vector (Re[E x H*]/2) > integrate over `EDGE` (normal outward) > edge_fraction = that / accepted power. Also read the radiated power and radiation efficiency from the far-field report.

**Model predictions** (ACTIVE; aperture current only, launch fraction and direct radiation NOT modelled):
| quantity | model |
|---|---|
| edge fraction of the launched surface wave | **39 %** (first-order alpha); 60 % (alpha x 0.56, Floquet/first-order ratio); **83 %** if alpha is 4.8x below the model as the S1 HFSS Q suggests; 90 % at alpha/8 |
| aperture efficiency (surface wave -> aperture) | 61 % (first order) |
| RHCP directivity of the aperture current | 19.3 dBi, AR 1.1 dB, peak at theta 0, SLL (outside 24 deg) -19.1 dB |
| RHCP gain upper bound (launch fraction 1) | 17.1 dBi |
| CONTROL | edge fraction ~100 % (alpha = 0 apart from loss), RHCP broadside gain low |

**Copy back to the reviewing agent** (both boards): 
| board | S11 (dB) at 10.05 GHz | accepted power (W) | EDGE flux / accepted (%) | radiation efficiency | RHCP realised gain at theta 0 (dBi) | LHCP at theta 0 | AR at theta 0 (dB) | peak RHCP (dBi, theta, phi) | passes |
|---|---|---|---|---|---|---|---|---|---|
| ACTIVE | | | | | | | | | |
| CONTROL | | | | | | | | | |
The pass/fail reading: ACTIVE edge fraction near 40 % = the model is right; near 80-90 % = the cell leaks about 5-8x less than modelled (see S1).
