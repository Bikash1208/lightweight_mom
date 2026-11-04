# Round 3 — retune, pass-state match, ring with the real pass state, HFSS gate test (MODEL_ONLY)

## 1. Retuned dog-bone (mesh 32 screen of 36 variants, mesh-64 check; `sidequest_gate_retune.py`, `gate_cell_v2_raw.json`)
New dimensions: end plates **1.6 mm** long x **6.0 mm** high, centre strip **1.2 mm** wide, cell length **6.0 mm** (neighbour gap 0.4 mm), diode gap 0.4 mm, bridge/footprint 0.8 mm
(was 1.2 / 6.0 / 0.8 / 5.6 / 0.4 / 0.8). Every edge on the 0.1 mm grid. Y_xx (ON, mesh 64, 10.54 GHz) = 100.1 + 68.0j mS, (OFF) = 0.001 + 2.5j mS; Y_yy = 6.6j mS both states.
**"No bound TM0 mode" band (the criterion of the request): mesh 64 10.25-10.95 GHz (10.15 bound, 11.0 bound on mesh 32); mesh 32 10.2-10.9/11.0 GHz; centre ~ 10.6 GHz,
margin about +-0.35 GHz around 10.54 GHz (>= +-0.15 GHz met on both meshes); mesh shift ~ 0.05 GHz.** ON R +-20 %: still inside (mesh 32 and 64); OFF C +-20 %: b_xx 2.3-2.7 mS, beta/k0 1.026-1.027.
**This band is NOT a stopband** (see the correction in `GATE_CELL.md`): inside it the damped quasi-TEM mode carries the wave with 2.1-2.9 dB loss over 3 rows.

## 2. Pass-state match with transition rows (`sidequest_gate_match.py`)
One transition row each side (ring = 1 + 3 + 1 rows) of sheet **b_t = 5.64 mS** (a square patch of about 5.4 mm; mesh-64 patch 5.4 mm = 5.44 mS), optimised over the ring's
angular range (x-cells for |psi - axis| <= 45 deg, y-cells beyond; psi_rel = 0-45 deg): **RL 21.9-53.9 dB, IL 0.01-0.03 dB** (target RL >= 15, IL <= 0.5 met). Bare 3 rows: RL 7.2, IL 0.93.
Tolerance: b_t +-10 % -> RL >= 15.5 dB, IL <= 0.13 dB; gate OFF b +-15 % -> RL >= 20.5 dB. 2 transition rows: no gain (optimum went unphysical). The IL excludes metal and diode loss.
The transition rows are static uniform patches (isotropic); the hologram average reactance near the gate must be brought to 5.6 mS (library cell with that mean b near the ring: not checked
against the Phase 26 cell library; the hologram's own cells are slotted circles, the test board uses plain squares).

## 3. State logic (accepted)
State A = gate B ON (block), gate A OFF (pass): quadrant A radiates. State B = gate A ON, gate B OFF. 34-35 diodes ON per state (35 on the test half-ring; 68-70 on the full ring),
the OFF half of the ring stays unbiased. **"All OFF" is not a radiating state:** both gates pass, both holograms are fed (two-beam state); the default after power-up must be commanded.
Each half-ring needs two bias groups (x-axis and y-axis diodes) or one if the bias network joins them.

## 4. Ring front radius and the real pass state (`sidequest_gate_ring.py`, `gate_ring_real_pass.json`, `gate_ring_real_block.json`)
Gate front face r_g = **15.05 mm** (centre of the step-1 window 13.8-16.3), inner transition row 8.65-15.05 mm, 3 gate rows to 34.25 mm, outer transition to 40.65 mm = start of the hologram.
Block state ideal (|Gamma| = 1, T = 0, phase convention of step 1): **S11 = -13.9 dB** (ideal pass -13.5), open-quadrant power 0.745, **gain 20.75 dBi** (ideal pass 20.8; D = 24.56 dBi, SLL -14.5 dB);
window S11 <= -10.6 dB for r_g = 12.0-20.1 mm; the real pass state costs < 0.1 dB.
**Real (leaky) block state** (ON stack S21 -2.9 to -4.9 dB, S11 -10.5 to -7.2 dB, all roots): feed S11 -21.3 dB, open-quadrant power 0.50, blocked-side power 0.23, **side isolation 3.3 dB**
(the unwanted quadrant would carry a beam only ~3 dB below the wanted one). Architecture A with this gate therefore does not work; the ideal-block numbers are an upper bound that the cell does not deliver.
Block phase (hence the radius optimum) is an assumption of the model; the HFSS test measures it.

## 5. Package: `hfss_exports/sidequest_gate_test/` (see its HFSS_SETUP.md); 6. Lossy-block fallback: documented only (step 1); 7. Literature: unchanged, bibliographic-only.
