# 3-port HMSA — HFSS verification (regenerated 2026-10-06)

10.54 GHz · Rogers RT/duroid 5880 1.575 mm · 4 mm cells · surface-wave (eigenmode) method · square patches 0.5–3.5 mm.
Board 416 × 120 mm (104 × 30 cells, 2688 patches), centred at (0, 0); three 120 mm sub-apertures with 28 mm gaps.
An 8 × 8 mm keep-out (no patches) surrounds every monopole.

| Port | Monopole (x, y) mm | Beam (θ, φ) | Prediction |
|---|---|---|---|
| 1 | (−148, 0) | (25°, 180°) RHCP | (24°, 182°), 20.5 dBi |
| 2 | (0, 0) | (25°, 90°) RHCP | (24°, 89°), 20.5 dBi |
| 3 | (148, 0) | (25°, 0°) RHCP | (24°, 0°), 20.5 dBi |

STEP (mm): `Substrate` → RT/duroid 5880; `Patches` (zero-thickness sheets) → Perfect E; `Ground` → Perfect E.
Setup: three copies of your monopole at the positions above; air box ≥ 15 mm around and ≥ 30 mm above (above the
monopole tips); Radiation/PML; driven at 10.54 GHz; far field one port at a time (Edit Sources).
Check: RealizedGainRHCP cuts (port 1 at θ = −25° in the φ = 0° cut, port 2 in the φ = 90° cut, port 3 at +25° in
φ = 0°), axial ratio < 3 dB at each peak, a maximum (not a null) at each beam, and S21/S31/S32 for isolation.
