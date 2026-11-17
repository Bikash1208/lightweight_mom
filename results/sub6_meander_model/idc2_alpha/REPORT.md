# Edge-to-edge two-tooth IDC: HFSS-beta inversion (MODEL_ONLY)

Input β/k₀ values were supplied by the user for the 5.8 GHz HFSS cells. No HFSS instance or geometry solve was used. The grounded RT5880 slab is εr=2.2, tanδ=0.0009, h=3.175 mm over PEC; the 0.254 mm layer below the continuous ground is omitted. The separate x/y real TM phase constants are inverted with the same grounded-slab transverse-resonance relation used in `sheet_model.py`. The mean sheet beta is solved again from the mean susceptance before calling the first-order tensor leakage model at θ=30°.

| Cell | β/k₀ x/y input | B_sheet x/y (mS) | B_mean / delta (mS) | β_mean/k₀ | M_t | α (Np/m) | α/k₀ | Ideal converted fraction R150 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MN_R5p6_W0p4_N2_T1p2_D3p0_E | 1.0830/1.0675 | 7.5884/6.4035 | 6.9960/0.5925 | 1.07495 | 0.0880 | 0.0997 | 0.000820 | 2.95% |
| ST_R5p6_L6p72_W3p6 | 1.0990/1.0776 | 8.6679/7.1959 | 7.9319/0.7360 | 1.08792 | 0.1188 | 0.1380 | 0.001136 | 4.06% |

N2_E / ST alpha ratio = **0.723**. N2_E alpha is **0.043×** the 50% target (2.3105 Np/m, alpha/k0≈0.0190). Thus this model does not support 50% conversion at R150. The ideal converted fraction assumes uniform leakage over all 150 mm and no other loss.

Input rounding alone gives N2_E alpha approximately 0.0932–0.1065 Np/m and ST 0.1317–0.1445 Np/m. These are sensitivity ranges, not HFSS uncertainty bounds. Forward inversion reproduces both input beta pairs to <1e-8.

The first-order alpha assumes a rotating-axis tensor modulation. A uniform periodic IDC cell has no such designed radiating harmonic. TM beta alone cannot determine sheet conductance, conductor loss, or a TE eigenmode. The observed flat resonance at p=14.5 mm and r=7.1 mm cannot be predicted or scaled in radius from these two 5.8 GHz beta pairs; that requires frequency-dependent dispersion at multiple radii for this IDC geometry. The earlier seven-leg meander model's branch end near 6.05 GHz is a different cell and is not evidence for the IDC radius trend. No radius dependence was fitted or invented.

## HANDOFF

- STATUS: Priority-one alpha estimate complete from user-supplied HFSS beta pairs.
- Cell JSON: `result.json`; IDC coordinates are user-supplied, no geometry solve.
- β/k₀ along/across: input x/y in table; mean beta solved from mean sheet susceptance.
- α, M_t, loss ratio: table; N2/ST alpha ratio 0.723. Sheet conductance/loss unavailable from beta alone.
- HFSS or MODEL_ONLY: HFSS beta inputs; derived alpha MODEL_ONLY.
- Next request: HFSS beta(f,r) or eigenmode data at two radii to test resonance motion; slot-width/pitch screen after this gate.
