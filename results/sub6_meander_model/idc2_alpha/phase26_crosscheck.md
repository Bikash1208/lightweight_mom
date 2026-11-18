# Phase 26 leakage-chain cross-check (MODEL_ONLY)

RT5880 (epsilon_r 2.2, tan_delta 0.0009), thickness 1.575 mm; 10.54 GHz. Phase 26 slotted circular cell: pitch 6.4 mm, disc diameter 6.0 mm, slot 4.8 x 1.6 mm, mesh 24/cell. Its x/y TM Bloch beta pair was recalculated with `phase26_ramp_locality.single_cell` and `tensor_sheet.bloch_beta`: beta/k0 = **1.040518276 / 1.105163033**.

The `sub6_idc2_hfss_beta_alpha.estimate` chain gives sheet Bx/By = **5.107206 / 10.827212 mS**, mean/delta = **7.967209 / 2.860003 mS**, beta_mean/k0 = **1.066668426**, and first-order alpha at theta 30 deg = **3.616854 Np/m**. Forward inversion reproduces the beta pair.

The matching Phase 26 eight-cell rotating-slot supercell (period 51.2 mm, K/k0 = 0.555533) reports a complex-root alpha of **3.179144 Np/m** total, or **3.152250 Np/m** after subtracting the uniform-cell dielectric baseline. The first-order/leakage-only ratio is **1.1474** (total ratio **1.1377**): **within x1.3, not x10**. Beta_mean/k0 - K/k0 = 0.511135 corresponds to a radiated angle of about 30.74 deg, close to the 30 deg first-order evaluation.

This is a same-family, same-mesh MODEL_ONLY check of the chain at one point. It does not validate the IDC cell's actual leakage, account for conductor loss, or make this complex-root solver an HFSS result.

Sources: `results/phase26/supercell_table_4p8x1p6.json` (n_cells=8), `scripts/research/phase26_ramp_locality.py`, `scripts/research/sub6_idc2_hfss_beta_alpha.py`.

## HANDOFF

- STATUS: Cross-check complete.
- Cell JSON: Phase 26 D6.0/L4.8/W1.6, n=8 rotating-slot row.
- beta/k0 along/across: x 1.040518276, y 1.105163033.
- alpha: first-order 3.616854 Np/m; complex-root leakage 3.152250 Np/m; ratio 1.1474.
- M_t: 0.39993; loss ratio: root total/leakage 1.0085.
- HFSS or MODEL_ONLY: MODEL_ONLY.
- Next request: compare IDC leakage against an IDC-specific complex-root or measured result.
