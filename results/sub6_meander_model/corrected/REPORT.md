# Corrected seven-vertical-leg meander study (MODEL_ONLY)

The user's baseline is seven slot legs at x = -3.6, -2.4, -1.2, 0, 1.2, 2.4, 3.6 mm, each spanning y = -1.5 to +1.5 mm, 0.4 mm wide, with alternating end bridges. The disc radius is 5.6 mm; period 11.6 mm. The continuous cut clears the rim and its interior contains no copper. The substrate is grounded RT5880, 3.175 mm. `along_x` means propagation along global +x (across the leg pitch); `across_y` is global +y (parallel to the legs). The old horizontal-run screen is superseded.

The earlier β used a TM-like hybrid tensor root, not a pure TM mode. Below are separate projected TM and TE roots of the bound-mode transverse-resonance equations. Off-diagonal tensor terms couple them in the physical mode; see each row's coupling ratio and hybrid character. A blank TE entry means no bound TE root in the searched interval.

| Geometry | TM β/k₀ x | TM β/k₀ y | TE β/k₀ x | TE β/k₀ y | M_t principal | α rad potential (Np/m) | Ideal converted fraction R150 | Formal M_t for 50% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline_7leg_w0p4 | 1.0666 | 1.2950 | — | — | 3.887 | 2.368 | 50.9% | 3.839 |
| 7leg_w0p6 | 1.0666 | 1.2747 | — | — | 2.460 | 2.408 | 51.4% | 2.409 |
| 9leg_w0p4 | 1.0661 | 1.0267 | — | — | 0.284 | 0.603 | 16.6% | 0.555 |
| 5leg_w0p4 | 1.0657 | 1.0781 | 1.0384 | — | 0.071 | 0.068 | 2.0% | 0.414 |

Converted fraction uses 1-exp(-2 α R) with R=0.150 m and the first-order *rotating-axis modulation potential* α. The required α for 50% is ln(2)/(2R)=2.3105 Np/m. The M_t numbers scale tensor anisotropy at fixed mean sheet and β; they are local mathematical thresholds, not design guarantees. M_t can exceed 1 when principal reactances have opposite signs. No finite-conductor loss, feed loss, spatially varying coupling, sector boundary, or HFSS validation is included.

The R150 two-sector ideal array-factor worst local maximum is RHCP at θ≈2°, φ=180° and its symmetric LHCP partner at φ=0°, -13.62 dB. With ρ^-1/2 exp(-αρ) illumination, α from the corrected baseline, a 5.8 mm feed-floor radius and a 10 mm tanh sector transition, the worst local maximum is -13.32 dB near boresight; directivity is 18.84 dBi per beam, AR≈0.90 dB and XPD≈25.7 dB. The soft boundary is one explicit modeling choice, not a fabricated feed design. Full array-factor data are in `aperture_r150_rerun.json`.

## TM/TE dispersion

| Frequency (GHz) | TM β x (rad/m) | TM β y (rad/m) | TE β x (rad/m) | TE β y (rad/m) |
|---:|---:|---:|---:|---:|
| 5.30 | 116.21 | 122.44 | 163.75 | — |
| 5.60 | 124.13 | 140.03 | — | — |
| 5.80 | 129.66 | 157.42 | — | — |
| 6.00 | 135.51 | 180.82 | — | 127.57 |

Finite-difference v_g/c by frequency interval is recorded in `analysis.json`. A flattening claim requires a continuous tracked TM branch and small positive v_g/c; missing or multiple roots are flagged, not interpolated through.

## HANDOFF

- STATUS: Corrected geometry and 5.8 GHz screen complete; dispersion/mesh status in result files.
- Cell JSON: `finalists.json`.
- β/k₀ along/across: projected TM and TE table above; hybrid values in row JSON.
- α, M_t, loss ratio: α and M_t above; loss comparison remains first-order MODEL_ONLY.
- HFSS or MODEL_ONLY: MODEL_ONLY; no switching claim.
- Next request: HFSS cell dispersion and board-level CP check after mesh and TM-branch review.
