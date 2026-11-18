# A conversion/gain and B-OFF anisotropy sensitivity (MODEL_ONLY)

This update uses the same R150, 11.6 mm lattice, **30 mm** column-stratified stochastic boundary, eight dither seeds, radial ray attenuation, and ideal CP array factor as `REPORT.md`. No cell alpha was requested or calculated; all alpha values below were specified or previously accepted. The 516 radiating cells exclude the centre feed cell.

## A's own conversion and gain

Conversion is assigned to A or B by integrating each ray's power loss through its A/B cells: for a path segment of length ds, dC = P_in[1-exp(-2 alpha ds)]. The table gives A's share of launched surface-wave power, **not the board total**. Ideal A gain is the array-factor directivity toward A plus 10 log10(board radiated conversion), referred to the same launched surface-wave power. It excludes the feed, absorber, finite conductor/dielectric loss and hardware mismatch. Values are medians over eight seeds; the one-beam case keeps B-cell leakage alpha_OFF=0.13 Np/m and a coherent residual B phase.

| alpha_ON (Np/m) | A conversion, two → one beam | A ideal gain, two → one beam | A target-field change |
|---:|---:|---:|---:|
| 1.16 | 14.35% → 14.75% | 13.65 → 13.79 dBi | +0.13 dB |
| **1.33** | **16.08% → 16.60%** | **14.13 → 14.30 dBi** | **+0.15 dB** |
| 1.48 | 17.54% → 18.18% | 14.49 → 14.68 dBi | +0.16 dB |

A gains only about 0.14–0.19 dB of ideal gain because switching B off releases a little more radial surface-wave power for A cells. Board-wide conversion at alpha_ON=1.33 falls from 31.7% to 18.4%; that fall is mostly B's removed radiation and should not be mistaken for a fall in A's own conversion.

## Rerun with OFF alpha = 0.136 Np/m

The same eight 30 mm dither layouts were rerun at the requested alpha_ON values. Suppression compares B-target intensity in the two-beam and one-beam states at **equal launched surface-wave power**. The B-OFF phase is kept coherent and aimed at B, so this is a conservative residual-beam assumption.

| alpha_ON | A conversion, two → one | A ideal gain, two → one | B suppression median (8-seed range) | Two-beam SLL median |
|---:|---:|---:|---:|---:|
| 0.42 | 5.75% → 5.80% | 9.74 → 9.78 dBi | **4.73 dB** (4.72–4.74) | -17.03 dB |
| 1.09 | 13.61% → 13.96% | 13.43 → 13.56 dBi | **8.48 dB** (8.45–8.50) | -16.93 dB |

Neither ON alpha reaches 10 dB B suppression with the assumed coherent OFF anisotropy. If the OFF coherent B coupling amplitude is scaled by m while OFF attenuation stays at 0.136 Np/m, the suppression gains `-20 log10(m)` dB. A 10 dB target would require **m <= 0.545** at alpha_ON=0.42 or **m <= 0.839** at alpha_ON=1.09, relative to the assumed OFF coupling.

One additional **alpha_ON=1.09, OFF anisotropy=0** case sets the coherent B radiation amplitude to exactly zero while retaining 0.136 Np/m OFF attenuation on B paths. In this idealized case the B-target LHCP intensity is exactly zero and suppression is **unbounded**, not a finite predicted dB value. The OFF attenuation's destination is unspecified; it is not counted as coherent B radiation. A's own conversion remains 13.96%, and modeled A gain is about 13.57 dBi. This case shows mathematical headroom only; real OFF cells may retain other scattering and polarization leakage. The model cannot infer OFF anisotropy or phase from alpha alone.

Inputs and outputs: updated `dither_grid.json` for the accepted 1.16–1.48 range; `../dither_grid_newalpha/dither_grid.json` for alpha_ON=0.42/1.09 and alpha_OFF=0.136; `../dither_grid_zero_anis/dither_grid.json` for the zero-anisotropy case. Reproduction: `scripts/research/sub6_sector_dither_grid.py`. No HFSS or new cell solve was run.

## HANDOFF

- STATUS: A's conversion/gain separated; requested B suppression cases complete.
- Cell JSON: same 516-cell dither layouts; no geometry or cell solve.
- beta/k0 along/across: unchanged and not used in this aperture calculation.
- alpha: ON 0.42/1.09; OFF 0.136; accepted anchor range retained for A comparison.
- M_t, loss ratio: OFF coherent anisotropy scale 1 and 0 tested; physical OFF M_t unknown.
- HFSS or MODEL_ONLY: MODEL_ONLY.
- Next request: wait for the capacitance chain, then refine the cell and validate OFF phase/coupling in the crop testbed.
