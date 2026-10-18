# Phase 20 R2 — PIN diode in the loaded MoM (labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET)

## Diode values (item 1) — `src/ris/diodes.py`
Skyworks SMP1345 Series datasheet 200046U (2018-06-19), SMP1345-079LF (SC-79). Each value carries its provenance in the code.

| Quantity | Value | Where | Status |
|---|---|---|---|
| R_s ON | 1.5 ohm typ (I_F = 10 mA), 3.5 ohm typ (1 mA); F = 100 MHz | Table 3, p. 2 | PROVISIONAL_DATASHEET |
| C_T OFF | 0.19 pF typ (1 V, 1 MHz); C_T at 0 V not given, may be higher | Table 3, p. 2 | PROVISIONAL_DATASHEET |
| L_s | 0.7 nH (SC-79) | Table 1, p. 1 | PROVISIONAL_DATASHEET |
| R_p (OFF) | not in the datasheet | — | NOT_AVAILABLE_FROM_CURRENT_DATA |
| R_p sweep 10 k / 3 k ohm | sensitivity only | — | ASSUMED_SENSITIVITY_NOT_DEVICE_DATA |

The datasheet covers 10 MHz–6 GHz; every value here is used at 8–12 GHz, an **extrapolation** with no measured de-embedding. The OFF
state resonates (L_s, C_T) at 1/(2 pi sqrt(L_s C_T)) = 13.8 GHz, within 2 GHz of the band, so OFF is sensitive to C_T and L_s.

## HFSS representation of a series diode circuit (item 3)
A PIN diode is a series circuit (ON: R_s + L_s; OFF: C_T [|| R_p] + L_s). HFSS "Lumped RLC" boundaries do not give an unambiguous
series combination: the one source found says only that one boundary can represent "any combination" and that HFSS computes from the
entered R, L, C (Simutech, *How to Model Lumped Components in Ansys HFSS*,
<https://simutechgroup.com/resources/blog/how-to-model-lumped-components-in-ansys-hfss/>); it does not say whether R, L, C act in
parallel or series, and I have not verified this in HFSS. Representation that does not depend on it:

**A chain of single-type boundaries joined by metal islands** (`hfss_chain`): `R_s` | island | `C_T` (or `C_T` with `R_p` set in the SAME
boundary) | island | `L_s`. Each boundary carries ONE element type (pure R, pure L, pure C), except the OFF boundary with C and R_p,
which relies on the parallel reading. Geometry used in the MoM check: 3 boundaries/elements of 0.25 mm along x, islands 0.5 mm, strip width 0.5 mm.

**Your one-line check before use:** assign one boundary with C and R both set across a 1 mm gap on a short line and confirm its impedance
is 1/(1/R + jwC) (parallel). If HFSS combines them differently, R_p cannot be added as its own series element: run OFF lossless
(R_p = infinity, the nominal case) in HFSS and treat the R_p sweep as MoM-only.

### MoM check: HFSS chain vs the ideal series model (one load of total Z_L over the same 1 mm span), 10 GHz
`scripts/research/phase20_r2_hfss_representation.py` → `hfss_representation.{json,md}` (meshes 40 and 80, theta = 0, 30, TE and TM).

| Strip metal each side of the 1 mm span | ON: max d phase, max d abs(S11) | OFF: max d phase, max d abs(S11) |
|---|---|---|
| none (elements touch the 6 mm arms) | 0.21 deg, 9e-4 | **67 deg, 0.18** |
| 0.5 mm of 0.5 mm strip (**pads**) | 0.13 deg, 6e-4 | **6.0 deg, 0.011** (worst case R_p = 3 kΩ; lossless OFF: 5.9 deg, 0.000) |

* **ON agrees to 0.2°** (a 3.5 / 1.5 ohm resistor and 0.7 nH inductor split into two boundaries is the same as one load).
* **OFF differs by ~6°, not a mesh effect** (mesh 40 vs 80 change it by 0.4°): where the capacitor and the inductor sit along the
  span matters because the OFF circuit (−40 ohm net) is a difference of two large opposite reactances (−84 + 44). Treat ~6° as the
  representation uncertainty of OFF in this cell; the HFSS-vs-ideal comparison should not be expected to do better.
* **Pads are required.** With the boundaries touching the wide arms the 2-D current crowding makes the OFF chain differs by up to 67° from
  the ideal model (ON is insensitive); with 0.5 mm of strip on each side it is 6°; with 1 mm it is 1.4° (separate probe).
