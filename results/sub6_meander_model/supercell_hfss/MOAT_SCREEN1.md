# MOAT SCREEN 1 — ideal ON unit cell

STATUS: STOPPED at the Step 1 anisotropy gate. No OFF or supercell solves were run.

GATES:
- Geometry/readback PASS: OFF 6 copper bodies, 58.769725 mm²; ON 4 bodies, 59.891219 mm². Two E/W sector joins add 1.121494 mm². Ring outer/inner radii 5.25/4.20 mm, island radius 2.75 mm, paddles r 3.20–3.85 mm and ±12°; all inside r 5.65 mm. Radial gaps are 0.45/0.35 mm in OFF. Four ON rows passed design validation, boundary/material readback, and adaptive convergence.
- ON dispersion FAIL: first branch at 5.8 GHz has βx/k0 = 1.07302, βy/k0 = 1.07495. Split across minus along = −0.00193, below the approximately 0.03 stop threshold. Mean 1.07399 is also below 1.08. Both branches are above the light line, so the stop is due to weak anisotropy, not a light-line failure.
- OFF symmetry and supercell gates: NOT RUN because Step 1 failed. The permitted one retune was not used: the explicit stop instruction applies when split is below approximately 0.03.

DATA: HFSS 2025.2, eigenmode, grounded RT5880 3.175 mm, 35 µm PEC copper, 40 mm PEC lid above top copper (z = 40.035 mm), periodic 11.6 mm cell. `across` sets $ph_x; `along` sets $ph_y. β/k0 at 5.8 is from linear interpolation of phase versus the first-mode frequency between 85° and 95°, with free-space phase 80.8° at 5.8 GHz. The first mode is isolated around 5.8 GHz; the next mode is at least 6.091 GHz at 85°.

| Direction | Phase | Mode 1 f GHz | Q | Mode 2 f GHz | Passes | Final tets | Wall s |
|---|---:|---:|---:|---:|---:|---:|---:|
| across x | 85 | 5.717888 | 3160.5 | 6.091276 | 9 | 17530 | 144.7 |
| across x | 95 | 6.200817 | 2097.1 | 6.661028 | 11 | 23839 | 187.3 |
| along y | 85 | 5.708410 | 3068.8 | 6.091282 | 8 | 13112 | 128.4 |
| along y | 95 | 6.201930 | 2106.2 | 6.654972 | 10 | 20746 | 171.1 |

Own prior anchor PEC-lid rows available for comparison: along 85° 5.742396 GHz (Q 3680.8), along 95° 6.129953 GHz (Q 1762.8), across 95° 5.929067 GHz (Q 1594.7). They used a 20 mm lid, so their frequency differences are contextual and are not a matched-lid dispersion comparison. An own anchor PEC across 85° row does not exist.

FLAGS: The first build attempt stopped before solve because the script's analytic paddle-area formula was half the correct sector area. The formula was corrected; both body-count and area readbacks then passed. A 48-logical-core local AEDT session was launched fresh per row and closed afterward; only own PIDs/projects were used. The geometry uses 48 straight segments per 24° paddle sector. No vias, no diode model, no retune, and no switching claim. The moat ON split is about 16× below the 0.03 stop threshold.

FILES: `MOAT_SCREEN1_rows.csv`; four `MOAT_ON_PEC_*_screen1_modes.json` rows and their convergence/profile/build/status files; ON/OFF build-only readback JSONs. Local AEDT projects and temporary files remain in private `hfss_work/` and were not pushed.

HANDOFF: STOPPED | cell: ring 5.25/4.20, island 2.75, four paddles 3.20–3.85 ±12°, E/W joined ON | βx/k0 1.07302, βy/k0 1.07495 | α not measured | anisotropy −0.00193 | loss ratio not measured | HFSS 2025.2 | next request: none; Step 1 stop gate was reached.
