# phase29_twofeed: does the shared-aperture superposition lose SLL? (HFSS check, Phase 29, MODEL_ONLY)

Shared 150 mm-circumradius hexagon (vertex at +/-y), RT/duroid 5880 1.575 mm, 10.54 GHz, 1425 slotted-circle cells (D 6 mm, 6.4 mm lattice), TWO monopoles on the x axis 38.4 mm = 1.35 lambda0 apart
(port 1 at +x, beam tilted 3.65 deg towards +x; port 2 at -x, beam tilted towards -x): two pencil beams one beamwidth (7.3 deg) apart, the geometry of Frey et al.
(IEEE OJAP 2023, 10.1109/OJAP.2023.3261903) Section III. **Control board** `build_board_control_1feed.py`: same hexagon, ONE centre feed, broadside: the single-feed reference of the same cell family.
All numbers are MODEL_ONLY predictions of the OLD model (radial reflection cell); `predictions_new_model.json` (angle-dependent tensor model) is added later. S11 and port coupling: NOT_AVAILABLE_FROM_CURRENT_DATA.

## Files
`build_board.py` (2 ports), `build_board_control_1feed.py` (1 port), `set_weights.py` (CASE 0..3: port 1 alone, port 2 alone, sum, difference), `predictions.json`, `pattern_cuts_old_model*.csv`, `patterns_3deg_old_model*.csv` (RHCP/LHCP gain, hemisphere 3 x 6 deg), `layout`: in the scripts (CELLS).

## Run order
1. AEDT: Tools > Run Script > `build_board.py` (and later the control script in a new project). Messages are printed; projects are saved as `phase29_twofeed.aedt` / `phase29_twofeed_control.aedt`.
2. Validate, Analyze Setup1 (10.54 GHz single frequency; ADD_SWEEP = True at the top adds 10.0-11.1 GHz). One solve gives both excitations: do NOT re-solve per port.
3. Read S11, S22, S21 first (the coupling at 1.35 lambda0).
4. Far field (infinite sphere, theta 0-180, phi -180..180, 1 deg): port 1 alone, port 2 alone (Results > Edit Sources: 1 W on one port, 0 W on the other; or `set_weights.py` CASE), then sum and difference. Export Realized Gain RHCP/LHCP and Axial Ratio.
5. Control board: the same far-field export for its single port.
Compare the SLL of the two-feed beams with the control's: the loss is the point of this check.

## Predictions (old model; RHCP realised gain per 1 W accepted; SLL = highest RHCP lobe outside 1.5 HPBW of the peak, relative)
| Case | Peak gain (dBi) | Peak (theta, phi) | HPBW (deg) | SLL (dB) | AR (dB) |
|---|---|---|---|---|---|
| port 1 alone | 20.6 | (5, 0) | 6.4 | -11.5 | 0.2 |
| port 2 alone | 20.6 | (5, 180) | 6.4 | -11.4 | 0.2 |
| single-feed control | 23.1 | (0, 346) | 7.1 | -24.2 | 0.2 |

Crossover of the two beams (each at the midpoint of the pair, re its own peak): -7.5 / -7.5 dB. Sum (in phase): 17.8 dBi;
difference: 17.5 dBi (see `predictions.json`). What would falsify the model: a two-feed SLL near the control's (the loss is absent in HFSS), or a gain drop far from the 2.5 dB predicted between control and two-feed ports.

<!-- STEP-SECTION-BEGIN -->
## Manual HFSS setup from the STEP files (no AEDT script)

Files (gitignored, copies in `C:\hmsa_sync\data\phase29_twofeed\`): `twofeed_board.step` (2 feeds), `control_board.step` (1 feed). mm. Bodies (names are in the STEP product tree): **Substrate** (solid, 1.575 mm, pin holes r 0.635),
**Ground** (zero-thickness sheet at z = 0, coax opening r 2.05 mm at each feed), **Metal_pattern** (all 1425 / 1430 slotted circles, D 6.0 mm with slot, zero-thickness sheets at z = 1.575), **Pin_1, Pin_2** (solid cylinders r 0.635 mm, z 0 .. 7.575 mm;
top 6.0 mm above the substrate), **Port_1, Port_2** (annulus sheets r 0.635 .. 2.05 mm at z = 0). Zero-thickness sheets are what the AEDT scripts build (no copper thickness). Import: Modeler > Import (STEP), heal off, units mm.
Feed positions (mm): control (0, 0); two-feed port 1 (19.199, 0), port 2 (-19.199, 0).

| Body | Assign |
|---|---|
| Substrate | new material RT/duroid 5880: eps_r 2.2, tan delta 0.0009 |
| Ground, Metal_pattern | Boundaries > Perfect E (sheets) |
| Pin_k | material pec (solve inside off) |
| Port_k | Excitations > Lumped Port, 50 ohm, integration line from (x_k + 0.635, y_k, 0) to (x_k + 2.05, y_k, 0) mm (radially across the annulus, z = 0); Port_1, Port_2 named so |
| Air box (create it) | x, y: board + 15 mm each side, z from -15 mm to 1.575 + 30 mm; Boundaries > Radiation on the box (all faces) |

Solution setup: Driven Modal, 10.54 GHz single frequency, max delta S 0.02, <= 15 passes, lambda refinement on (an interpolating sweep 10.0-11.1 GHz is optional). Far field: Radiation > Insert Far Field Setup > Infinite Sphere, theta 0-180 deg,
phi -180..180 deg, step 1 deg. One solve serves every excitation: set the sources afterwards (Results > Edit Sources: power in W, phase in deg; port 1 alone, port 2 alone, sum = 0.5 W / 0 deg each, difference = 0.5 W / 0 deg and 0.5 W / 180 deg).

Export for each board: the S-matrix (S11, S22, S21; control: S11), and for each excitation the **Realized Gain RHCP and LHCP** (and Axial Ratio) as a function of theta, phi (CSV, 1 deg): control board: its single port; two-feed board: port 1 alone, port 2 alone, sum, difference.
Compare with the prediction tables above (`predictions.json` old model, `predictions_new_model.json` angle-dependent model): the quantity of interest is the SLL (highest RHCP lobe outside 1.5 HPBW of the peak) of the two-feed ports versus the control.
<!-- STEP-SECTION-END -->
