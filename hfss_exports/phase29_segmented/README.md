# phase29_segmented (option A, MODEL_ONLY)

150 mm hexagon, 10.54 GHz, same cell/stack-up/settings as `phase29_twofeed` and `phase29_sevenport`: 7 equal-area regions (central hexagon + 6 sectors), one port per region (P1 centre, broadside; P2..P7 sector axes 0, 60, ... 300 deg, beam theta 30 deg). Cells of a region come from the single-feed synthesis; target strength tapered (cosine) over 12.8 mm inside each region boundary.

- `segmented_board.step` (gitignored; copy in `C:/hmsa_sync/data/phase29_segmented/`), `build_board.py` (AEDT), `ports_and_excitations.csv` (port table; predictions raw and calibrated: gain -2.5 dB, SLL +3.5 dB from the two-feed HFSS).
- `steering_lut.csv`: per direction (theta 0..45 step 3, phi step 15): best single port, switched and blended gain, and the 7-port MRC weights (power W and phase deg, sum of power = 1 W) for Results > Edit Sources. Ports taken as uncoupled.
- Manual STEP setup: STEP section of `../phase29_twofeed/README.md` (ports 1..7).
