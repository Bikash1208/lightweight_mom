# phase29_sevenport (MODEL_ONLY, HFSS-checked margin)

Same method and settings as `hfss_exports/phase29_twofeed/`: 150 mm hexagon (270 x 310 mm board), 1.575 mm RT5880, 6 mm circles, pins r 0.635 mm to z = 7.575 mm, one fixed port per tile (P1 centre, P2..P7 ring), 10.54 GHz single frequency, Driven Modal, one solve then Edit Sources per port.

- `sevenport_board.step` (gitignored; copy in `C:\hmsa_sync\data\phase29_sevenport\`): bodies Substrate, Ground, Metal_pattern, Pin_1..7, Port_1..7. Manual setup: the STEP section of `../phase29_twofeed/README.md` (feeds/ports 1..7 instead of 2).
- `build_board.py`: AEDT script, same geometry.
- `ports_and_excitations.csv`: port positions, design beam, predicted gain / peak direction / SLL raw and calibrated (gain -2.5 dB, SLL +3.5 dB, from the two-feed HFSS), and the excitation per case (port k alone 1 W).
- `predictions.json`: same plus the all-in-phase case. Inter-port coupling not modelled.
