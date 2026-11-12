# Feed fix 1 (centre probe, 5.8 GHz) - AEDT 2025.1, feed-only model

**STATUS:** F1 and F2 pass on the feed-only model. The feed block (`feed_block.json`) is ready for the new exporter.

## GATES
| Gate | Result |
|---|---|
| F1 no DC contact pin-ground | PASS: pin r 0.635 mm, ground opening r 2.05 mm: radial gap 1.415 mm (numeric, from the model read-back). On the board STEP the first non-pin conductor in the ground plane is the ground itself at r = 2.05 mm, and the first other top-layer conductor is at r = 6.0 mm (5.365 mm from the pin surface). |
| F2 Zin not near 0 / open; accepted power >= 50 % | PASS: Zin(5.8 GHz) = 37.98 - j7.40 Ohm, S11 = -15.93 dB, accepted power 97.4 % |

## DATA (feed-only model, 5.0-6.5 GHz, 0.1 GHz discrete steps, see `feed_only_sweep.csv`)
- Stack: RT5880 3.175 mm / ground 35 um (opening r 2.05) / RT5880 bias layer 0.254 mm; pin r 0.635, z -3.175 ... +8.925; 60 mm diameter disc; radiation boundary on a box 15 mm beyond the disc and the pin; no cells.
- Port: Lumped Port 50 Ohm on the annulus r 0.635-2.05 at z -3.175, integration line (0.635,0,-3.175) -> (2.05,0,-3.175), no de-embedding.
- S11 is below -10 dB from 5.2 GHz; at 5.8 GHz -15.9 dB; best -27.4 dB at 6.4 GHz. Accepted power 85 % at 5.0 GHz, 97.4 % at 5.8 GHz.
- Solve: 81.6 s for the 16-point sweep (8 cores). Design validation (`validate_simple`) passed.

## What the earlier short-like result was NOT
The earlier V_OFF result Zin = 0.043 + j7.57 Ohm is not caused by a pin-ground contact: the pin and the ground are separate bodies with a 1.415 mm radial gap, nothing else is within 5.3 mm of the pin on the board STEP, and the same pin/port block gives 38 - j7 Ohm in the feed-only model. The collapse therefore comes from the loaded board (cells, vias, FR4/bias stack) or from that project's solution (4 passes, delta S 0.083 against a 0.02 target), not from the feed block. Final matching must be done on the full board.

## FLAGS
1. The saved V_OFF project itself is not on this server (it lives under `C:\lwmom_labB\stage1_run3_work` on the earlier run's machine); the pin / annulus / opening values come from `results/sub6_planA/stage1/feed_diag/inspection.json` and `RESULT.md` and were re-checked numerically on `BOARD_N240_RF_full_win.step` imported into AEDT 2025.1 (1244 objects, `board_step_feed_clearance.json`). The clearance probes have 0.01 mm resolution.
2. The feed-only substrate needs a pin-shaped hole; without it the validation says "Parts PIN and SUBSTRATE intersect" and the solver fails.
3. `pyaedt` `validate_full_design` calls `Desktop.GetMessages` internally. Earlier controller builds in this session used it (the run rule says never call GetMessages). It is replaced by `validate_simple` from this point on; the feed-only runs after the first failed attempt use `validate_simple`. The first feed-only attempt (`feed_only_base`) and the builds of the already-solved Round 1 / screen configs used the old call.
4. "Top launcher": the board STEP has only the pin cylinder above the substrate (r 0.635, up to z = 8.925); no separate launcher body. The feed block reproduces exactly that.
5. The 12 x 12 mm SMA window is not part of the feed-only model; it is in `feed_block.json` for the board exporter.

## FILES
`feed_block.json`, `feed_only_result.json`, `feed_only_sweep.csv`, `board_step_feed_clearance.json`, `RESULT.md`.
