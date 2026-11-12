# Stage 1 V_OFF feed diagnosis — hard stop

## STATUS

STOPPED. The 5.8 GHz solved V_OFF project was read without changing its geometry. A discrete 5.0–6.5 GHz / 0.25 GHz sweep was started on a private copy and its controller and AEDT process disappeared during the solve. No completed sweep data exist. T2/T5/T1/T8 remain stopped.

## GATES

| Gate | Result |
|---|---|
| Port geometry and vertex readback | PASS |
| Existing 5.8 GHz S/Z readback | PASS |
| Existing-mesh discrete sweep | ABORTED: controller death during solve; no final S11 samples |
| Feed fixed | NO: proposal only; no model edit |

## DATA

| Item | V_OFF project readback |
|---|---|
| PORT sheet | Annulus, z = −3.175 mm, inner radius 0.635 mm, outer radius 2.05 mm; vertices (0.635, 0, −3.175) and (2.05, 0, −3.175) mm |
| Pin | Separate `METAL_PEC` body, radius 0.635 mm, z = −3.175…+8.925 mm; vertex (0.635, 0, −3.175) touches the PORT inner edge |
| Ground | Separate `METAL_PEC_1` body, z = −3.210…−3.175 mm, opening radius 2.05 mm; vertex (2.05, 0, −3.175) touches the PORT outer edge |
| Ground short | No: pin and ground are separate imported objects; at their common z plane the pin ends at r = 0.635 mm and ground starts at r = 2.05 mm, leaving the annular 1.415 mm radial gap occupied by PORT |
| Integration line | (0.635, 0, −3.175) → (2.05, 0, −3.175) mm; 50 Ω port |
| S11 at 5.8 GHz | −0.9535523178 + j0.2956301650; magnitude 0.9983282112, −0.014533 dB |
| Zin at 5.8 GHz | **0.042789 + j7.572952 Ω** |
| Power per 1 W incident | Reflected 0.996659 W; accepted **0.003341 W** (0.3341%) from 1−|S11|² |
| Requested S11 discrete sweep | 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5 GHz; no final samples. Sweep status `ABORTED_CONTROLLER_DIED`. |

The 10.5 GHz EXPORT_STANDARD board spec defines the same annulus and pin: [sub6_final_rf_std.py:5-8,149-161] and the Stage 1 port construction [sub6_planA_anchor_build.py:136-138]. The STEP puts the port at **z = −3.175 mm**, the ground top and RT5880 bottom, across the ground opening. The bias-layer bottom is z = −4.010 mm. The saved project agrees with this port geometry; the failure is a severe 5.8 GHz input mismatch rather than a missing port-edge contact.

The V_OFF far-field grid is 2° in theta. An AEDT request for theta 25° returned the same grid value as theta 24°. Future results should report theta 24° and 26° explicitly, including azimuth 0° and 180°.

## FLAGS

1. **Controller death during sweep:** controller PID and its own AEDT PID were absent when checked while `sweep.json` still said `SOLVING_SWEEP`. The last own license heartbeat was 09:58:22 IST on 2026-11-12. The controller exited without writing a traceback or final S-parameter samples. The requested hard-stop rule applies; no sweep restart or further board solve was launched.
2. **Proposed feed fix, not applied:** retain the validated port annulus/contact geometry; redesign or match the 10.5 GHz probe feed for 50 Ω at 5.8 GHz. The 12.1 mm pin extending from z = −3.175 to +8.925 mm and its coupling to the board should be optimized in a separate feed-only model, then verified by Zin and accepted power before any full-board far-field comparison. No particular new pin length can be justified from the one completed frequency point.
3. The V_OFF four-pass adaptive solution did not meet its 0.02 delta-S target (final 0.083473); the impedance is diagnostic, not a converged final design metric.

## FILES

`inspection.json` contains AEDT vertex/boundary/5.8 GHz solution readback. `sweep.json` records the attempted discrete sweep and abort. `before.prof` is the original adaptive solve profile exported from the private diagnostic copy. No `.py` files were copied or pushed with this report.
