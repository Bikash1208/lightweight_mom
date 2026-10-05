# vortex_l2

Charge-2 vortex with the sign rule measured here: for RHCP-compensated beams the realised copol winding equals the requested charge only with l=-2 (measured winding -2, copol null -32 dB). The mirror choice, l=+2 with RHCP, is NOT realised (winding 0, power moves to the cross component); LHCP with l=+2 works instead (winding +2). Cross-polar power at boresight is about -10 dB for both working choices.

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization rhcp, oam_charge -2

## Measured

- null_depth_db: -32.197
- winding_number: -2.0
- axial_ratio_db_at_beam: 1.293
- e_r_db_at_beam: -9.618
- e_l_db_at_beam: -32.197
- main_beam_theta_deg: 16.0
- main_beam_phi_deg: 345.205
- directivity_dbi: 21.314
- copol_component: E_L
- xpd_db_at_beam: -22.579
