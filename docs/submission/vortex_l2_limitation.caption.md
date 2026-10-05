# vortex_l2_limitation

LIMITATION, not a clean charge-2 vortex: a requested l=+2 (RHCP) puts the boresight power into the opposite (cross) circular component with zero winding, so the charge is not realised as requested in this model. Reported as measured.

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization rhcp, oam_charge 2

## Measured

- null_depth_db: -22.13
- winding_number: -0.0
- axial_ratio_db_at_beam: 1.366
- e_r_db_at_beam: -0.027
- e_l_db_at_beam: -22.13
- main_beam_theta_deg: 20.0
- main_beam_phi_deg: 0.0
- directivity_dbi: 27.926
- copol_component: E_L
- xpd_db_at_beam: -22.103
