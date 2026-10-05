# pencil_rhcp

The same beam with per-beam RHCP compensation: a pencil beam. The RHCP component is the main curve, and the LHCP (cross) component is the one shown dashed.

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization rhcp, oam_charge 0

## Measured

- null_depth_db: -0.03
- winding_number: 0.0
- axial_ratio_db_at_beam: 1.442
- LHCP_db_at_beam: -21.67
- RHCP_db_at_beam: -0.03
- main_beam_theta_deg: 20.0
- main_beam_phi_deg: 0.0
- directivity_dbi: 27.924
- copol: RHCP
- xpd_db_at_beam: 21.64
