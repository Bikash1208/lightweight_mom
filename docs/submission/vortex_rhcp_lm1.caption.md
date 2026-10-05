# vortex_rhcp_lm1

Clean RHCP vortex, charge l=-1. With a monopole feed the sign of l must be opposite to the handedness in this code's convention (RHCP: l<0; LHCP: l>0). Measured: RHCP winding -1 with an on-axis RHCP null; LHCP (cross) leakage stays well below the null.

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization rhcp, oam_charge -1

## Measured

- null_depth_db: -31.421
- winding_number: -1.0
- axial_ratio_db_at_beam: 12.758
- LHCP_db_at_beam: -35.492
- RHCP_db_at_beam: -31.421
- main_beam_theta_deg: 16.0
- main_beam_phi_deg: 0.0
- directivity_dbi: 23.398
- copol: RHCP
- xpd_db_at_beam: 4.072
