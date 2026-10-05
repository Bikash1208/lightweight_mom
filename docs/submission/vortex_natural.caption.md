# vortex_natural

Normal (natural, no compensation) beam at theta=20 deg, phi=0 deg: the monopole's radial surface current is not compensated, so the aperture radiates a vortex (null on boresight, nonzero OAM winding).

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization natural, oam_charge 0

## Measured

- null_depth_db: -40.103
- winding_number: 1.0
- axial_ratio_db_at_beam: 9.995
- e_r_db_at_beam: -45.796
- e_l_db_at_beam: -40.103
- main_beam_theta_deg: 20.0
- main_beam_phi_deg: 350.137
- directivity_dbi: 25.653
- winding_rhcp_lhcp: [1.0000000000000002, -1.0]
