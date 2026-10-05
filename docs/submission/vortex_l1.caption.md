# vortex_l1

An intentionally requested OAM vortex, charge l=+1 (RHCP request). It reproduces the natural beam: the uncompensated monopole current already carries winding +1, so the natural beam is an l=+1 vortex.

## Inputs

- field: inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)
- feed (x0, y0) mm: (0.0, 0.0)
- cells: 8 mm pitch, 37 x 37 reference grid
- substrate: RT/Duroid 5880, eps_r 2.2, tan_delta 0.0009, 1.575 mm, PEC backing
- method: Formulation B (phase-difference hologram) with the square-patch library (fast mesh)
- scripts: examples/make_submission_figures.py
- beam: theta 20 deg, phi 0 deg, desired_polarization rhcp, oam_charge 1

## Measured

- null_depth_db: -40.103
- winding_number: 1.0
- axial_ratio_db_at_beam: 9.995
- e_r_db_at_beam: -45.796
- e_l_db_at_beam: -40.103
- main_beam_theta_deg: 20.0
- main_beam_phi_deg: 350.137
- directivity_dbi: 25.653
- copol_component: E_L
- xpd_db_at_beam: 5.692
