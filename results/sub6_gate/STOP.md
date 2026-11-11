# SUB6 gate test: HARD STOP (wave port touching lattice faces)

Driven-modal strip model built and validated (design validation PASSED, both wave ports and both lattice pairs
accepted at assignment time). At SOLVE time HFSS 2025.1 (v251) refused the excitation:

    [error] Port refinement, process hf3d error: Matching boundary LatticeY_air_Master1 touches non-Floquet port P1.
    Please contact Ansys technical support.

Per instructions no other excitation was substituted. No rows were solved; results/sub6_gate/rows is empty.
Setup at the stop: REF_H20, wave ports P1/P2 (YZ sheets, ground top -> lid, full 11.6 mm cell width, 1 mode,
integration line vertical), lattice pairs on the y faces (phase 0deg; one pair on AIR faces, one on RT5880 faces).
Code: hfss_exports/sub6_gate/ (gate_spec.py, gate_hfss.py, test_gate_spec.py; tests never launch AEDT).
