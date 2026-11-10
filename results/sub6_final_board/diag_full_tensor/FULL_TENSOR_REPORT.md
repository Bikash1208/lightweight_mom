# PARALLEL TEST: full Bloch tensors vs reflection law. MODEL_ONLY / NOT_HFSS_VALIDATED. Report only.
Method: bloch sheet tensor (same MoM path/mesh as the scans: static mesh 64, island mesh 96, forder 128), kt = beta k0 of the table (per state), direction theta' = 0/22.5/45/67.5/90 deg to the slot axis (cell unrotated, kt rotated; lattice orientation not varied). 10 static rows used on the board + island ISL_W2.9_q0.8_Li0.5_c0.6 OFF/ON. Symmetry check (island OFF theta'=-45): Yuu, Yvv even, Yuv odd, confirmed to 1e-12.
Data: tensors.json, board_full_tensor.json. Scripts: diag_full_tensor_lib/run/board.py.
Board: J = E0 e^{-jP} R(axis) Y_c(theta') e_rho, E0 = a/(|ImDelta|/2) (model's own alpha, W, P, axes, flips); no sign choice.
Results: see chat RESULT (tables derived from board_full_tensor.json).
