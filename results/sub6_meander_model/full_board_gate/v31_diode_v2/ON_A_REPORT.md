# Corrected B-ALL-DIODE-V2 — ON-A full-board HFSS 2026 R1

**STATUS:** SOLVED and extracted, one full-board state at 5.68 GHz. The source was a copy of the validated V2 build; Zr=3.333ohm and Zx=−232.667ohm. Setup1 was Direct, min 6/max 12 adaptive passes, max ΔS 0.02. No other geometry or state was solved in this step.

**GATES:** Normal completion and adaptive convergence in six passes. Delta-S per pass: pass 1 n/a; pass 2 0.18802; pass 3 0.02226; pass 4 0.00490; pass 5 0.00207; pass 6 **0.00170**. Requested/granted cores 45/45; maximum reported tets 2,004,442; peak HFSS memory 54.8 GB; wall time 1552 s.

**DATA:** S11 **−17.263 dB**; accepted power **98.122%**. RHCP peak **8.449 dBi** at (30°,0°), first-null SLL **−9.191 dB**. LHCP peak **6.363 dBi** at (26°,180°), first-null SLL **−5.858 dB**. LHCP at (36°,180°) **−1.836 dBi**. The FFD uses the same 2° theta/phi grid and circular-polarization convention as the stored v3.1 BOTH results. The residual LHCP peak has moved from the BOTH direction (36°,180°) to (26°,180°).

**FILES:** `ON_A_solve.prof` is the raw HFSS profile; `ON_A_solve_status.json` has setup/readback and per-pass convergence; `ON_A_extract/metrics.json`, `metrics_summary.json` and `ON_A_extract/farfield/GATE_B_ALL_DIODE_V2/element_pattern/exportfield_0.ffd` contain the extracted gains and field. The solved AEDT project remains local in `hfss_work/full_board_gate_v31_diode_v2/solve_ON_A_001/`.
