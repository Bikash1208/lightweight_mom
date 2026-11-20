# Corrected B-ALL-DIODE-V2 — OFF-D full-board HFSS 2026 R1

**STATUS:** SOLVED and extracted, one full-board OFF-D state at 5.68 GHz. The source was a fresh copy of the same validated V2 build used for ON-A; Zr=2.667ohm, Zx=+48.667ohm. No geometry changes or other state solves were made in this step.

**GATES:** Setup1 Direct, min 6/max 12 adaptive passes, max ΔS 0.02. Normal completion and convergence in six passes. Delta-S per pass: pass 1 n/a; pass 2 0.14626; pass 3 0.02259; pass 4 0.00587; pass 5 0.00411; pass 6 **0.00162**. Requested/granted cores 45/45; maximum reported tets 2,004,442; peak HFSS memory 49.8 GB; wall time 1438 s.

**DATA:** S11 **−17.591 dB**; accepted power **98.259%**. RHCP peak **8.655 dBi** at (32°,0°), first-null SLL **−8.849 dB**. LHCP peak **1.657 dBi** at (36°,178°), first-null SLL **−3.997 dB**. LHCP at (36°,180°) **1.602 dBi**. Compared with ON-A, the global LHCP peak decreases by **4.706 dB**; these two impedance states alone do not establish 10 dB suppression.

**FILES:** `OFF_D_solve.prof`, `OFF_D_solve_status.json`, `OFF_D_extract/metrics.json`, updated `metrics_summary.json`, and `OFF_D_extract/farfield/GATE_B_ALL_DIODE_V2/element_pattern/exportfield_0.ffd`. The solved AEDT project remains local in `hfss_work/full_board_gate_v31_diode_v2/solve_OFF_D_001/`.
