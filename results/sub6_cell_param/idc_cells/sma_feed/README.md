# SMA feed for every design board (standing rule, 2026-11-21)

Files: `sma.step` (user-supplied SMA connector, SHA-256 4c8c69271dc4715d7a7f2834e35c62da8a69d28475ceb9b062841713d974c01f), `full_board_sma_export.py` (board STEP, `--feed sma`), `full_board_sma_hfss.py` (HFSS build/solve, any machine via `--aedt-version`), `BOTH_oldfeed_251_1pass.aedtz` (machine-check project, see below).

Place both scripts in `scripts/research/` of a checkout that also has `scripts/research/sub6_step_exact.py` and `src/transmitarray/export.py` (from compute/sub6-meander-model; S4 already has them), and use the v31 layout CSV (`results/sub6_meander_model/full_board_gate/v31/layout_R150_th30_v31.csv` on that branch).

    python full_board_sma_export.py layout_v31.csv step_out --state BOTH --feed sma      # CadQuery environment
    python full_board_sma_hfss.py --step step_out/BOARD_R150_GATE_V31_BOTH_SMA.step --sma sma.step --out run_dir --state BOTH --aedt-version 2025.1 --cores 8 [--pin-tip 8.925] [--max-passes 1 --min-passes 1] [--sweep]

Geometry (SMA frame, mm): pin r 0.3048 protruding 7.794 above the flange face, coax PTFE r 0.635..2.0574, outer conductor r 2.667, flange 10.16 x 10.16 x 1.676, wave port on a coax annulus r 0.635..2.0574 at z -7.62 built by the script (the STEP's Circle2 disc r 2.6926 is deleted: HFSS 2026.1 rejects a port edge not backed by PEC), with a PEC cap disc behind it. Placed with the flange face on the board ground bottom (z -3.21; the board bias layer already has a 12 x 12 mm window); pin tip is `--pin-tip` (board z, default 8.925 = 12.1 mm above the ground like the old probe); ground anti-pad r 0.9017; air box extends 25 mm below the board.

Machine check (old lumped-port feed): `BOTH_oldfeed_251_1pass.aedtz` is the fresh BOTH project built on the lab-B PC in HFSS 2025.1 (one adaptive pass, results stripped, no history purge; SHA-256 adbecbb255112e73a70db8e484411a234ca01c737c260415fb7f2c19f73fe6f8). On that PC it reads S11 0.0890 - j0.0540 (about -19.7 dB), Z11 59.4 - j6.5 ohm. Solving this exact file on another machine (min = max passes = 1) tests the machine, not the model.

Design decisions (user, 2026-11-21): the bias laminate has a 12 x 12 mm window under the SMA, so the flange sits directly on the ground plane (no via fence in the model). The user runs the full-board solves; every script here must run unchanged on their machine (build natively on the machine that solves; do not move projects between HFSS versions: a migrated project made the TAU mesher fail on S4 and fall back to a coarse Classic mesh). Each machine keeps its own converged SMA BOTH as its reference; the old probe-feed references (-16.76 dB, 9.03/9.00 dBi) are history only. One pin tip for all machines (frozen from the converged SMA BOTH on the lab-B PC).
