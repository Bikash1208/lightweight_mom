# Nylon cell spot-check set, 28 GHz (optional)

Six single cells spanning the library's 360-degree window (pitch 5.0 mm, total thickness 25.0 mm, 1.0 mm base
plate, flats on the x/y axes).  One STEP per cell, z = 0 is the bottom (exit) face, z = 25.0 the top (feed-side) face.
`rcwa_reference.csv` has the model's complex transmission (project convention exp(+jwt)) at the sampled angles.

HFSS setup: nylon eps_r 2.45, tan_d 0.0185 (single material); unit-cell boundaries (master/slave, 5 mm square) on the
four side faces, Floquet ports above and below, scan angle theta/phi from the CSV.  S21 phase must be referenced to
the **top face of the cell for the incident wave and the bottom face for the transmitted wave** (port de-embedding of
the air above/below, no extra offset).  `t_pp` is p (in-plane E) polarisation, `t_ss` is s; at theta = 0, p = x.

Targets: phase within +-3.6 deg and |S21| within +-0.1 dB of the CSV.  HFSS comparison status:
NOT_AVAILABLE_FROM_CURRENT_DATA (model-only until you run it).

| cell | kind | radius [mm] | RCWA phase [deg] | RCWA |S21| [dB] (theta=0, p) |
|---|---|---|---|---|
| 1 | pillar | 0.000 | -142.0 | -0.64 |
| 2 | pillar | 1.239 | -213.4 | -0.84 |
| 3 | pillar | 1.772 | -294.1 | -1.11 |
| 4 | hole | 2.153 | -366.7 | -1.24 |
| 5 | hole | 1.848 | -430.9 | -1.90 |
| 6 | hole | 1.315 | -511.3 | -1.79 |
