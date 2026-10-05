# multiport

Multi-port board: four single-beam sub-apertures and one quad-beam sub-aperture, each with its own feed; 5 lambda_0 nominal sub-apertures with a 28 mm guard gap (about 1 lambda_sw). Per-port patterns use the vector model; directivity is reported in dBi, and the directivity-based aperture efficiency is an UPPER BOUND against the port's own uniform-aperture maximum: it excludes residual surface-wave power, edge diffraction, feed loss and mismatch. Measured boards in the literature reach about 10-15% (Mishra et al., TVT 2026). The isolation flag is a first-order edge-distance check in lambda_sw, not a full-wave result. Unreachable: a cell is unreachable when its closest library patch has a phase error above 30 deg. Those cells still carry the closest library patch in the STEP, and their phase error is included in the patterns. Formulation B uses the normal-incidence reflection-phase library, which is approximate for surface-wave feeding. Port 4 steering limit: requested 60 deg, realised main beam 52 deg, XPD 5.0 dB at the request. Horizontal surface currents near a ground plane limit steering above about 45 deg; this is a model-level limit, and HFSS confirmation is pending. Port 5 beams were equalised with 2 iterations of w <- w * 10^(-level/20); final levels -0.6, -0.0, -0.1, -1.0 dB relative to the port peak, a spread of 1.0 dB (target 2-3 dB).

## Inputs

- field: analytic TM0 monopole, phase-only, per port (no per-port .fld in this app version)
- k_swg: tm0_beta = 224.98 rad/m
- layout: 2x2 small ports (142 mm nominal, 114 mm active after a 28 mm guard gap) plus one large port
- beams: ports 1-4: one beam each (15,0), (30,90), (45,180), (60,270); port 5: all four beams
- polarization: RHCP compensation on every beam (pencil beams); XPD = RHCP over LHCP at each beam
- cells: 8 mm pitch; formulation B with the square-patch library
- scripts: examples/make_submission_figures.py

## Measured

- port 1 requested beams (theta, phi deg): (15.0, 0.0)
- port 1 main beam realised (theta, phi deg): (14.0, 0.0)
- port 1 directivity (dBi): 20.35 (port maximum 23.05 dBi)
- port 1 directivity-based aperture efficiency (upper bound): 0.538
- port 1 unreachable cells: 31.1%
- port 1 cells: 14 x 14
- port 1 XPD, RHCP over LHCP at each requested beam (dB): [32.34]
- port 1 level of each requested beam relative to the port peak (dB): [-0.03]
- port 2 requested beams (theta, phi deg): (30.0, 90.0)
- port 2 main beam realised (theta, phi deg): (28.0, 88.8)
- port 2 directivity (dBi): 19.94 (port maximum 23.05 dBi)
- port 2 directivity-based aperture efficiency (upper bound): 0.489
- port 2 unreachable cells: 32.7%
- port 2 cells: 14 x 14
- port 2 XPD, RHCP over LHCP at each requested beam (dB): [14.42]
- port 2 level of each requested beam relative to the port peak (dB): [-0.16]
- port 3 requested beams (theta, phi deg): (45.0, 180.0)
- port 3 main beam realised (theta, phi deg): (42.0, 182.5)
- port 3 directivity (dBi): 19.01 (port maximum 23.05 dBi)
- port 3 directivity-based aperture efficiency (upper bound): 0.395
- port 3 unreachable cells: 31.1%
- port 3 cells: 14 x 14
- port 3 XPD, RHCP over LHCP at each requested beam (dB): [10.56]
- port 3 level of each requested beam relative to the port peak (dB): [-0.52]
- port 4 requested beams (theta, phi deg): (60.0, 270.0)
- port 4 main beam realised (theta, phi deg): (52.0, 271.2)
- port 4 directivity (dBi): 18.01 (port maximum 23.05 dBi)
- port 4 directivity-based aperture efficiency (upper bound): 0.313
- port 4 unreachable cells: 38.3%
- port 4 cells: 14 x 14
- port 4 XPD, RHCP over LHCP at each requested beam (dB): [4.97]
- port 4 level of each requested beam relative to the port peak (dB): [-1.38]
- port 5 requested beams (theta, phi deg): [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
- port 5 main beam realised (theta, phi deg): (44.0, 182.5)
- port 5 directivity (dBi): 19.04 (port maximum 30.08 dBi)
- port 5 directivity-based aperture efficiency (upper bound): 0.079
- port 5 unreachable cells: 34.2%
- port 5 cells: 31 x 31
- port 5 XPD, RHCP over LHCP at each requested beam (dB): [23.81, 13.68, 10.19, 5.64]
- port 5 level of each requested beam relative to the port peak (dB): [-0.63, -0.02, -0.06, -1.04]
- board size: 540.0 mm x 256.0 mm (guard gap 28 mm)
- isolation flags (first-order, edge gap in lambda_sw): 0 of 10 port pairs flagged
- unreachable threshold (deg): 30.0
- port 5 equalisation history, beam levels (dB) per iteration: [[0.0, -2.01, -5.27, -9.31], [-1.85, -2.59, -2.69, -0.31], [-0.63, -0.02, -0.06, -1.04]]
- port 5 final beam-level spread (dB): 1.02
