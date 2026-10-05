# multiport

Multi-port board: four single-beam sub-apertures and one quad-beam sub-aperture, each with its own feed; 5 lambda_0 nominal sub-apertures with a 28 mm guard gap (about 1 lambda_sw). Per-port patterns use the vector model; directivity is reported in dBi, and the directivity-based aperture efficiency is an UPPER BOUND against the port's own uniform-aperture maximum: it excludes residual surface-wave power, edge diffraction, feed loss and mismatch. Measured boards in the literature reach about 10-15% (Mishra et al., TVT 2026). The isolation flag is a first-order edge-distance check in lambda_sw, not a full-wave result.

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
- port 1 XPD, RHCP over LHCP at each requested beam (dB): [np.float64(32.34)]
- port 1 level of each requested beam relative to the port peak (dB): [np.float64(-0.03)]
- port 2 requested beams (theta, phi deg): (30.0, 90.0)
- port 2 main beam realised (theta, phi deg): (28.0, 88.8)
- port 2 directivity (dBi): 19.94 (port maximum 23.05 dBi)
- port 2 directivity-based aperture efficiency (upper bound): 0.489
- port 2 unreachable cells: 32.7%
- port 2 cells: 14 x 14
- port 2 XPD, RHCP over LHCP at each requested beam (dB): [np.float64(14.42)]
- port 2 level of each requested beam relative to the port peak (dB): [np.float64(-0.16)]
- port 3 requested beams (theta, phi deg): (45.0, 180.0)
- port 3 main beam realised (theta, phi deg): (42.0, 182.5)
- port 3 directivity (dBi): 19.01 (port maximum 23.05 dBi)
- port 3 directivity-based aperture efficiency (upper bound): 0.395
- port 3 unreachable cells: 31.1%
- port 3 cells: 14 x 14
- port 3 XPD, RHCP over LHCP at each requested beam (dB): [np.float64(10.56)]
- port 3 level of each requested beam relative to the port peak (dB): [np.float64(-0.52)]
- port 4 requested beams (theta, phi deg): (60.0, 270.0)
- port 4 main beam realised (theta, phi deg): (52.0, 271.2)
- port 4 directivity (dBi): 18.01 (port maximum 23.05 dBi)
- port 4 directivity-based aperture efficiency (upper bound): 0.313
- port 4 unreachable cells: 38.3%
- port 4 cells: 14 x 14
- port 4 XPD, RHCP over LHCP at each requested beam (dB): [np.float64(4.97)]
- port 4 level of each requested beam relative to the port peak (dB): [np.float64(-1.38)]
- port 5 requested beams (theta, phi deg): [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
- port 5 main beam realised (theta, phi deg): (15.0, 0.0)
- port 5 directivity (dBi): 22.59 (port maximum 30.08 dBi)
- port 5 directivity-based aperture efficiency (upper bound): 0.179
- port 5 unreachable cells: 33.0%
- port 5 cells: 31 x 31
- port 5 XPD, RHCP over LHCP at each requested beam (dB): [np.float64(27.27), np.float64(19.3), np.float64(10.73), np.float64(5.33)]
- port 5 level of each requested beam relative to the port peak (dB): [np.float64(0.0), np.float64(-2.01), np.float64(-5.27), np.float64(-9.31)]
- board size: 540.0 mm x 256.0 mm (guard gap 28 mm)
- isolation flags (first-order, edge gap in lambda_sw): 0 of 10 port pairs flagged
