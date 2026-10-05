# multiport

Multi-port board: four single-beam sub-apertures and one quad-beam sub-aperture, each with its own feed, with a 28 mm guard gap (about 1 lambda_sw) between sub-apertures. Per-port patterns use the vector model; directivity is reported in dBi, and the directivity-based aperture efficiency is an UPPER BOUND against the port's own uniform-aperture maximum: it excludes residual surface-wave power, edge diffraction, feed loss and mismatch. Measured boards in the literature reach about 10-15% (Mishra et al., TVT 2026). The isolation flag is a first-order edge-distance check in lambda_sw, not a full-wave result.

## Inputs

- field: analytic TM0 monopole, phase-only, per port (no per-port .fld in this app version)
- k_swg: tm0_beta = 224.98 rad/m
- layout: 2x2 small ports (88 mm nominal, 60 mm active after a 28 mm guard gap) plus one large port
- beams: ports 1-4: one beam each (15,0), (30,90), (45,180), (60,270); port 5: all four beams
- polarization: RHCP compensation on every beam (pencil beams)
- cells: 8 mm pitch; formulation B with the square-patch library
- scripts: examples/make_submission_figures.py

## Measured

- port 1 requested beams (theta, phi deg): (15.0, 0.0)
- port 1 main beam realised (theta, phi deg): (13.0, 355.1)
- port 1 directivity (dBi): 15.52 (port maximum 17.48 dBi)
- port 1 directivity-based aperture efficiency (upper bound): 0.637
- port 1 unreachable cells: 26.5%
- port 2 requested beams (theta, phi deg): (30.0, 90.0)
- port 2 main beam realised (theta, phi deg): (24.0, 93.7)
- port 2 directivity (dBi): 15.04 (port maximum 17.48 dBi)
- port 2 directivity-based aperture efficiency (upper bound): 0.571
- port 2 unreachable cells: 36.7%
- port 3 requested beams (theta, phi deg): (45.0, 180.0)
- port 3 main beam realised (theta, phi deg): (37.0, 177.5)
- port 3 directivity (dBi): 14.55 (port maximum 17.48 dBi)
- port 3 directivity-based aperture efficiency (upper bound): 0.510
- port 3 unreachable cells: 34.7%
- port 4 requested beams (theta, phi deg): (60.0, 270.0)
- port 4 main beam realised (theta, phi deg): (43.0, 271.2)
- port 4 directivity (dBi): 13.98 (port maximum 17.48 dBi)
- port 4 directivity-based aperture efficiency (upper bound): 0.447
- port 4 unreachable cells: 38.8%
- port 5 requested beams (theta, phi deg): [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
- port 5 main beam realised (theta, phi deg): (13.0, 0.0)
- port 5 directivity (dBi): 18.06 (port maximum 25.32 dBi)
- port 5 directivity-based aperture efficiency (upper bound): 0.188
- port 5 unreachable cells: 34.3%
- board size: 324.0 mm x 148.0 mm (guard gap 28 mm)
- isolation flags (first-order, edge gap in lambda_sw): 0 of 10 port pairs flagged
