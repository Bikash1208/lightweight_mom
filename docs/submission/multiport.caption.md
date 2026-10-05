# multiport

Multi-port board: four single-beam sub-apertures and one quad-beam sub-aperture, each with its own feed. Per-port patterns use the vector model; gain is reported as directivity and as aperture efficiency against the port's own uniform-aperture maximum. The isolation flag is a first-order edge-distance check in lambda_sw, not a full-wave result.

## Inputs

- field: analytic TM0 monopole, phase-only, per port (no per-port .fld in this app version)
- k_swg: tm0_beta = 224.98 rad/m
- layout: 2x2 small ports (60 mm) plus one large port, gap 0 mm
- beams: ports 1-4: one beam each (15,0), (30,90), (45,180), (60,270); port 5: all four beams
- polarization: RHCP compensation on every beam (pencil beams)
- cells: 8 mm pitch; formulation B with the square-patch library
- scripts: examples/make_submission_figures.py

## Measured

- port 1 requested beams (theta, phi deg): (15.0, 0.0)
- port 1 main beam realised (theta, phi deg): (14.0, 4.9)
- port 1 directivity (dBi) / aperture efficiency: 15.70 / 0.664 (port maximum 17.48 dBi, aperture 36 cm^2)
- port 2 requested beams (theta, phi deg): (30.0, 90.0)
- port 2 main beam realised (theta, phi deg): (24.0, 88.8)
- port 2 directivity (dBi) / aperture efficiency: 15.13 / 0.583 (port maximum 17.48 dBi, aperture 36 cm^2)
- port 3 requested beams (theta, phi deg): (45.0, 180.0)
- port 3 main beam realised (theta, phi deg): (36.0, 182.5)
- port 3 directivity (dBi) / aperture efficiency: 14.57 / 0.512 (port maximum 17.48 dBi, aperture 36 cm^2)
- port 4 requested beams (theta, phi deg): (60.0, 270.0)
- port 4 main beam realised (theta, phi deg): (45.0, 266.3)
- port 4 directivity (dBi) / aperture efficiency: 13.76 / 0.425 (port maximum 17.48 dBi, aperture 36 cm^2)
- port 5 requested beams (theta, phi deg): [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
- port 5 main beam realised (theta, phi deg): (15.0, 0.0)
- port 5 directivity (dBi) / aperture efficiency: 16.42 / 0.196 (port maximum 23.50 dBi, aperture 144 cm^2)
- isolation flags (first-order, edge gap in lambda_sw): 8 of 10 port pairs flagged
