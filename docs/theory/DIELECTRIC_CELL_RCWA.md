# Dielectric unit-cell RCWA (Phase 19, `src/dielectric_cell/`)

Stack: air | z-invariant layers (octagonal nylon pillar in air; optional uniform base plate) | air.
Rectangular lattice, any incidence (theta, phi). Independent of the periodic-MoM core.

## Convention
Project convention `exp(+jwt)` (CONVENTIONS.md). Internally the solver uses `exp(-jwt)`
(`eps = eps_r(1 + j tan_d)`, fields `exp(+j(kx x + ky y + kz z))`) and conjugates its outputs.

## Formulation
With `z' = k0 z`, `H~ = -j eta0 H`, harmonics `Kx = kx/k0`, `Ky`:
`d/dz'[Ex;Ey] = A [Hx;Hy]`, `d/dz'[Hx;Hy] = B [Ex;Ey]`,
`A = [[-Kx e^-1 Ky, Kx e^-1 Kx - I],[I - Ky e^-1 Ky, Ky e^-1 Kx]]`,
`B = [[-KxKy - E_yx, Kx^2 - E_yy],[E_xx - Ky^2, E_xy + KyKx]]`,
`e = [eps]` (Laurent, E_z is tangential to pillar walls), `E` = in-plane effective permittivity.
Layer modes: `A B W = -gamma^2 W`, `V = j B W Gamma^-1` (downward-propagating).
Interface matching of a layer embedded in air gives `S11 = (A - XBA^-1XB)^-1 (XBA^-1XA - B)`,
`S21 = A^-1 X (A + B S11)` with `A = W^-1 W0 + V^-1 V0`, `B = W^-1 W0 - V^-1 V0`, `X = exp(j gamma k0 d)`;
layers combine by Redheffer star product (stable for thick, lossy layers).

## Factorisation
`laurent`: `E = [eps] (x) I`. `normal` (default): `E = [eps] P_t + [1/eps]^-1 P_n`, `P_n = n n^T`,
`n` = gradient of the distance to the pillar boundary (periodic images), `P_t = I - P_n`
(Schuster et al., 2007). 1-D limit reduces to Li's rule.

## Reference planes and polarisation
S11 at the top face of the first layer, S21 at the bottom face of the last; no air de-embedding
(an empty cell is an air layer: `S21 = exp(-j k0 cos(theta) h)`). s = E perpendicular to the plane of
incidence, p = in-plane; transmission uses identical unit vectors, reflection the tangential-E sign
(`r_pp = r_ss` at normal incidence).

## Validation (tests/test_dielectric_cell.py)
- Uniform slab vs independent Macleod characteristic-matrix formula, TE/TM, 0/25/40 deg, lossy and
  lossless: agreement ~1e-14 (asserted 1e-9).
- Power balance lossless octagon: 1e-14 in the phi=0 plane, <1e-8 off-plane (NVM operator is not exactly
  Hermitian). C4 TE/TM equality ~1e-10; mirror, 90-deg azimuth rotation, reciprocity (|t_sp|=|t_ps|, 6e-5 rel).
- Wiener bounds for a sub-wavelength pillar.
- Convergence (25 mm nylon pillar, 5 mm pitch, 28 GHz): normal-vector method phase moves <0.1 deg between
  M=6 and M=12 (M = harmonics per axis, N = (2M+1)^2); the Laurent rule is still ~2 deg off at M=12.
  Production order **M = 8**; the slow test asserts M=8 vs M=14 within 0.2 deg and 0.01 dB.
