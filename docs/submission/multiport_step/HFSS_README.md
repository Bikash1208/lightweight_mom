# HFSS setup notes for this export (lessons from Phase 12, PHASE12D_ADDENDUM.md)

- Floquet port >= 18 modes, air-box top at 5 mm above the pattern plane, no de-embed (matches this
  project's own MoM reference plane).
- PEC ground and PEC pattern reproduce the MoM reference to ~0.3 deg mean phase off-resonance
  (G1); a 2-mode-port / copper-ground setup was off by up to 4.7 deg -- do not use it.
- Patch thickness: this export is zero-thickness in the solver's own model
  (`ZERO_THICKNESS_METAL_MODEL`); a 35 micron drawn thickness in the CAD file changes the resonant
  phase by ~3.7 deg near resonance (measured, G4) -- expect a corresponding HFSS/MoM offset if you
  keep the drawn thickness.
- Mesh convergence at resonance is first-order in cell size on this project's own MoM (unresolved
  edge-current singularity); do not assume a coarse-mesh MoM comparison point is converged near
  resonance (see `docs/phases/PHASE13_WEBAPP.md`).
- 10-15 GHz, 401 points reproduces this project's own validated frequency sweep.

## Checking the surface-wave dispersion table (Phase 14) against HFSS

If you want a spot check of the eigenmode method's auto-built dispersion table (`docs/phases/
PHASE14_EIGENMODE.md`) against a real HFSS eigenmode solve, for one representative patch size:

- Master/slave periodic boundaries on the unit cell's four side walls, with a linear phase shift
  (not a Floquet port -- this is a source-free eigenmode solve, not an excited one) swept to find
  the frequency where a bound mode exists at the design frequency, or equivalently sweep the phase
  shift at the fixed design frequency and locate the eigenfrequency crossing it.
- Eigenmode solver (not the driven-modal/terminal solver the S-parameter exports above use).
- Air region above the patch tall enough that the bound mode's evanescent tail (decay length
  ~1/alpha, alpha = sqrt(beta^2 - k0^2)) has decayed well before the outer boundary -- a PML or an
  absorbing boundary a few decay lengths up is typically needed; a plain radiation boundary too
  close will perturb the found beta.
- Record beta (or the effective index n_eff = beta/k0) at the same patch size and frequency this
  project's own table used, and compare via the CSV upload path (Step 2 Advanced, "Upload an HFSS
  eigenmode table") or `dispersion.hfss_table.read_dispersion_csv` directly.

## Checking polarization and vortex/OAM (Phase 15) against HFSS

This project's own vector pattern model (`pattern.vector_aperture_pattern`,
`docs/phases/PHASE15_POLARIZATION_VORTEX.md`) predicts polarization and OAM winding from this same
board's own current; to check a real far-field solve against it:

- **Polarization / axial ratio / XPD**: HFSS's Far Field Report supports Ludwig-3 co/cross-pol
  directly ("rE_CoPolar"/"rE_XPolar" report quantities, or "Realized Gain" split by the same), and
  circular components as "rE_RHCP"/"rE_LHCP" (or "RealizedGainRHCP"/"RealizedGainLHCP") -- plot
  either against theta/phi and compare peak direction, XPD (co/cross ratio in dB) and axial ratio
  against this project's own per-beam table (Step 3, or `pattern.vector_aperture_pattern.
  copol_component`/`circular_components` directly).
- **Vortex / OAM**: plot the PHASE (not magnitude) of "rE_RHCP" or "rE_LHCP" on a small circle
  around the beam direction (a theta/phi or u/v cut a few degrees wide, centered on the beam) --
  a flat phase is a pencil beam; N full 360 deg cycles around that circle is an OAM vortex of
  charge N, the same measurement `pattern.vector_aperture_pattern.measure_vortex_winding` makes on
  this project's own predicted pattern. A deep on-axis null in that same RHCP/LHCP magnitude (tens
  of dB below the pattern's own peak) is the other half of the same signature -- "Natural"
  polarization on this monopole-fed board is EXPECTED to show both (the vortex this phase's brief
  is about); a compensated polarization (HP/VP/slant/RHCP/LHCP, Part A) is designed to show neither.

## Feed locations (PHASE15 Part C, multi-feed)

- Feed 0: (-213.000, 71.000) mm
- Feed 1: (-71.000, 71.000) mm
- Feed 2: (-213.000, -71.000) mm
- Feed 3: (-71.000, -71.000) mm
- Feed 4: (142.000, 0.000) mm

## Multi-port isolation (PHASE15 Part C)

Assign one coax (or lumped) port to each monopole feed below, then solve the full 5x5 S-parameter matrix -- the off-diagonal terms S_ij (i != j) are the inter-port isolation this project's own `hologram.ports.isolation_estimate` cannot compute (it only flags tight edge-to-edge DISTANCE between sub-apertures, in surface-wave wavelengths, as a first-order heuristic; a true isolation number needs exactly this full-wave solve).

- Port 1: monopole at (-213.000, 71.000) mm
- Port 2: monopole at (-71.000, 71.000) mm
- Port 3: monopole at (-213.000, -71.000) mm
- Port 4: monopole at (-71.000, -71.000) mm
- Port 5: monopole at (142.000, 0.000) mm
