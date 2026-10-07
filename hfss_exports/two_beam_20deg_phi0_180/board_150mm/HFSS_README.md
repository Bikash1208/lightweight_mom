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
  resonance (see `docs/history/PHASE13_WEBAPP.md`).
- 10-15 GHz, 401 points reproduces this project's own validated frequency sweep.

## Checking the surface-wave dispersion table (Phase 14) against HFSS

If you want a spot check of the eigenmode method's auto-built dispersion table (`docs/history/
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
