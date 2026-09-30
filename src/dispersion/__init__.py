"""Part 3: surface-wave dispersion for HMSA (holographic metasurface antenna) design.

The brief's own physics point (PHASE13_WEBAPP_PROMPT.md "HMSA synthesis"): the per-cell control
quantity for a surface-wave hologram is the unit cell's surface reactance / surface-wave index
(dispersion) -- NOT normal-incidence S11 phase, which is what every earlier phase's library
actually stores. Three sources:

  1. `hfss_table.py` -- import of an HFSS eigenmode dispersion table (geometry parameter vs beta or
     X_s at the design frequency; CSV) as a library type (`labels.DISPERSION_FROM_HFSS`).
  2. `sheet_model.py` -- Phase 14 Part A: this project's OWN computed dispersion table, built from a
     homogenized patch-sheet-admittance model (Luukkonen et al. 2008) combined with the closed-form
     grounded-slab transverse-resonance condition -- `labels.DISPERSION_FROM_SHEET_MODEL` +
     `EXPERIMENTAL`. Validated against `analytic_tm0` (small patch -> bare-slab beta) and checked
     for a monotonic beta(L). This is the one Formulation A actually builds a table from when no
     HFSS table is uploaded.
  3. `mom_scaffold.py` -- a direct MoM surface-wave-eigenmode search over the periodic array's own
     Z(kx, ky), kept as a documented research item: it robustly finds the BARE-slab TM0 pole, but
     Phase 14's own investigation found it cannot resolve a loading-DEPENDENT mode at any
     interactively-affordable mesh size (see that module's own docstring for the finding and why).
     Not used to build Formulation A's dispersion table.

`analytic_tm0.py` is the one piece of real, closed-form physics this package DOES ship: the TM0
surface-wave dispersion of the bare grounded dielectric slab (no patches), used as the sanity-check
baseline the brief asks the field-import page to compare an imported reference field's own
estimated beta(r) against.
"""
