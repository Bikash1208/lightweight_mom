"""Part 3: surface-wave dispersion for HMSA (holographic metasurface antenna) design.

The brief's own physics point (PHASE13_WEBAPP_PROMPT.md "HMSA synthesis"): the per-cell control
quantity for a surface-wave hologram is the unit cell's surface reactance / surface-wave index
(dispersion) -- NOT normal-incidence S11 phase, which is what every earlier phase's library
actually stores. Two sources:

  1. `hfss_table.py` -- import of an HFSS eigenmode dispersion table (geometry parameter vs beta or
     X_s at the design frequency; CSV) as a library type (`labels.DISPERSION_FROM_HFSS`).
  2. `mom_scaffold.py` -- Part 3D: a real MoM surface-wave-eigenmode search (`find_surface_wave_mode`,
     `build_dispersion_table_from_mom`), labelled `DISPERSION_FROM_MOM` + `EXPERIMENTAL`. Built in
     this project's own Z(kx, ky) framework (no new solver core code); validated against
     `analytic_tm0` for a lightly-loaded cell. Prefer an HFSS table when you have one -- this is
     younger and self-validated on fewer boards.

`analytic_tm0.py` is the one piece of real, closed-form physics this package DOES ship: the TM0
surface-wave dispersion of the bare grounded dielectric slab (no patches), used as the sanity-check
baseline the brief asks the field-import page to compare an imported reference field's own
estimated beta(r) against.
"""
