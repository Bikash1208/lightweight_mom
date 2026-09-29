"""Part 3: surface-wave dispersion for HMSA (holographic metasurface antenna) design.

The brief's own physics point (PHASE13_WEBAPP_PROMPT.md "HMSA synthesis"): the per-cell control
quantity for a surface-wave hologram is the unit cell's surface reactance / surface-wave index
(dispersion) -- NOT normal-incidence S11 phase, which is what every earlier phase's library
actually stores. Two sources, in the order the brief asks for:

  1. `hfss_table.py` -- import of an HFSS eigenmode dispersion table (geometry parameter vs beta or
     X_s at the design frequency; CSV) as a library type (`labels.DISPERSION_FROM_HFSS`). Usable now.
  2. `mom_scaffold.py` -- a MoM surface-wave-eigenmode capability. NOT implemented (a documented,
     tested scaffold only, per the brief's explicit "if it can't be finished cleanly in this phase,
     ship (1) and leave (2) as a documented, tested scaffold -- don't block the app on it").

`analytic_tm0.py` is the one piece of real, closed-form physics this package DOES ship: the TM0
surface-wave dispersion of the bare grounded dielectric slab (no patches), used as the sanity-check
baseline the brief asks the field-import page to compare an imported reference field's own
estimated beta(r) against.
"""
