# Mushroom-type gate cell, baseline HFSS eigenmode (user run, 2026-11-04)

Cell `hfss_exports/sidequest_mushroom_cell/` (patch 5.6 mm, pad 0.8, diode gap 0.4, via r 0.2, RT5880 1.575 mm, pitch 6.4 mm), PEC lid at 30 mm, lid modes removed. Numbers as reported by the user (no raw files received).
- OFF (C = 28.16 fF), TM surface band at 30/60/90/120/150/180 deg: 3.86 / 7.66 / 10.98 / 12.54 / 13.17 / 13.36 GHz; beta/k0 ~ 1.06 at 10.54 GHz. Model (first version): 3.84 / 7.60 / 10.76 / 13.07 / 15.13 / -, beta/k0 1.080: OFF agrees to 2 % up to 90 deg.
- ON (R = 2.44 ohm): lowest band flat at 4.9-5.1 GHz for 90-180 deg (f(180) = 4.92 GHz); second band crosses 10.54 GHz between 60 deg (7.74 GHz) and 90 deg (11.14 GHz), nearly identical to OFF.
  The first model (wire-medium slab, L_eff 2.18 nH) predicted the ON stopband 10.15-10.98 GHz: wrong by about 2.2x in resonance frequency; falsifiers hit (f(180) outside 9.7-10.6; a band crosses 10.54 GHz).
Status: HFSS_CHECKED (one cell, eigenmode, numbers from the user's report).
