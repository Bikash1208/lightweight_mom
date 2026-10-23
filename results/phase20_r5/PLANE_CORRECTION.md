# Reference-plane (de-embedding) correction: one method, applied to every HFSS data set

**Method (the one already used for the run-2 delta table, now written down and shared by all files).** The run-1 plane `d1 = 1.236 mm` is fitted from the bare grounded-slab export against the analytic reference (`fit_reference_plane_mm`, max residual 0.51 deg). For every later export the plane is `d = d1 - shift`, where `shift` is the plane move of that file relative to run 1 *measured from the data*: the median over the TE / TM co-polar entries with |S| > 0.9 at the common (state, f, theta) of phase(b / a) / (2 k0 cos theta) (`run_shift_mm`). A plane move is a pure phase slope in k0, so this does not depend on the cell physics as long as both runs model the same cell, and it needs no number typed in by hand.

| data set | shift vs run 1 (median, 16-84 %) | d used |
|---|---|---|
| run 1 (0.5 GHz steps) | 0 | 1.236 mm |
| run 2 (0.05 GHz steps) | 1.609 mm (1.569..1.688) | -0.373 mm (delta table `hfss_delta_table_run2.npz`) |
| tolerance run, C_off = 0.025 pF rows | 1.557 mm (1.553..1.569) | -0.321 mm (tolerance band) |

**Your note (-0.482 mm for run 2)** comes from the settings change (the port plane moved 1.718 mm against the 1.236 mm needed). The measured shift is 1.609 mm, i.e. 0.109 mm less: the data are the better statement of what the S-parameters contain, which is why the data-driven value is used, and the tolerance file turns out to carry a plane 0.052 mm different from run 2 (shift 1.557 vs 1.609 mm), so each file gets its own value instead of a copy of run 2's.

**Impact:** a plane error is a common phase for the four states (2 k0 d cos theta: 0.109 mm = 2.7 deg at 10.5 GHz, 0.052 mm = 1.3 deg), so it changes no step, amplitude, opposite-handed level or pattern; it enters only through cos theta between the theta = 0 and 30 deg data (< 0.4 deg). Re-running the R5-3 beam table with the tolerance plane changed from -0.374 to -0.321 mm changed no table entry and the sensitivity entries by <= 0.03 dB / 0.4 dB (opposite-handed). If you prefer -0.482 mm for run 2, the delta table would be rebuilt with that value; results would not change beyond these amounts.
