# Existing island-cell calibration method

Located in the repository:

- scripts/research/sub6_final_cal.py, calibrate(state): f(phi) fitting,
  beta/k0 at 5.8 GHz, group-slope/Q loss conversion.
- results/sub6_final_board/hfss_data.json: the actual user-returned
  ISL_vias2 eigenfrequencies/Q and mode-identification notes.
- jobs/sub6-board-island-hfsscal.md: original calibration request.
- Commit 21fe444, hfss_exports/sub6_island_choke/HFSS_SETUP.md:
  archived ISL_vias2 setup, +90-degree frame, 0.8-mm FR4, two bias stubs,
  30-mm lower air extension and radiation bottom.
- Commit b98be34, hfss_exports/sub6_island_cell/HFSS_SETUP.md:
  original ISL_novia reference package.

The recorded alpha_x/beta calibration data were ISL_vias2, not a no-via
measurement. The archived no-via package used a solid ground and a PEC
bottom at the ground bottom; the via package's distant radiation bottom
was introduced to remove the artificial bottom-cavity effect. The new
packages preserve that radiation-bottom method, with the user's revised
0.6-mm anti-pads and specified 5.8-GHz Impedance sheet values.

The method is a PERIODIC EIGENMODE dispersion/Q extraction, not a driven
transmission-line/S21 attenuation extraction:

beta(phi)=phi*pi/(180*p), p=0.0116 m.
v_g=d(2*pi*f)/d(beta)=2*pi*p*df/d(phi radians).
alpha_x=omega/(2*Q*v_g)=f*(d beta/d f)/(2*Q), Np/m.

Follow the same top-surface TM mode across phase, obtain phi at f=5.8 GHz,
and use that branch's Q and local dispersion slope. Record f, Q, mode index
and Ez identification per phase/state/direction. Exclude the low-Q via
resonance, isotropic FR4 underside mode and air modes recorded in hfss_data.

Historical fitting details: OFF had two HFSS points per direction; ON
used a single HFSS 95-degree anchor and retained the MoM branch slope
because the 80-degree HFSS points were ambiguous. OFF 95-degree Q was
missing, so 85-degree Q was reused. These are historical limitations;
they must not silently carry over to the new OFF/ON data. The new manual
requests >=3 nearby HFSS points to obtain the local slope without relying
on the now mesh-unreliable ON MoM branch. This is the same eigenmode/Q
measurement method with a directly measured slope, not a new HFSS solver.

Quantity distinction: alpha_x is uniform-cell modal attenuation from Q.
The earlier retune table's alpha is the first-order tensor hologram
radiation leakage, and was explicitly retained as MoM-only by sub6_final_cal.
The two alphas are not interchangeable. Perturbative directional no-via
loss estimates in model_reference.csv are additional labelled MODEL_ONLY
sheet estimates: alpha_loss=-Re(Y_total)/(d Im(Y_total)/d beta), derived
by linearizing Y_total(beta-j alpha)=0. All roots and extrapolations are
reported; no negative/strong-loss estimate is presented as a qualified
forward-wave alpha_x. The weak-loss screening threshold is alpha/beta<0.01;
a negative estimate can indicate a backward branch or an interpolation
artifact and is retained as a raw value with an unreliable-reference flag.

Both historical diode values and boundary representation differ from the
new requested Impedance sheets; do not transfer old calibrated numbers
to these packages as new predictions or new HFSS results.
