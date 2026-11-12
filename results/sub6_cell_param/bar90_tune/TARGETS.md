T1–T5 definitions and ranking. Save this verbatim as `results/sub6_cell_param/bar90_tune/TARGETS.md` and push it.
All targets are evaluated on the hologram branch at 5.8 GHz.

* Hologram branch: a mode below the light line with SW slope (> 0.02 GHz/deg), not LL-HUG (β/k₀ ≥ 1.002), and not GEOM_INVARIANT, tracked across phases by frequency continuity.
* If two candidates remain in a direction, label them (a, b), evaluate T1–T4 for each OFF/ON pairing, and report all pairings.
* β/k₀ at 5.8 GHz comes from linear interpolation in (f, β) between the two phase points that bracket 5.8 GHz.
   * If 5.8 GHz isn’t bracketed, extrapolate from the nearest two points and flag it `EXTRAP`.
   * If the branch tops out below 5.8 GHz, mark it `NO_BRANCH_AT_5p8`. That direction fails T3.

Targets (pass thresholds):

* T1 (OFF anisotropy): |across − along| ≥ 0.14 in OFF. Report the sign. A negative sign only means the slot axis is rotated 90° in the hologram mapping, so it doesn’t fail the test.
* T2 (ON isotropy): |across − along| ≤ 0.02 in ON.
* T3 (bound waves): β/k₀ ≥ 1.00 in both states and both directions. ON along between 0.99 and 1.00 passes, flagged `T3_FLAG`.
* T4 (no nearby resonance): no RES (slope < 0.005 GHz/deg) or HYBRID (0.005–0.02 GHz/deg) mode between 5.5 and 6.1 GHz, in either state or direction. The hologram branch itself and GEOM_INVARIANT modes are excluded. Modes outside 4.0–6.5 GHz are ignored.
* T5 (report only, not a pass/fail): for the hologram branch near 5.8 GHz, OFF and ON:
   * Q;
   * slope in GHz/deg;
   * implied attenuation α ≈ ω / (2 Q v_g), with v_g taken from the slope, in dB/m.

A config passes if T1, T2, T3 and T4 all pass.
Ranking:

1. Number of T1–T4 passed (more is better).
2. Among full passes:
   * a. smaller ON |across − along| (T2 margin);
   * b. larger OFF |across − along| (T1 margin);
   * c. lower ON attenuation (T5);
   * d. lower OFF attenuation (T5);
   * e. larger minimum via-pad margin.
3. Among non-passes: smallest T1 shortfall, then smallest T2 excess.

Report a ranking table with every config and all values. Don’t drop configs that fail.
