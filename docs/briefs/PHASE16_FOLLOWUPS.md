# Phase 16 brief — HMSA follow-ups from the 2026-10-06 live review (original flow only)

Smallest correct changes; no MoM-core changes; fast suite after each item; commit per item; do not move the
`submission-2026-10-06` tag. (The unit-cell/FSS tool idea is parked in `_parked/`; do not implement it.)

1. **Multi-beam saturation in Formulation A (bug).** `formulation_a_control_map` builds
   Ψ = Re{E_ref · conj(Σ_b w_b E_obj,b)} with no normalisation, so with N beams |Ψ| reaches N. Live run: port 5 of the
   2×2+1 board (four beams, 4 mm cells) was 51 % unreachable, clipped to the extreme patch sizes. Normalise the
   summed object wave by its maximum magnitude over the aperture (keep relative weights); report the factor.
   Test: a four-beam single port gives < 5 % unreachable and all four beams within 3 dB.
2. **Stepper navigation (bug).** The `#stepper` links only change the URL hash; they never call `showStep`, so from
   Step 3 the only way back is Home, which wipes the state. Make steps 1–3 clickable once reached, keeping all
   inputs; going back to Step 2 and re-running must reuse them.
3. **Per-port beam cards.** List each requested beam with achieved direction, level, point XPD and integrated XPD
   (port 5 currently shows only its overall peak).
4. **Port side = active aperture.** The guard gap goes outside the stated side (today the preset subtracts the gap,
   so "120 mm" gives 92 mm of hologram). Update preset, UI label and the 1×3 RHCP 25° test (120 mm active, 28 mm gap,
   4 mm cells) back to a normal passing test (< 1.5° error).
5. **Optional pointing correction** (off by default). Small apertures squint toward broadside: the wide array factor
   (≈0.9λ/(L cosθ)) times the falling element pattern (cosθ obliquity, grounded-substrate response) shifts the peak;
   error ∝ 1/L². When enabled, re-aim each beam by its predicted miss, 1–2 iterations; report before/after.
6. **Feed markers in the 26 GHz paper board STEP (bug).** `examples/make_paper_board_26ghz.py::_board_feeds_mm`
   returns board-corner coordinates while `export_multi_port_array` writes a centred board, so the monopole markers
   and the HFSS README positions are offset (≈ half the board). Use the ports' own feed coordinates; regenerate
   `docs/submission/paper26/step`; add a test that every feed marker lies inside its own port rectangle in the
   exported STEP frame (apply the same check to `export_multi_port_array` callers generally).
7. Live-verify Steps 2 and 3 in a browser after 1–4 (the previous round could not run the UI).
