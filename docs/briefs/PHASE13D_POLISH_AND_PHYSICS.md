# Phase 13D brief — findings from the live walk-through of the 3-step app (2026-09-29, 22:25)

Walk-through: light mode, 800 px pane. Step 1 with a monopole field in the user's exact format (151×151, 2 mm, NaN
at feed, β = 224.96 rad/m) at 10.54 GHz → Step 2: 1 beam θ=20°, φ=0°, RT/duroid 5880 1.575 mm, cell 4 mm,
defaults otherwise → Step 3 run → export. Screenshots: `docs/screenshots/` (step 1 also in dark mode).

## What works
Step 1 phase map + the one-line check ("β = 225 rad/m · matches the bare substrate at 10.54 GHz ✓"); the stepper,
cards, segmented beam-count control, preset substrate list; synthesis resampled to 74×74 cells; the main beam lands
at θ = 20.0°, φ = 0.0° as requested; export produced a STEP with every cell; More tools page is tidy.

## 1. Physics / results that are wrong or misleading (fix first)
1. **The phase-lookup method can't realise this hologram at small cells, and the app doesn't say why.** With 4 mm
   cells the achievable phase span is ≈ 33° (layout colour bar −22°…+11°), so the layout is effectively binary:
   **4006 of 5476 cells can't reach their value, mean/max error 73°/162°.** The beam still forms (binary holograms
   work) but efficiency/side lobes are poor. This is the known limitation: a surface-wave hologram's control
   variable is surface impedance, not normal-incidence reflection phase. Required:
   - In Step 3, when > 20 % of cells are unreachable, show one plain sentence explaining the cause and the fix
     ("Small cells give only ~33° of phase range with the phase-lookup method. Upload an HFSS dispersion table to use
     the impedance method, or increase the cell size (fewer cells per fringe).").
   - Show the achievable range in Step 2 as soon as board + cell size + frequency are set (the auto-library knows it).
   - **Implement the MoM surface-wave dispersion capability** (currently a scaffold): for a grounded periodic cell,
     find β(L) by locating the TM surface-wave pole with real kx > k0 in the existing Floquet framework; validate
     against the analytic TM0 of the bare slab (β = 224.96 rad/m at 10.54 GHz for this board — also measured from the
     user's HFSS monopole field) and, where available, against an HFSS eigenmode table. Then make the impedance method
     the default in Step 3 when dispersion data exist (MoM or HFSS). This is the most valuable physics item left.
   - The auto-library ("fixed 6-point mesh", 6 sizes?) — sample sizes densely enough to resolve the curve and use an
     aligned mesh; report its phase span and resolution.
2. **Side-lobe level −1.5 dB is wrong.** The φ = 0 cut shows the highest side lobe ≈ −10 dB; the −1.5 dB comes from
   the bright ring at the horizon in the u–v map (≈ −5…−8 dB near θ → 90°). The aperture-integration model lacks the
   element/obliquity factor (e.g. cos θ for a tangential aperture field, or the proper surface-wave-antenna element
   pattern); add it (documented), and compute SLL excluding the main lobe region correctly. Re-check that the
   end-fire rise (−8.6 dB at 90° in the cut) disappears or is physically justified.
3. **"Done." appears before export finishes.** Export of 5476 cells took 126 s and produced a **103 MB STEP** with no
   progress indication; the card appeared only afterwards. Show progress ("Building CAD… 40 %"), estimated time, and
   don't say Done until everything is. Make the STEP optional or offer a lighter format first (DXF / Gerber of the
   patch layer — a PCB fab wants Gerber, not a 100 MB STEP); consider a single-solid or tiled STEP for large arrays.
4. **Zip contents**: brief asked for pattern PNGs, layout CSV + DXF, STEP, HFSS notes, project record. Only STEP, CSV
   and notes are listed; add the PNGs, DXF (or Gerber) and the project record (reloadable JSON).

5. **Patch-layout PNG (user request).** Step 3 and the export zip must include a true-geometry drawing of the board:
   every patch drawn at its actual L × W and position on the substrate (not a phase/colour map) — **patches bright
   copper (#E8893A, thin #C8651B edge), substrate light grey (#D9DADC, #B8B9BC outline), white page** — substrate
   outline, feed point marked, mm axes, a title with cell count / pitch / beams, and a legend or histogram of the patch
   sizes used. High-resolution (≥ 2000 px) PNG, plus the same as SVG/PDF for print. Reference rendering made from this
   walk-through's `layout.csv`: `docs/screenshots/patch_layout_example.png` (74 × 74, 4 mm pitch; it also makes the
   §1.1 problem visible — only 6 sizes, almost all 0.6 mm or 3.6 mm, i.e. a near-binary hologram).

## 2. Wording and polish
- Brand: replace "lwmom" in header/title with a product name (e.g. "HMSA Designer"); keep `lwmom` as the package.
- Replace every " -- " in UI text with an em dash or rewrite ("phase -- magnitude is optional" → "phase; magnitude
  optional"; "Aperture size -- from your uploaded file, [mm]" → "Aperture size" + hint "From your file").
- Units: "θ [deg]" → "θ (°)"; "[GHz]" → "GHz" suffix inside the field; plot titles "Cut at phi=0deg" → "φ = 0° cut",
  axis "theta (deg)" → "θ (°)"; chips "162deg" → "162°" (the list already uses °).
- Step 1 copy: "This is the surface wave your board will hologram against." → "We use it as the reference wave."
- Magnitude file uses the browser's grey "Choose File" button; the dispersion upload too. Use the same drop-zone /
  link style as the phase file.
- The TM0 line says "the bare substrate" before a substrate was chosen — name it ("RT/duroid 5880, 1.575 mm") and
  re-check when the board changes in Step 2.
- Step 2: hide mean reactance / modulation depth when the phase-lookup method is in use; add the "Accuracy
  (fast / accurate)" control the brief asked for; default cell size should follow the fringe rule (the app already
  computes cells-per-fringe — suggest a value and show the achievable phase span next to it).
- Step 3: before running, show a one-line summary of the design and a time estimate; button "Run synthesis" →
  "Design antenna". "Experimental." chip needs its reason ("Phase-lookup method is experimental for surface-wave
  holograms").
- Stepper sits tight under the header — add ~16–24 px top margin.
- More tools descriptions still read technical ("S11/S21", "MoM", "Touchstone (.sNp)") — fine for that page, but
  shorten ("Check one cell's reflection").

## 3. Done when
Items in §1 fixed with regression tests (the SLL case and the unreachable-explanation case use this walk-through's
inputs); §2 applied; fresh screenshots (light + dark) replace the ones in `docs/screenshots/`; fast suite green;
commit per section.
