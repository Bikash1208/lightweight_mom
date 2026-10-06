# Phase 16 brief — Unit-cell plane-wave S-parameter tool, with Periodic / Radiation boundary options

## Why
On 2026-10-06 a user STEP (square-loop FSS, 10 mm period, 8.8/6.6 mm loop, RT5880 1.0 mm, no ground) was solved
with the existing periodic MoM **unchanged**: a duck-typed `SquareLoop` geometry (`in_feature`, `corners`,
`polarity`, `feature_area`), `compute_s_parameters(unit_cell, f, MeshSettings(100,100,100), ReferencePlanes(0,0,h))`,
backing "air". Result: S21 null 9.64 GHz, −10 dB stop band 7.62–11.68 GHz, |S11|²+|S21|² = 0.999, mesh
0.2→0.1 mm shift +0.2 %, Floquet 60→100 shift 0.003 GHz. Data: `hfss_exports/fss_square_loop_10mm/`.
Everything below productises that workflow. No change to the MoM core.

## Part A — geometry
1. New duck-typed shapes in `src/geometry.py` (same interface as `AnnularSlotPatch`): square loop (outer, inner),
   circular ring (R_out, R_in), cross / Jerusalem cross (arm length, width, end-cap), split ring (gap), and
   "polygon" (list of outer + hole polygons, point-in-polygon `in_feature`). Rotation + centre offset for all.
   Existing rectangle/patch/dipole/slot/annular-slot stay as they are.
2. **STEP / DXF import:** read the metal pattern (cadquery for STEP: take the solid(s) whose bbox is thin in z, or the
   top face of the thinnest body; ezdxf for DXF), detect the period from the repeated copies (or the user enters it),
   take ONE copy, and build a "polygon" geometry. Report: period, metal area, thickness found (warn that metal is
   modelled zero-thickness PEC), substrate thickness if a dielectric body is present. Never trust a material name in
   the file (the user's file said PLA) — the user picks the substrate.
3. Preview PNG of the cell (copper on grey) with the mesh grid overlaid.

## Part B — mesh that fits the geometry
1. **Exact-grid helper:** collect every edge coordinate of the pattern relative to the cell boundary (x and y),
   find the largest dx = period/n that hits all of them within 1 µm (gcd search over n ≤ 200). Report it
   ("loop edges need 0.1 mm → 100×100"). If n is impractical, offer the nearest n and report the worst edge error.
2. Cost estimate before running (unknowns, ~seconds per frequency from a 1-frequency timing run) and a
   memory guard.

## Part C — the run
1. Inputs: substrate (preset or εr/tanδ/h), backing (air → S11,S21,S12,S22; PEC → S11 only), frequency sweep
   (start/stop/step + optional dense band), polarisation x/y (both computed anyway), reference planes (default: metal
   plane and substrate faces; optional offsets to match HFSS de-embedding).
2. Background job with progress, resumable checkpoint per frequency (`runner.py` pattern), parallel over
   frequencies when cores allow.
3. **Adaptive refinement:** after the coarse sweep, find nulls/peaks/phase jumps (|S21| or |S11| local minima, ±180°
   flips) and add points there automatically until the null frequency is located to < 0.1 % (linear complex
   interpolation of S across the bracket, as done by hand on 2026-10-06).
4. **Convergence report (automatic, small):** (a) Floquet order 0.6× vs 1× at the key frequency; (b) mesh: the
   chosen exact grid vs the next coarser exact grid of a geometry-compatible proxy if needed, Richardson first-order
   estimate of the remaining shift. Print "null 9.64 GHz, converged estimate 9.66 ± 0.02 GHz".
5. Sanity: energy check |S11|²+|S21|² ≤ 1 (lossless → 1) per frequency, flag violations; warn when higher Floquet
   orders propagate (period > λ0) or substrate-trapped orders approach (period·√εr ≈ λ0).
6. Oblique incidence: NOT in this phase — keep the "normal incidence only" message; leave a TODO hook.

## Part D — outputs
CSV (dB, deg, re, im; S11, S21, S12, S22, xx/yy/xy/yx), **Touchstone .s2p/.s4p** (per polarisation, or 4-port
x/y), magnitude + phase PNG (dataviz colours; null/peak annotated, −10 dB band marked), a summary card (null/peak
frequencies, −10 dB band, fractional bandwidth), and an **HFSS setup note**: unit cell with Master/Slave walls 0°,
Floquet ports top/bottom de-embedded to the stated planes, mode count, air-box height, and the warning that a finite
N×N board with radiation boundaries is NOT the same thing.

## Part E — app
"Unit cell S-parameters" under More tools: Step 1 geometry (shape form or STEP/DXF upload + preview),
Step 2 boundary (Periodic / Radiation-estimate / Radiation-HFSS export, with N×M or panel size, and the two measurement planes: size and distance), substrate/backing/sweep/accuracy (exact-grid suggestion, cost estimate), Step 3 results (plot, table,
downloads, convergence card). Plain language, same look as the HMSA pages.

## Part F — boundary condition option (Periodic vs Radiation)
A user-facing "Boundary" choice. Be explicit in the UI about what each one can and cannot give.

1. **Periodic — infinite array (default, exact).** Floquet / Master–Slave equivalent, as in Parts A–E. Outputs true
   S11/S21/S12/S22 at Floquet ports. Normal incidence.
2. **Radiation — finite panel, plane-defined S-parameters (estimate).** The user gives N×M cells (or the panel size
   from the STEP: e.g. 11×11, 110 mm), the incident plane wave (normal incidence, x or y), and two **measurement
   planes**: size (default: the central 80 % of the panel; options 50–150 %) and distance d (default λ0/2 at the
   centre frequency; user-editable) in front of and behind the panel.
   - Currents: take the periodic solution's cell current at each frequency and place it on every cell of the panel
     (local periodicity, uniform illumination; infinite substrate). State this in the UI.
   - Fields on the planes: FFT the panel current to its plane-wave spectrum, apply the layered-medium spectral
     response (`green_function.spectral_dyad`, the same response the MoM uses), propagate to z = ±d with
     exp(−j k_z d), inverse FFT. Zero-pad the FFT enough that the plane fields are not aliased (pad ≥ 4× panel).
   - **Plane S-parameters:**
     S21_plane = (1/A) ∫_A E_total · ê* dA / E0 · exp(+j k0 d_back)   (total field behind = incident transmitted
     through the unpatterned slab + scattered field from the panel currents),
     S11_plane = (1/A) ∫_A E_scat · ê* dA / E0 · exp(+j k0 d_front)   (scattered field only in front).
     The exp(+j k0 d) factors move the phase reference to the substrate faces so the numbers line up with the
     Floquet S11/S21. Also give the cross-pol projections.
   - Report sensitivity: S11/S21_plane vs plane size (50/80/100/150 %) and vs d (λ/4, λ/2, λ, 2λ) at the key
     frequencies, as a small table, and always plot the periodic (Floquet) S11/S21 on the same axes.
   - Also output monostatic RCS (dBsm, with the PEC flat-plate reference 4πA²/λ²) and the bistatic pattern, from the
     layered far-field code.
   - Label: "finite-panel estimate — edge diffraction of the currents is not modelled (uniform cell currents),
     finite substrate edges not modelled; values depend on plane size/distance; approaches the periodic result as
     the panel grows."
3. **Radiation — finite panel, full-wave (HFSS).** Generate the HFSS model instead of solving it: STEP of the whole
   panel (pattern + finite substrate N·p × M·p × h), plus a setup note: air box ≥ λ/4 at the lowest frequency
   (or PML), Radiation/PML on all faces, incident plane wave (θ, φ, E-pol), far-field setup for RCS, and the **same
   two measurement planes** as rectangles (sizes and distances written out), with the Field Calculator recipe for
   S11_plane / S21_plane exactly as defined above (scattered vs total field, projection on ê, area average,
   phase de-embedding). Add an **import** for the HFSS plane-integral / RCS CSV and overlay it on option 2 and on
   the periodic result, so the user sees the edge-effect difference directly.
4. Not in this phase (record as Phase 17 candidate): a true finite-array MoM (non-periodic layered Green's function
   with Sommerfeld integration, finite metal on an infinite substrate). For 11×11 loops at the exact 0.1 mm grid this
   is ~2×10^5 unknowns — needs ACA/FFT acceleration; a dense solve is not feasible. Say so in the UI tooltip.

Tests for Part F: (a) a large panel (e.g. 31×31) with the 80 % plane gives S11_plane and S21_plane within 0.3 dB / 3° of the periodic
S11/S21 away from the null, and the gap shrinks monotonically 11×11 → 21×21 → 31×31;
(b) the PEC-plate limit (pattern = full metal) reproduces the analytic flat-plate RCS 4πA²/λ² within 0.2 dB;
(c) the HFSS export of the user's 11×11 square-loop panel opens with the right sizes (110×110×1 mm, 121 loops).

## Tests
- Square loop above reproduces 9.64 ± 0.03 GHz and the 7.62–11.68 GHz band (slow test; fast test at coarse mesh
  checks shape/null ordering).
- STEP import of the user's file gives period 10 mm, loop 8.8/6.6 mm, area 33.88 mm².
- Exact-grid helper returns 0.1 mm for that loop.
- Energy check on a lossless case = 1 within 1e-6.
- Grounded patch case reproduces the Phase 12 G1 reference.

## Also in this phase (from the 2026-10-06 review)
- Ports: "port side" means ACTIVE aperture; the guard gap goes outside it. Update preset + UI label + the 1×3
  120 mm test (expect pass with 120 mm active).
- Optional "pointing correction" (off by default): re-aim each beam by the predicted miss, 1–2 iterations,
  report before/after. Cause documented: small-aperture squint toward broadside (wide array factor × falling
  element pattern), error ∝ 1/L².
- `examples/make_paper_board_26ghz.py` `_board_feeds_mm` writes feeds in the board-corner frame while the STEP is
  centred → monopole markers offset. Use the ports' own feed coordinates; regenerate the 26 GHz STEP; add a test
  that every feed marker lies inside its port rectangle in the STEP frame.
