# Fresh server session — Phase 30: new unit cell for the multi-port HMSA

You are working in the `lightweight_mom` repo on the Windows server (our own MoM solver + the "HMSA Designer" app;
HFSS/AEDT is used only to verify). Start by reading `CLAUDE.md`, `docs/STATUS.md`, `docs/ROADMAP.md` and the latest
Phase 29 notes, then `git fetch` and branch `phase30-new-cell` from the latest `phase29-seven-hex` (last commit 7bf9529
or later). Session scope: **solver + design only** (no fabrication, bias, FPGA or hardware work).

## Repo rules (unchanged)
- One commit per item; fast test suite green before each commit; physics tests for anything new.
- Solver changes are **additive** and bit-identical for every existing case (use `tests/_platform_tol.py` tolerances).
- Honest labels: MODEL_ONLY until HFSS-checked; never call a fitted surrogate "physics".
- Heavy compute through `jobs/` with manifests; no ML in the solver core.
- Literature: **journals only** (no conference papers, no arXiv), cite by DOI, verify every DOI exists before citing;
  no PDFs committed. The user has IP-based journal access on this machine.
- Never move/delete tags; delete nothing without approval. Utf-8 + LF for all new text files.
- **Token discipline:** the last long run used 25 % of the weekly limit. Be terse, no exploratory sweeps beyond what is
  asked, coarse grids first, short reports.

## Why (HFSS evidence, user-run, 10.54 GHz design frequency)
Current cell: Phase 26 tensor cell, 6.4 mm lattice, RT5880 (εr 2.2, tan δ 0.0009) 1.575 mm, coax-probe monopole feeds.
- Two-feed board: P1 17.5 dBi / SLL −8.9 dB, P2 16.6 dBi / −8.1 dB (model was ~19.4 dBi / ~−12 dB).
- 7-port shared-aperture board: per-port 11.3–13.4 dBi, SLL −1.2 to −6 dB; the **raw** forward model matched HFSS
  gain (P1 13.4 vs 13.6 dBi). **Drop the −2.5 dB / +3.5 dB calibration margin** in
  `scripts/research/phase29_calibration.py` from all reporting (keep the file, mark it superseded).
  Blended steering works to about θ ≈ 12° at 10–13 dBi.
- compact7 board: 7.0–10.5 dBi per port; blended steering to θ ≈ 30° at 5–9 dBi.
- **Common failure on every board:** 48–68 % of each port's power is radiated near the board surface (θ 60–90°),
  ~98 % Eθ, with 6-/12-fold azimuthal ripple → the TM surface wave reaches the hexagon edge un-leaked. Only ~7 % of
  power lands in the main-beam cone. Opposite-hand power in the beam ≈ 25 %, AR at beam peaks 0.5–5.3 dB.
  Radiation boundary was ≥ λ/2 from the board (not an artefact).
- Root cause: the current cell's impedance swing is too small (anisotropy only ~1.2–1.5 : 1, M range ~0.16–max
  insufficient), so α (leakage per length) is too low for a 150 mm (~5.3 λ) board.

## Goal
Choose a new unit cell (geometry + substrate) that, on the 150 mm hexagon at 10.54 GHz, can:
1. leak the surface wave fast enough that **< 10 % of launched power reaches the edge** (needs roughly α/k ≥ 0.03,
   M up to ~0.5 with a taper);
2. realise a **full tensor** (Xρρ, Xρφ, Xφφ) by size + rotation with anisotropy ≥ ~3 : 1, for CP with AR < 3 dB;
3. stay well-behaved for surface waves arriving from different directions (multi-feed), ideally **directionally
   selective** (Frey-style) to reduce cross-hologram sidelobes;
4. give a smooth curve from M ≈ 0 to M_max;
5. be single-layer, printable, pitch ≤ ~λg/5.

## Candidates to screen (journal sources)
Geometries:
- (a) elliptical patch ("grain of rice"), size + axial ratio + rotation;
- (b) elliptical patch with centre slot — Frey et al., IEEE OJAP 2023, doi:10.1109/OJAP.2023.3261903;
- (c) patch with cross slot (and patch with single slot).
Substrates (≈ thin, high-εr, as in the efficient published designs):
- RT5880 1.575 mm (baseline), RO3006 0.635 mm (εr 6.15), RO3010 0.635 mm (εr 10.2).
Background: Faenzi, Minatti, Maci et al., Sci. Rep. 2019, doi:10.1038/s41598-019-46522-z (cell families, high-εr
substrates, M-taper, ~70 % efficiency); Minatti et al., IEEE TAP 2015 (space MTS antennas); Fong et al., IEEE TAP 2010
(scalar/tensor holographic impedance surfaces); Teniou, Sauleau, Ettorre et al., IEEE TAP 2017 (tensorial metasurfaces
for aperture-field implementation). Verify and record the TAP DOIs yourself. Optional (only with a journal source):
via-loaded (mushroom) cell for extra reactance range.

## Work items (one commit each; ■ = stop and report)
**0. Plan + docs.** ROADMAP/STATUS: open Phase 30 (design change: new cell), record the HFSS evidence above, mark the
calibration margin superseded. Docs only.

**1. Literature table (journals only, DOIs verified).** For the candidate cells: frequency, substrate, pitch, achieved
reactance range / M, anisotropy, polarisation, aperture size, measured efficiency, leakage/taper method. 5–10 lines of
conclusions. Also: how published boards suppress edge leakage on small apertures (taper, absorber rim).

**2. Cell screening with our MoM (MODEL_ONLY).** If the eigenmode/dispersion solver needs a new geometry type
(elliptical patch with slot, cross slot), add it as an additive cell family with physics tests (its own commit).
For each geometry × substrate, coarse grid first, at 10.54 GHz:
- tensor reactance (Xρρ, Xφφ, Xρφ) vs size and rotation; average X̄ and swing;
- maximum achievable M and the M(size) curve near zero;
- α/k vs M (leakage rate), and β/k;
- anisotropy ratio;
- sensitivity to surface-wave arrival angle 0–90°;
- required pitch;
- edge-power fraction on a 150 mm hexagon with an ideal M taper at the best M_max;
- dielectric/conductor loss estimate.
Heavy runs as a server job via `jobs/`.
**■ Stop and report:** a ranked table, the top 2 candidates with trade-offs, and a recommendation. No board design and
no full library before the user approves.

**3 (after approval).** HFSS eigenmode check package for the chosen cell (predictions first: β, and α where possible,
for 3–4 representative cells), plus a short README for the user to run it manually.

**4 (after HFSS agrees).** Full cell library + angle-dependent LUT on the server, held-out validation as in R4.

**5 (after 4).** Redesign the 7-port shared-aperture board: P1 broadside, six outer beams on an even ring at θ ≈ 17°
(~18° wide), M taper for < 10 % edge power, CP target. Forward-model per-port gain/SLL/AR and blended coverage to
±25°, then a STEP package (STEP file + port table + steering LUT) like `hfss_exports/phase29_sevenport/`.

## Do not
- Design boards or build the full library before the item-2 stop is approved.
- Touch the parked RIS / Phase 28 track.
- Cite anything you have not verified; no conference/arXiv sources.
