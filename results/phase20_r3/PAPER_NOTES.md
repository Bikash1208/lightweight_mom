# R3 paper notes — full texts read 2026-10-19 (IEEE Xplore, IIT Indore access)

Read in full by the advising session through the user's institutional access; these are design-relevant summaries in
our own words (no figures or verbatim text copied). PDFs are not in the repo (copyright; cite by DOI). Treat numbers
below as **published simulation/measurement values**, not our results.

## CP-1 — Crossed-bowtie 2-bit CP RRA, IEEE TAP 2023, DOI 10.1109/TAP.2023.3291773 (IEEE doc 10179269)
- Mechanism: element-rotation (PB) — CP phase = 2 × rotation; four rotation states at 45° steps from diode switching.
- Cell: four pairs of bowtie patches arranged annularly (two interlaced crossed bowties), **8 diodes** (MACOM
  MADP-000907-14020), but only **two bias voltages per cell**: diodes in each group are oriented so that one pair is
  forward- and the other reverse-biased for a given polarity (±V on two lines → four states). Outer diode ends grounded
  through **blind vias placed at current minima** (shown RF-negligible); inner ends fed through two plated through-holes
  to a bias layer under the ground with microstrip low-pass filters (isolation to ~40 dB).
- Stack: **multilayer** (five layers; a 1.52 mm core plus 0.2 mm layers), period 15 mm. Not a single-substrate planar
  model; the RF-active part is the bowtie layer over ground.
- Anisotropy: activating one bowtie pair gives ≈ 180° x/y reflection-phase difference over a wide band (the PB
  condition). Bandwidth comes from a very small gap between interlaced bowties (≈ 0.2 mm); 1-bit single bowtie is narrow.
- Element (sim., 30° incidence): cross-pol CP ≤ −13 dB from 8.5–10.3 GHz, four parallel phase curves 90° apart,
  insertion loss mostly ≤ 1.5 dB; cross-pol ≥ 10 dB down up to ~40° incidence.
- Array (16×16, 240 mm, offset 30° RHCP horn, F/D 1.02): measured 24.5 dBic at 9.9 GHz, 35 % peak aperture efficiency,
  3 dB AR bandwidth 22 % (8.6–10.7 GHz), SLL ≈ −20 to −22 dB, scanning ±60° with ≈ 3.2 dB scan loss; 2-bit
  quantisation costs ≈ 1 dB vs ideal.
- Sensitivity the authors found: the inter-patch gap (0.21 mm nominal) and the via location dominate the measured shift.
- Fit for us: **partial** — topology modelable on one substrate if the thin layers are folded into one equivalent
  substrate, vias as labelled approximation; the 0.2 mm gap will be mesh-fragile (our R26 lesson: < 0.35 mm). The
  **two-line bias trick (anti-parallel diode pairs)** is worth reusing for any 4-state cell.

## CP-2 — Multibit PB CP RRA element, IEEE AWPL 2025, DOI 10.1109/LAWP.2025.3606542 (IEEE doc 11152537)
- Mechanism: PB; one of **four diodes ON**, three OFF, emulates rotation by 0/45/90/135° → 0/90/180/270° CP phase.
- Cell: **regular octagonal patch** (r ≈ 3.7 mm) with a **centre blind via to ground** (DC return; at the patch centre),
  **four phase-delay lines at 45° spacing** around the patch, each with one diode (MADP-000907-14020, modelled as a
  series circuit); period 12 mm at 11.6 GHz. Bias: four lines, chip inductors (2 × 3.2 nH) as RF chokes instead of
  radial stubs. Upper substrate Arlon AD255 (εr 2.55) **1.5 mm** — close to our RT5880 1.575 mm; bias board below ground.
- Anisotropy check: with one diode ON, x′/y′ reflection-phase difference 180° ± 10° over 11.15–12.2 GHz, both
  magnitudes ≥ −1 dB.
- Element (measured in a waveguide simulator, 2×2 cells): co-pol CP insertion loss 0.42–0.65 dB (sim. 0.33–0.42 dB),
  < 1 dB over 11.1–12.5 GHz (≈ 12 %); cross-pol ≤ −15 dB; state-to-state steps 83–99° (S1–S2), 185–210° (S1–S3),
  275–295° (S1–S4). Oblique stability shown to 25° (cross-pol ≤ −15 dB).
- 3-bit extension: eight delay lines + eight diodes on a circular patch, 45° steps.
- Array (16×16 sim.): 24.1 dBic, 39 % aperture efficiency, AR < 3 dB over 22 %, ±60° scanning.
- Fit for us: **good** — planar single patch layer over one substrate of almost our permittivity/thickness; centre via
  at the voltage null = our `VIA_OMITTED_APPROXIMATION`; delay lines and diodes are planar loads. Scaling to 10 GHz on
  RT5880: expect period ≈ 13–14 mm (inside our 6–15 mm window), patch and lines larger; redesign with SMP1345 values,
  not by scaling dimensions.

## Same element, LP mode — Quad-polarisation RRA, IEEE AWPL 2025/26, DOI 10.1109/LAWP.2025.3646699 (IEEE doc 11309746)
- Same authors' octagonal 4-diode PB element (above). New finding: under **linear** incidence (feed rotated 22.5°),
  specific diode codes give **1-bit LP→co-LP** and **1-bit LP→cross-LP** (≥ −1 dB, ≈ 180° steps over ≈ 11.2–11.8 GHz
  and 10.9–12.25 GHz). Under LP incidence the PB phase acts oppositely on the two CP halves, so one CP handedness can be
  selected from an LP feed at an inherent 3 dB cost.
- Array (16×16, 192 mm, one fixed LP horn, FPGA control): LP 21.0 dBi / cross-LP 20.8 dBi (1-bit, SLL ≈ −11 dB,
  cross-pol only ≈ −9 dB); LHCP 21.5 / RHCP 21.4 dBic (≈ 20 % aperture efficiency, AR < 3 dB over 14 %); ±60° scan.
- Relevance to the user's "LP and CP as a user choice": **one cell can serve both** — 2-bit CP, but only **1-bit LP**.
  A 2-bit LP mode still needs its own cell (LP-1/LP-2) unless 1-bit LP is accepted.

## Consequences for R3 (advising-session recommendation)
1. CP-2 (octagonal patch + four delay-line diodes + centre via) is the lead CP candidate: closest to our board and
   solver, measured, 4 diodes. CP-3 (four-slot circle) stays as our own alternative; CP-1 is a reference/benchmark
   (multilayer, 0.2 mm gap).
2. Decide with the user whether LP may be 1-bit from the shared PB cell, or must be 2-bit (separate LP cell).
3. None of the papers reports diode-tolerance (±20 % C_T, L_s) analysis; our tolerance objective remains new work.
4. All three use MADP-000907-14020 (flip-chip, lower parasitics) rather than SMP1345 (SC-79, L_s 0.7 nH); expect the
   SC-79 package inductance to be harder at 10 GHz. Consider registering MADP-000907 as a second diode model (values
   from its datasheet, cited) so the cell can be compared with the published element directly.
