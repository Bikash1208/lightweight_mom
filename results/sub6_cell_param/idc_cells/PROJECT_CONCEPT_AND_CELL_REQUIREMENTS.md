# SUB6 switchable dual-circular-polarisation holographic board: concept, required outcome, required unit-cell properties

Status date: 2026-11-19. Author: executor S2 (lab B). Labels: MEASURED = full-wave result, MODEL_ONLY = analytic or first-order chain, PREDICTED = derived requirement, not yet demonstrated.

## 1. Project concept

A flat, circular holographic metasurface antenna (HMSA), **150 mm radius**, at about **5.8 GHz** (the operating frequency is free, see section 3), on RT5880 3.175 mm over a ground plane with an RT5880 bias layer.

- **How it radiates.** A centre feed launches a bound surface wave across the board. The surface is a lattice of small copper cells (pitch 11.6 mm, about 0.22 wavelength, about 525 cells). Each cell is anisotropic: its surface impedance differs along x and y. Modulating that impedance leaks the surface wave into a beam. Rotating a cell changes the phase of the circularly polarised leakage, so the cell rotation is the hologram.
- **Dual circular polarisation, static.** Two sectors of the board are synthesised: sector A radiates RHCP at (30 deg, 0 deg) and sector B radiates LHCP at (30 deg, 180 deg). The default is a two-quadrant board. Four-quadrant boards are parked.
- **Switching.** Each cell gets diodes. With all cells ON, both beams radiate (two-beam state). With sector B's cells switched OFF, only beam A radiates (one-beam state). Beam A's gain must not drop in the one-beam state. The model says it rises slightly (14.13 to 14.30 dBi at ON leakage 1.33 Np/m, MODEL_ONLY).
- **Novelty.** The user originally asked for an interdigital or meander-type cell. As of 2026-11-19 interdigital is **not** required.

## 2. Required outcome (board level)

- Two beams of equal gain, one RHCP and one LHCP, at the same polar angle, azimuths 0 and 180 deg.
- Beam angle chosen from {25, 30, 35, 40} deg. Gain cost relative to 30 deg is small (+0.2, -0.24, -0.53 dB for 25, 35, 40 deg, projected aperture only).
- Sidelobe level as low as possible, target -15 dB or better. The soft (dithered) sector boundary of about 30 mm gives -18.5 dB at no conversion loss (MODEL_ONLY, at 1.33 Np/m).
- Two states: both beams ON, and beam A alone (sector B suppressed, at least about 10 dB suppression of B).
- Bandwidth is not a constraint. Operating frequency and beam angle are chosen together at the end.
- Fabrication limits: minimum feature 0.3 mm, every cell outline inside radius 5.65 mm (so any rotation fits the 11.6 mm lattice), RT5880 only for the cell and bias layers.

## 3. Required unit-cell properties (PREDICTED, from our own chain and measurements)

These are the numbers a candidate cell has to meet. The units are per metre of surface-wave travel, with R = 0.15 m.

**ON state (anisotropic):**
- Leakage constant **alpha at least 1.2 Np/m** at the common operating frequency (about 30 percent of the surface-wave power radiated over the board; 50 percent needs 2.31 Np/m). Equivalent to alpha/k0 of about 0.010 or more. Conversion = 1 - exp(-2 alpha R).
- The x/y wave numbers of the cell must differ. The anchor (MEASURED, unit-cell eigenmodes) has beta/k0 of about 1.134 along one axis and about 1.067 along the other: a **split of about 0.065 with a mean near 1.10**. A candidate needs a comparable split (at least about 0.05) and a mean at least about 1.08.
- The anisotropy must come from a slot resonance just **above** the operating frequency (the anchor's resonance sits about 0.15 GHz above its crossing). A resonance below the band leaves the wave barely loaded.
- The bound wave must be **single-moded** at the operating frequency, with its mode identity trackable across at least 3 phase points and continuous in frequency and Q. Crowded multi-mode spectra (B1p5, D3p6) make the leakage ambiguous.
- Leakage is extremely sensitive to beta (0.005 in beta_y changes alpha by about 28 percent), so the supercell (full-wave) value is the reference. The 5.8 GHz first-order chain is a screen only, valid only where delta/|B| is at most about 0.4.

**OFF state (switching):**
- Leakage **alpha_OFF at most alpha_ON / 12** at the same frequency (about 0.1 Np/m for alpha_ON of 1.2). That corresponds to a field ratio of about 0.3, about 10 dB suppression of beam B.
- The OFF cell must be **isotropic**: beta_x equal to beta_y within mesh noise, so rotating the cell changes nothing. A cell with four-fold or higher rotational symmetry has an isotropic surface-impedance tensor (homogenised). At about 81 degrees of phase per cell, spatial dispersion can leave a small residual, so a supercell OFF check is still required.
- A single slot made anisotropic or isotropic only by a diode cannot meet this: the ideal open twin of D3p6 still leaked at least 1.17 Np/m (about 1.3 times contrast, MEASURED in the supercell).
- Real diodes (BAR90-02ELS, OFF sheet about 22.23 - j335 ohm, crop testbed -10 V gives 21.6 - j297.5 ohm) give less contrast than an ideal open circuit.

**Realisation constraints:**
- Each diode needs two DC-separate nodes, each with a bias path. A diode across a closed slot in one connected disc cannot be biased.
- Four-fold OFF symmetry needs matched parasitics on all four arms (identical diodes, or matched caps of about the diode OFF capacitance).
- Shorted ground vias make a mushroom surface and change the dispersion. Bias vias must be RF-isolated (choked) and included in the cell model.
- The cell must stay inside radius 5.65 mm and 0.3 mm features, with room for diode lands.

## 4. What this means in one line

We need a cell that behaves like the anchor (a closed slot tuned just above the operating frequency, alpha at least 1.2 Np/m) in the ON state and like an isotropic cell (alpha at most 0.1 Np/m) in the OFF state. The anchor does the first, symmetry does the second, and no tested cell yet does both.

## 5. CORRECTION 2026-11-19 (user): diodes only in key areas, not per-cell switching

User statement: all board cells are the same kind of static cell, no half-static / half-switchable split; the diodes sit only in a few key areas (groups or individuals) whose only job is to turn the second beam (LHCP, beam B) on or off. This supersedes sections 1 and 3 where they say every cell carries diodes, and retracts "essentially every B cell must switch" and the per-cell alpha_OFF <= alpha_ON/12 requirement (those held only for per-cell in-place switching of all of sector B).
- All hologram cells are anchor-type (MEASURED alpha 1.11 / 1.45 / 1.83 Np/m). Sector A and B are only regions of the hologram.
- Two tracks to be studied (user: "lets try both"):
  1. **Gate band** near the feed on the B side: one or two rows of diode cells that pass the surface wave to the B region (two beams) or block it (one beam). Gate cell needs: pass state beta/k0 within about 0.02 of the anchor mean (about 1.10), single-moded, >= 90 percent transmission; block state no bound mode at f_op in any in-plane direction (check 0 and 45 degrees too), driven transmission <= 10 percent. Constraint that decides it: B cells between feed and gate keep radiating (anchor alpha 1.45: r < 50 mm carries about 40 percent of B's radiated POWER, about 4 dB in power; the 10 dB criterion is on the PEAK, and a small inner aperture gives a broad weak beam: crude scaling peak suppression about 30*log10(R/rho), about 14 dB at 50 mm, MODEL_ONLY, to be verified by the board study; the residual power becomes a broad lobe, check wide-angle level and SLL), so the gate must sit very close to the feed or the B side inside it must be non-leaking cells. Gate must span the B half plus the dither zone. First cell candidate: disc with a diode-switched via to ground (mushroom stopband when connected).
  2. **Selected cells** (groups/individuals inside the board): power-conserving model needed (switching a cell off does not remove its power; the bound wave carries it on). Study beam-B suppression vs switched fraction and placement. Still needs an in-place ON/OFF cell unless the selected cells act as a blocking barrier (then it is track 1 at another radius).
- Risks to check: blocked power reflects toward feed and sector A (feed match, beam A ripple/gain), edge scatter into sidelobes, wave slipping round the gate ends.
- Both studies are board-model questions first; the board model that gave 14.13/14.30 dBi and -18.5 dB lives on the other agent's branch (scripts/research), not locally.
- **Update (user, 2026-11-19): three states.** BOTH beams, RHCP only, LHCP only. In each one-beam state the other beam is at least 10 dB below its own two-beam peak and the surviving beam's gain must not drop. Gate track therefore needs two independently switched gates (one per half); selected-cell track is symmetric.
- **Update (user, 2026-11-19):** (1) The board is judged on the FULL 150 mm aperture only (no cropped/small-board tests; a small aperture cannot convert enough into surface wave). Reading: the gate must not sit inside the feed's launch zone where the wave has not yet become a bound surface wave. (2) The gate cell passes, blocks and preferably ABSORBS in the block state (reflected power hurts feed match and spills into sidelobes). (3) The earlier "transmission <= 10 percent" gate target is too loose: beam-B peak in the one-beam state is about (rho/R)^1.5 + t(1-(rho/R)^1.5) (t = field transmission, MODEL_ONLY); for -10 dB: rho 40 mm needs power transmission <= about 4 percent, 50 mm <= 2.3 percent, 60 mm <= 0.7 percent, rho >= about 70 mm unreachable even with a perfect gate. Usable window: launch-zone radius to about 65 mm; if the launch zone ends beyond about 65 mm the gate concept cannot reach 10 dB. Several rows needed (about 3 dB block per row at 5 rows, <= 0.1 dB pass per row, beta within 0.02 of the anchor). Pass state should be the anchor itself (rotated with the hologram) with the diode switching loss, not a plain disc.
- **Standing rule (user):** every engineering decision is taken by the advisor; ask the user only about requirements and pushes.

## 6. Requirement change 2026-11-20 (supersedes three-state wording in section 5)
Two states only: both beams on; LHCP off with RHCP unchanged (polarisation same). Diodes only on the LHCP (B) half. See handoff section 36.

**Update 2026-11-21:** diode count is now at most 150 (250 not allowed); see HANDOFF section 59 (plain-disc partition model: 150 diodes give about 5-7 dB suppression, 8.0 dB with perfect diodes; 9 dB needs about 233).
