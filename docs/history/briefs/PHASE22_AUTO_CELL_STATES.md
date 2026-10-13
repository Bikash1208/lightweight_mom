# Phase 22: automatic coded-cell states (the user enters only the cell size)

## Why
Coded (1-bit/N-bit) design in the app needs the reflection coefficient of each cell state at the patch plane. Today there are three ways to supply it, and none is automatic for the common case:
- **"Enter two measured values".** The user must run HFSS or the Unit cell tool, add the 2·k0·5 mm de-embedding by hand, and convert to real/imag.
- **"Best pair from the square-patch library".** This needs a pre-built library on a coarse grid. Its screening ("clear of patch-plane resonance") also tends to refuse exactly the near-resonant square that gives a 180° pair. For the reference board that square is 5.85 mm at 10.54 GHz on RT5880 1.575 mm with a 6.4 mm cell: −42.0° against the empty cell's +138.1°, a separation of 179.9°.
- **"Square + annular slot".** This also needs libraries.

Goal: the user picks the substrate, cell size, frequency and number of bits. The app finds the cell geometries and their Γ itself, with mesh-converged accuracy, and shows them. The user can still override with HFSS values.

## Reference case (acceptance test)
- **Inputs:** RT5880 (εr 2.2, tanδ 0.0009, 1.575 mm), PEC ground, cell 6.4 mm, 10.54 GHz, 1-bit, square-only family.
- **Expected:**
  - state 0 = empty, Γ = 1∠+138.14° (exact grounded slab, patch plane);
  - state 1 = square L ≈ 5.85 mm, Γ ≈ 0.996∠−42° ± 5°;
  - separation 180° ± 5°.
- **Hand reference:** `scripts/adhoc/square4/` and `scripts/adhoc/fss`-style runs. The n=64 (0.1 mm) MoM values were:

  | L (mm) | Phase at patch plane (n=64) |
  |---|---|
  | 5.0 | 94.8° |
  | 5.4 | 58.8° |
  | 5.6 | 25.1° |
  | 5.8 | −23.4° |
  | 6.0 | −76.9° |

  The 0.2 → 0.1 mm Richardson shift is −4 to −6°, which puts the converged crossing near L ≈ 5.85 mm.
- **End to end:** feeding the auto pair into the 2×2 grid (300 mm side, 0 mm gap, beams at 30°, φ 45/135/225/315, RHCP, model β) must reproduce the Phase 21 acceptance numbers. That means 2209 cells per port, θ within 1°, and about 27.7 dBi.

## Work items

### A. Backend: `src/hologram/auto_states.py`
1. **`empty_state(substrate, f, period)`**: the exact grounded-slab Γ at the patch plane, from the existing `scattering_from_fields` with zero current.
2. **`find_square_for_phase(target_phase_deg, substrate, f, period, L_grid_mm=0.05, ...)`**:
   - **Coarse scan:** L from 0.3·P to 0.97·P in 0.2 mm steps at a 0.2 mm mesh. Choose cells per side so that each L lands exactly on the grid, or pick the nearest exact mesh. Phase is monotonic through resonance, so unwrap the phases.
   - **Bracket the target and refine:** bisection or secant at the 0.1 mm mesh, then a Richardson correction from the 0.2/0.1 pair, as in Phase 12.
   - **Snap L to the fabrication grid** (`L_grid_mm`, default 0.05 mm). Report the snapped L's Γ, not the interpolated one.
   - **Report:** the phase sensitivity (°/0.01 mm), |Γ|, the mesh-uncertainty estimate (Richardson shift), and a `NEAR_RESONANCE` label when the sensitivity is above 1°/0.01 mm. This is a warning, not a refusal.
   - **Keep it fast.** Cache by (εr, tanδ, h, P, f, mesh, L), on disk under the app's storage, so re-runs are instant. Target under 2 min for a fresh 1-bit pair on a laptop. Run solves in the existing job runner with progress messages.
3. **`auto_states(n_bits, ...)`**:
   - **1-bit:** state 0 is empty and state 1 is the square closest to empty + 180°.
   - **2-bit:** 4 states at empty + 0/90/180/270°. **3-bit:** 8 states.
   - If a target is unreachable with squares, return the closest one with an `UNREACHABLE_LEVEL` label and its phase error. Do not fail silently.
   - Return the same `TwoStateChoice` (or N-state) objects the coded synthesis already consumes, with `geometry = {"kind": "rectangle", "L_mm": ..., "W_mm": ...}` so export draws it.
4. **Screening.** With beta_mode "model", the near-resonance screen in `screen_candidates` must not refuse the auto-found square. Keep the refusal for beta_mode "state", where it matters, and document the reason.

### B. UI (Design method → Coded)
1. **A new first option under Two-state cells: "Automatic (squares; computed from the cell size)",** default **on**. It uses the Board section's substrate, cell size and frequency; nothing else is needed.
2. **After the solve, show a table:** state, geometry (empty / L×L mm), |Γ| dB, phase at the patch plane, separation, sensitivity, and labels.
   - Add an **"Use these as uploaded values"** button. It copies the numbers into the "Enter two measured values" fields so the user can edit or override them with HFSS values.
3. **Advanced: an optional "fix L" field.** The user types L and the app computes that cell's Γ instead of searching, for when the user has a fabricated cell size in mind.
4. **Keep "Enter two measured values"** for HFSS overrides. Change its help text to: "phase at the patch plane (de-embed the Floquet port to the top face of the substrate)".
5. **2-bit/3-bit with Automatic:** use the same auto search for 4 or 8 states, instead of the library requirement.
6. **Unit cell tool (`unit_cell.html`):** add a "phase at the patch plane" column, the reported phase plus 2·k0·d_ref, so its output can be pasted directly.

### C. Two small fixes found while driving the app
1. **Step 1 lock.** Multi-port does not use the Step 1 .fld file, but Step 2 stays locked until a file is analysed. In multi-port mode, make Step 2 reachable without Step 1, and enable the Frequency field there. Single-port behaviour is unchanged.
2. **Multi-port results, two directivity numbers per port.** Each port shows "Directivity 27.8 dBi, θ 29.8°" in the beam table and "Directivity 26.7 dBi, peak at θ = 29.3°, φ = 132.2°" in the text under the pattern. These come from two pattern evaluations on different angular grids. Compute both from the same fine-grid pattern (0.1° pointing grid), or label the coarse one clearly as the coarse plot grid.

### D. Help text for "Surface-wave beta used for the phase"
Add a one-line explanation under each button:
- **Model (default):** β of the actual mixed surface of both states, solved by the app. This is the most accurate.
- **Paper:** Mishra et al.'s eq. 1 closed form, from εr and thickness only. Use it to compare with the paper.
- **Each state's own:** β is path-integrated cell by cell from each state's own β. Use it only when both states are far from resonance.

## Constraints
- Branch from `phase21-multiport-layouts-coded`. Commit per part with the fast suite green, no push, no tag moves.
- The existing options and results stay byte-identical, apart from C2's display change.
- Model-only labels as usual; nothing claims HFSS agreement.

## Report back
- The auto-pair table for the reference case versus the expected values above, with timings (fresh and cached).
- The 2-bit auto states for the same board.
- Screenshots of the new Coded panel and of the multi-port run reachable without Step 1.
