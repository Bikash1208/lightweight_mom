# S4 both-beams baseline track — PHASE 1 (stored data only, no HFSS solve)

Branch `compute/sub6-baseline-s4` (from 4c7f1c4). Reference data: `../lwmom_s4_ref` (read-only). Scripts here: `common.py`, `p1_ab_power_ring.py`, `p1_d_attribution.py`, `p1_c_pointing.py`, `p1_c3_joint_fit.py`, `p1_e_levers.py`, `model_core.py`.
Labels: **STORED_HFSS_2026R1** = re-scored stored far fields (no AEDT started). **MODEL_ONLY** = af_jones-physics model. SLL = `sidequest_sll_def.analyze` (the other agent's first-null definition); my BOTH numbers reproduce -8.43 / -9.52 dB exactly.

## Machine
AEDT/HFSS **2026 R1 (v261, build 26.1.0)** installed — same version as the stored reference. RAM 127.8 GB total, 110.3 GB free at start. 24 cores / 32 threads. Licence: **NOT_CHECKED** (no AEDT session started in Phase 1; no ANSYSLMD_LICENSE_FILE variable set). PyAEDT 1.9.0 in `.venv` (Python 3.10). Git and Python were not installed and were added with winget. Full-board peak 69 GB + 15 GB = 84 GB < 110 GB free, so Phase 2 memory is OK; solves will be slower than the other server (it requested 46 cores).

## a. Power budget (STORED_HFSS_2026R1; Prad = integral of the FFD, normalised to accepted power via the HFSS peak gain; checked: gain-scale + 10log10(Pacc) is constant over all 7 states, so FFD is incident-normalised and the scaling is consistent)
| state | Pacc % | radiated % of Pacc | RHCP / LHCP % | below horizon % | remainder % of Pacc |
|---|---:|---:|---|---:|---:|
| BOTH | 97.89 | 35.77 | 19.21 / 16.56 | 4.41 | 64.23 |
| A_ONLY | 98.04 | 40.28 | 21.30 / 18.98 | 4.81 | 59.72 |
| B_ONLY | 96.63 | 37.71 | 19.50 / 18.20 | 4.67 | 62.29 |
| A_ONLY_IDEAL | 98.10 | 31.06 | 18.80 / 12.26 | 5.07 | 68.94 |
| B_ALL_OFF_IDEAL | 98.26 | 30.92 | 18.50 / 12.42 | 4.98 | 69.08 |
| ALL_PLAIN | 98.36 | 25.59 | 12.76 / 12.82 | 5.59 | 74.41 |
| B_ALL_SHORT | 98.31 | 30.65 | 18.36 / 12.29 | 5.02 | 69.35 |

Remainder split: conductors are PEC in the model (0 loss); dielectric loss is bounded above by 2.4 % of accepted power (RT5880 tan d 0.0009, alpha <= 0.0795 Np/m over 150 mm, MODEL_ONLY bound); so >= ~62 % of accepted power in BOTH is absorbed in the absorber ring (r 156-186 mm, 30 mm tall, matched eps=mu=1.5 loss tangent 0.667) — mostly unconverted surface wave. Near-field r*Ez^2 profile (BOTH, 2 mm above copper) falls only about 2 dB from r~11 mm to r~145 mm, consistent with most of the wave reaching the edge.
**Where the radiated 35.8 % goes (BOTH):** the two beams (within 12 deg of each peak) carry only 3.45 % (RHCP) + 3.91 % (LHCP) = 7.4 % of accepted power. The structure/ring term (ALL_PLAIN) alone radiates 25.6 %, of which 5.6 % is below the horizon. Gain (9 dBi) is limited by this conversion, not by hologram quality alone. (User set no gain target; reported only.)

## b. Source of the theta 44-48 deg ring (STORED_HFSS_2026R1) — NOT SETTLED by stored data
Facts: ALL_PLAIN far field is axisymmetric (phi-std 0.2-0.3 dB at 44-60 deg), >99 % m=0, E_theta only (|Ephi/Etheta| 0.05-0.07), so RHCP = LHCP (-3.8 dBi each, phi-averaged total -1.7 dBi at 48 deg). It is broad (44-100 deg within about 3 dB) and 24 % of its radiated power is below the horizon. Its m=0 amplitude at 40-52 deg changes only -0.5..+1.5 dB between all 7 states (BOTH 0.0 dB), i.e. it does not follow the cell state. That does not separate the candidates, because even the edge-reaching surface-wave power changes by only ~0.4 dB in amplitude between states (beams take only 7 % of power). Order-of-magnitude check: the pin protrudes 8.9 mm above the top surface (12.1 mm above ground), a monopole of kh~1.4 rad; a 25 % radiated fraction (about 12 ohm of 50 ohm) is in the plausible range for it. A below-horizon share of 24 % and a plateau out to 100 deg also fits edge/absorber diffraction. Verdict: feed-direct is somewhat favoured but UNPROVEN.
**ONE diagnostic solve that would settle it (not run):** ALL_PLAIN geometry with the feed pin cut flush to the top surface (z=0, port unchanged), same absorber. Ring falls >= 10 dB -> feed-direct; ring unchanged within ~1 dB -> edge/absorber (then second step: move/extend absorber). Needs ~63 GB and ~20 min on the 46-core server; expect longer on 24 cores.
Why it matters: the ring alone is -12.8 dB relative to the 9 dBi beams, so -15 dB is impossible without lowering it >= 2.2 dB even if everything else were perfect (ring also adds coherently to cell lobes).

## c. Pointing: effective beta per hand (MODEL_ONLY fit to STORED_HFSS_2026R1 beam angles)
HFSS BOTH beams with a sub-grid estimator (same estimator applied to the model, so grid bias cancels): RHCP 31.58 deg (phi 359.2), LHCP 35.94 deg (phi 180.4); model nominal (beta_t = beta_d = 1.10 k0) gives 30.24 / 30.14.
Joint fit of the physical/effective beta_t per hand (A cells for RHCP, B cells for LHCP; layout kept at 1.10 k0): **beta_eff,A = 1.158 k0 (1.14-1.18 for +-0.4 deg), beta_eff,B = 1.286 k0 (1.28-1.29)**. Verified with the af_jones-physics model (same code path as `af_jones.py`, extended per hand): reproduces 31.58 / 35.94 exactly. 
**Layout beta_d per hand that puts both beams at 30 deg.** Two answers, both MODEL_ONLY (beta_t fixed at the fit):
- Rule beta_d = beta_eff (r-term cancels): **A 1.158, B 1.286 k0** -> model peaks 29.90 / 28.26 deg, SLL -15.5 / -16.0 dB. Not exact for B: the A cells' cross-hand (LHCP) background changes with beta_t,A and shifts the LHCP peak ~1.9 deg.
- Angle-exact least-squares: **A 1.182, B 1.225 k0** -> 30.24 / 30.14 deg (nominal) but SLL -15.2 / -13.9 dB; this pair chases the peak estimator and the other hand's thinned cells, and the +-0.4 deg intervals above do not include that jitter.
Phase 2 should try beta_d = beta_eff first (A 1.16, B 1.29), then adjust B by HFSS result (about 30 deg of LHCP peak per 1.0 k0 of beta_dB near here, so B 1.225-1.29 brackets 28-30 deg). B is fragile: the model beam moves fast above 1.28 k0. Do not fix B for a final layout until the diode-cell result is in (as instructed); A can be tried first.
Layout check: the model `beam` and `rot_deg` (mod 180) match `v31/layout_R150_th30_v31.csv` (516 rows; 0 beam mismatches, max rotation difference 0.0 deg), so the model is the layout HFSS solved. `af_jones.py` itself was not run unmodified: it has one beta and a coarse grid; its physics was reimplemented in `model_core.py` with per-hand beta and the 2-deg HFSS grid.
FLAG: fitted beta_eff 1.16 / 1.29 is far from the supercell beta/k0 ~1.10 (1.08-1.12); the angle error is therefore NOT a plain beta error — it is an effective per-hand correction that may be absorbing mutual coupling / amplitude effects. The model cannot explain why B is further off than A with one physical beta. The HFSS iteration will show whether the retarget transfers.

## d. Lobe attribution (STORED_HFSS_2026R1 complex fields, linear superposition: S = ALL_PLAIN structure, A = B_ALL_OFF_IDEAL - ALL_PLAIN, B = BOTH - B_ALL_OFF_IDEAL; interaction folded into A, B; ASSUMPTION)
Levels in dBi per component at the lobe direction:
| hand | (th,ph) | BOTH | S | A | B | S+A | S+B | reading |
|---|---|---:|---:|---:|---:|---:|---:|---|
| RHCP | (10,164) | +0.61 | -15.0 | -5.2 | -5.7 | -4.4 | -5.6 | A and B add coherently (+6 dB); structure negligible |
| RHCP | (40,330) | +0.23 | -8.6 | -6.7 | -8.0 | -1.6 | -4.4 | all three add |
| RHCP | (18,260) | +0.19 | -13.5 | -3.9 | -11.6 | -1.4 | -7.3 | A, B small help |
| RHCP | (46,156) | +0.05 | -5.3 | -11.6 | -11.1 | -2.0 | -2.2 | structure (ring) |
| RHCP | (34,290) | +0.01 | -15.6 | -6.8 | -8.2 | -4.3 | -5.2 | A+B |
| LHCP | (18,130) | -0.53 | -12.4 | -8.9 | -7.1 | -6.0 | -3.9 | B dominant, A adds |
| LHCP | (36,146) | -0.63 | -12.0 | -7.6 | -11.3 | -3.6 | -5.8 | A dominant |
| LHCP | (32,228) | -0.87 | -15.1 | -11.8 | -5.6 | -7.9 | -3.8 | B |
| LHCP | (48,80) | -1.17 | -3.9 | -25.0 | -9.8 | -4.7 | -0.6 | structure (ring) |
| LHCP | (44,270) | -1.53 | -5.3 | -18.7 | -13.7 | -4.0 | -2.5 | structure (ring) |
Both worst lobes (RHCP +0.61 at (10,164), LHCP -0.53 at (18,130)) are near-broadside **cell** lobes, not structure. Even with the B cells plain (B_ALL_OFF) RHCP SLL is -9.15 dB at (44,188): the ring/cell floor. B_ALL_SHORT and the 24 residual B cells (A_ONLY_IDEAL) contribute < -13 dBi at these lobes. SLL of component sums: LHCP S+B -8.5 dB (ring lobe (48,78)); B alone -11.2 dB; RHCP A alone -12.2 dB.
Scaling S in the stored fields (first-null SLL): S -3 dB -> R -8.5 / L -9.9; S -10 dB -> -8.6 / -10.4; S -20 dB -> L -10.5. At other scalings (-6 dB, -20 dB, removed) the first-null mask fails for RHCP (reads about -3 dB, see FLAGS). A cone metric (15-20 deg exclusion) moves SLL by about +-1.5 dB with mixed sign, and with S removed LHCP gets 1.4 dB worse. Honest range: **ring scaling alone moves SLL by about +-1.5 dB, not a large gain** at the current layout; ranking unchanged.

## e. Phase 2 fix ranking (expected SLL gain; MODEL_ONLY unless stated)
Model reference: nominal pointing -15.5 / -17.2 dB; as built (beta_eff, pointing error) -12.7 / -9.0 dB (LHCP -9.0 vs HFSS -9.5: good match; RHCP model -12.7 vs HFSS -8.4: the rest is ring + coupling).
1. **Per-hand beta retarget** (beta_d = beta_eff: A 1.16, B 1.29; or 1.18/1.23): model SLL -12.7/-9.0 -> -15.5/-16.0 (+2.8 / +7.0 dB) or -15.2/-13.9; also fixes the one-theta gate (to within the 2-deg grid only if B is tuned in HFSS). Highest value, layout-only. Iterate A first.
2. **Restore the 16 gate=1 discs to anchors (gate=0)** (own fix; changes the reference): model, only the 16 CSV gate=1 cells zeroed: nominal -13.8/-13.9 vs -15.5/-17.2 with them present (+1.7 / +3.3 dB); with beta retarget -14.1/-13.2 vs -15.5/-16.0 (+1.5 / +2.8 dB).
3. **Seed / dither width**: seed spread at W=30 mm is -15.0..-17.9 dB (R) and -15.0..-17.2 (L); W 20-30 mm best; ~1-2 dB, free but model-only and HFSS may differ.
4. **Feed treatment** (pin protrusion; only if the b-diagnostic says feed): <= ~1 dB now, but REQUIRED later: ring is -12.8 dB below the beams, so it caps SLL above -15 dB until lowered >= 3 dB (more if coherent with cell lobes).
5. **Edge/absorber treatment**: same effect as 4 if the ring is edge-driven; no evidence either way.
Expected honest outcome: 1+2+3 together may reach about -12..-14 dB in HFSS; -15 dB in both hands likely also needs the ring reduced (4 or 5). Not promised.
Phase 2 start requires the user's go.

## Flags
**Main-lobe shoulder:** both BOTH beams have a broadside-side shoulder about -2.7 dB below peak, about 16 deg from the peak ((16,0) RHCP, (20,184) LHCP), probably from the beta-mismatch defocus. The first-null mask can count it as main lobe. Phase 2 gating risk: a fix that opens a null there can make SLL read about -3 dB; check the pattern, not only the number.
**Licence:** server file `ansyslmd.ini` = `1055@10.100.105.13`; TCP port 1055 reachable; `lmutil lmstat` (read-only) answered HOST_NOT_FOUND for the server name `ansys-license-server` (name does not resolve on this PC). Checkout is NOT verified. Phase 2: one attempt; stop and report if it fails.
All attribution and ring-scaling numbers assume linear superposition of separately solved states. The model has no feed ring, absorber or mutual coupling (optimistic). No HFSS solve was run. `results/sub6_cell_param/idc_cells/HANDOFF_2026-11-19.md` not edited.

## Files
Current: `common.py`, `model_core.py`, `p1_ab_power_ring.py`, `p1_d_attribution.py`, `p1_c_pointing.py`, `p1_c3_joint_fit.py` (+ `p1_c3_fit.json`), `p1_e_levers.py`, `p1_*.json`. **Superseded (do not copy beta from them): `p1_c2_fit_levers.py`, `p1_c_fit.json` (1.141/1.253, single-hand fits).** A repo-local git identity was set (user.email from the account, name "S4 baseline track"); nothing global. The read-only reference worktree `../lwmom_s4_ref` is clean (no changes).
