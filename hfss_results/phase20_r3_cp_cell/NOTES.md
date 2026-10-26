# CP-2 cell — user's HFSS run (2026-10-23), first comparison by the advising session

Files: `hfss_bare_slab_theta0_30.csv` (no metal, no loads), `hfss_cp2_states1-4_theta0_30.csv` (states 1-4 one-hot,
MADP ON = R 5.2 ohm / OFF = C 0.025 pF via parallel-RLC sheets with negligible partner element, theta 0 and 30 deg,
phi 0, 8-12 GHz in 0.5 GHz steps). Columns: re/im of S(FloquetPort1:i, FloquetPort1:j), i received, j incident.
Comparison script `advisor_compare.py` and its output `advisor_compare_output.txt` (quick check, not the repo tool —
the session should rerun with `tensor_cp_compare`/R4 hook and confirm).

## Conventions (from the bare slab)
* Floquet mode 1 = TE, mode 2 = TM (matches `bare_grounded_slab_reference.csv` labels); no 180 deg TM flip.
* Reference-plane offset: HFSS phase = analytic − 2·k0·d·cos(theta) with **d = 1.236 mm** (fit over all f, both
  angles and modes, residual ≤ 0.5 deg). The de-embedding stopped ~1.24 mm short of the patch surface; corrected by
  multiplying every HFSS S by exp(+j2·k0·d·cos theta). (User: check the Floquet-port de-embed distance vs air height.)
* Cross terms: HFSS S(1,2) = S_TE_TM, S(2,1) = S_TM_TE, same sign as the MoM (sign +1 fits 3× better than −1).
* Circular terms rebuilt with the repo formulas (verified to 1e-14 against mom_reference.csv's own r columns):
  r_RR = ½(TMTM − TETE + j·TMTE + j·TETM), r_LR = ½(TMTM + TETE + j·TMTE − j·TETM), r_RL = ½(TMTM + TETE − j·TMTE + j·TETM).

## Result (HFSS vs MoM F=2n reference)
* The PB ladder works in HFSS: four states rotate the CP phase by ≈ 90 deg steps; low loss near resonance.
* **Best HFSS point is 10.5 GHz** (both angles): |r_RR| −0.55…−0.78 dB; steps 98/92/102/68 deg (θ=0),
  97/92/101/71 deg (θ=30); opposite-handed −16.8…−24.6 dB (θ=0), −11.9…−25.1 dB (θ=30).
* At 10.0 GHz HFSS is worse than the MoM: opposite-handed −8.4 dB (θ=0) and −5.7 dB (θ=30) vs MoM −13.5 / −14.4;
  |r_RR| down to −1.5 / −2.2 dB. Reading: the HFSS cell resonates ≈ 0.3–0.5 GHz higher than the MoM predicts.
* Phases: HFSS − MoM per state ≈ −12…−29 deg at 10 GHz; ladder max deviation from 90 deg: HFSS 21.9 deg (θ=0) /
  18.7 deg (θ=30) at 10.5 GHz vs MoM 12.8 / 9.8 deg. Persistent short step **state 4 → 1 ≈ 64–71 deg** (MoM 75–80):
  the axis-aligned (0, 90 deg) and diagonal (45, 135 deg) stubs are not equivalent on the square lattice.
* Amplitudes agree within ≈ 0.2 dB at 10.5 GHz; the model's loss picture is right.
