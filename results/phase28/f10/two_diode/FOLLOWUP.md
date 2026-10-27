# Two-diode cell redesign: mesh-64 follow-up (10.0 GHz, own loading, MODEL_ONLY)

Topology: rational-sector family, roles `MMGAMMBG` (diodes A, B; M = static metal neck, G = static gap), the best of the 3584-config mesh-32 screen. 'axis err' = largest deviation of a gap between neighbouring state axes from 45 deg; spreads = (max-min)/mean over the four states. Substrates: rt5880_1p575 is the board of all earlier phases; the others are LISTED options (none chosen).

| case | b_bar [mS] | delta [mS] | axis err [deg] | alpha [Np/m] | beta/k0 |
|---|---|---|---|---|---|
| T1_MMGAMMBG_rt5880_1p575 | 4.98..8.38 | 1.04..1.44 | 17.0 | 0.158..0.804 | 1.034..1.061 |
| T1_MMGAMMBG_rt5880_3p175 | 4.81..8.05 | 0.96..1.38 | 18.7 | 0.262..0.820 | n/a |
| T1_MMGAMMBG_ro4003c_1p524 | 8.70..15.93 | 1.84..3.34 | 11.8 | 0.652..12.443 | 1.173..1.396 |
| V1_hub1p0_MMGAMMBG | 5.31..9.34 | 1.25..1.66 | 12.0 | 0.343..0.865 | 1.036..1.072 |
| V2_diag3p0_MMGAMMBG | 6.80..11.26 | 0.75..2.40 | 20.3 | 0.207..1.597 | 1.047..1.098 |
| V3_hub1p0_diag3p0_MMGAMMBG | 7.87..13.43 | 1.22..3.52 | 37.9 | 0.317..4.639 | 1.056..1.129 |
| V4_neckw0p6_MMGAMMBG | 4.89..8.13 | 1.04..1.39 | 15.1 | 0.137..0.739 | 1.034..1.059 |
| LS_0.1nH_MMGAMMBG | 4.98..8.69 | 1.19..1.66 | 19.5 | 0.158..1.048 | 1.034..1.065 |
| LS_0.2nH_MMGAMMBG | 4.99..9.07 | 1.24..1.90 | 21.8 | 0.157..1.328 | 1.034..1.069 |
| LS_0.3nH_MMGAMMBG | 4.99..9.52 | 1.24..2.18 | 23.9 | 0.157..1.620 | 1.034..1.074 |

8-diode G (states 0 and 1) on the thicker boards, for reference:

| case | state | b_bar | delta | alpha | beta/k0 |
|---|---|---|---|---|---|
| REF_G8diode_rt5880_3p175 | 0 | 4.56 | 2.424 | 3.820 | nan |
| REF_G8diode_rt5880_3p175 | 1 | 4.30 | 2.154 | 3.309 | nan |
| REF_G8diode_ro4003c_1p524 | 0 | 7.51 | 4.803 | 26.628 | 1.135 |
| REF_G8diode_ro4003c_1p524 | 1 | 7.53 | 4.815 | 26.833 | 1.135 |
