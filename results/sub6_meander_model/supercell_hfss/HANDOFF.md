# Rotated-supercell leakage: calibration stop

Status: **STOPPED before the supercell**. All HFSS results here used a fresh local
AEDT 2025.2 session. The user authorized v252 after v251 was unavailable. Work
used only the `compute/sub6-meander-model` branch.

## Step 0 — open-top support (HFSS)

- `AssignRadiation` on the top face of a bridged D3p6 Eigenmode cell failed
  before a solve: AEDT reported `Radiation boundaries are not used when the
  solution type is Eigenmode`.
- `CreatePML` on that top face created `PML_AIR1`; a subsequent solve returned
  finite Q. At phase 85° across, the first four modes were 5.55 GHz / Q 930.06,
  5.81 / 1846.18, 6.09 / 363308.36, and 6.65 / 791.81. This establishes API
  and solver support, not a leakage measurement or mode assignment.
- The initial PML probe failed a PyAEDT boundary-list check because PML appears
  as a cover object, not a normal boundary in that list. A second fresh build
  read back `PML_AIR1` and solved. This was a readback-code correction, not a
  model change.

## Step 1 — lidded 85° across calibration (HFSS)

| Cell | Reference GHz | Local GHz | Error | Gate |
|---|---:|---:|---:|---|
| Anchor `ST_R5p45_L8p72_W2p91` | 5.5531686625 | 5.5623108979 | +0.164631% | Pass |
| `MN_R5p6_W0p4_N2_T1p2_D3p0_E_B2p2` | 5.6932301269 | 5.6932301269 | <0.000001% | Pass |
| `MN_R5p6_W0p4_N2_T1p2_D3p6_E_B2p2` | 5.6825625491 | 5.6825625491 | <0.000001% | Pass |
| `MN_R5p6_W0p4_N2_T1p2_D3p6_E` | 5.6124610222 | — | — | **Stopped at copper readback** |

The open-twin cut creates disconnected copper pieces. The builder's original
`DISC` handle returned zero volume. Summing all PEC solids also returned zero
on a single retry. The expected total is 87.08 mm². Because this is a geometry
readback failure, the open row was not solved and no rotated supercell was built.
Do not infer an alpha from the open-top unit-cell Q values.

## HANDOFF

- STATUS: STOPPED, geometry readback on open twin; no supercell results.
- CELL JSON: `SUPERCELL_SPEC.md` at 242f3fb is the source of cell parameters.
- β/k0 along/across: not measured in this run.
- α: not measured; no gate call.
- M_t: not measured.
- Loss ratio: not measured.
- Evidence: three lidded HFSS calibration solves; one PML unit-cell probe.
- Next request: verify the split-body copper readback using an AEDT object/face
  audit and rerun the open-twin calibration before starting the anchor supercell.
