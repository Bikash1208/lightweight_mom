# HFSS cost - K = 2 LP board (`hfss_exports/sidequest_board_k2_lp/`), ESTIMATE (2026-11-05, MODEL_ONLY)
**Nothing here is an HFSS result.** The repo records no tets / peak RAM / passes / time for the HFSS-validated Phase 26 150 mm single-feed board (`hex150_tensor_cp_single`, 1430 slotted circles, same stack RT5880 1.575 mm, 6.4 mm lattice): those are `NOT_AVAILABLE_FROM_CURRENT_DATA`. Phase 26's own order-of-magnitude note (`results/phase26/README.md` section L) gives "a few 10^7 tets, > 100 GB" for the 5684-cell 4-hexagon board (4-5x the single hexagon), i.e. about 0.7-1.5 x 10^7 tets for one hexagon. The second calibration is the measured 90 mm Phase 30 test board (572 cells, 5+ h; tets/RAM also NA).

| item | Phase 26 150 mm board | this board (K2 LP) |
|---|---|---|
| cells (6.0 mm discs) | 1430 slotted circles | 1430 rim discs, cross slit 0.4 x 5.6, rim 0.2 mm |
| lumped RLC sheets | 0 (1 port) | 720 (0.4 x 1.0, mesh 0.15 mm) + 1 port; 544 strapped cells |
| smallest feature | slot-end bridge 0.43 mm | rim 0.2 mm, slit 0.4 mm, sheet clearance 0.14 mm |
| STEP | 35 MB (1430 solids) | 97.6 MB (1432 metal lumps, 720 sheets) |
| tets (est.) | 0.7-1.5 x 10^7 (repo note, unmeasured) | **9-25 M** (x 1.3-1.7 finer features, + sheets, + 6 faces radiation box 14 mm over the probe) |
| RAM, direct (est., 10-15 GB per M tets) | NA | **100-300 GB** (lower end fits 256 GB; upper does not) |
| time, state A, 6 passes, delta S 0.02, one frequency | NA | **20-50 h** on 32-64 cores (from 90 mm: 5+ h x (1430/572)^1.3-1.7 = 17-24 h, x 1.2-2 for finer features) |
| state B (mesh reused) | | 7-16 h (about a third; a fresh mesh costs as much as A) |
| passes / delta S / basis | | 6 / 0.02 / first order (HFSS defaults to refine on the substrate and sheets) |
| licence | HFSS | HFSS + HPC Pack(s) for > 4 cores: 32+ cores needs 1-2 Packs; DDM only if RAM exceeds the node (total memory ~ direct) |

Levers if the upper bound is exceeded (all change the HFSS model, say so): drop the 0.15 mm sheet seeding to 0.25 mm; first-order basis is already the minimum; the 504-diode variant (analog statics, 13.6 dBi) removes 216 sheets but not a single cell, so it saves about 10 % only and needs un-validated analog static shapes - not recommended. Report tets, passes, final delta S, peak RAM and time of the first run so this table can be firmed up.
