# BAR90 island A — Step 0 audit

## STATUS

STOPPED before cell solves under the Step 0 code/table disagreement rule. The board cell and selected calibrated design method are identifiable: `sub6_final_design.py:16-18` defaults to `cal`, and `sub6_final_layout.py:368` calls it explicitly. The mesh-96 table's realized cell dimensions disagree numerically with the board exporter, and the directional HFSS calibration used FR4, bias stubs, and a different lower boundary. G2 would compare different physical models.

## GATES

| Gate | Result | Evidence |
|---|---|---|
| G1, geometry and reference readback | PASS | Requested numeric/source readbacks complete; the code/table discrepancy triggers the separate Step 0 stop rule. |
| G2, BAR50 anchor | NOT RUN | Step 0 ambiguity stop rule. |
| G3, clean branch | NOT RUN | No eigenmode solves. |
| G4, BAR90 differences | NOT RUN | No valid G2 anchor. |

## DATA

The board exporter reads `results/sub6_final_board/mesh96/fab/cell_table.csv` [sub6_final_rf_std.py:33,67]. All 240 switched rows identify `ISL_W2.9_q0.8_Li0.5_c0.6`, for example [cell_table.csv:12]. Its local board geometry comes from [sub6_final_rf_std.py:35-44,85-115,118-123,136-138,152-161] and [sub6_final_layout.py:21-31,67]:

| Item | Physical board cell, mm | Source |
|---|---|---|
| Period | 11.6 | `sub6_final_layout.py:21`; cell-table positions in multiples of 11.6 |
| Disc | radius 5.435; 35 µm Cu | `cell_table.csv:12`; `sub6_final_rf_std.py:35,85,94` |
| Slot | 8.696 × 2.9; local long axis x, rotated by `slot_axis_deg` on board | `cell_table.csv:12`; `sub6_final_rf_std.py:85-94` |
| Island | 4.348 × 1.3; local x = −2.174…2.174, y = −0.85…0.45 | `sub6_final_rf_std.py:44,97-101` |
| Gap/clearance | diode gap 1.00, island opposite-side clearance 0.60, end clearance 2.174 | `sub6_final_rf_std.py:45,87-100`; dimensions derived from those edges |
| Diode sheet | 1.00 × 0.35, local x = −0.175…0.175, y = 0.45…1.45, z = 0.0175; after +90° canonical rotation x = −1.45…−0.45, y = −0.175…0.175 | `sub6_final_rf_std.py:45,87-89,98,113`; `sub6_island_cell_packages.py:34-36,65-69` |
| Via | barrel radius 0.15, pad radius 0.30; local island (−1.78,−0.53), disc (−3.63,−3.63). After +90° rotation: (0.53,−1.78), (3.63,−3.63). Each board row rotates these positions by its own `slot_axis_deg`. | `sub6_final_layout.py:22,77-81`; `sub6_final_rf_std.py:42,102-112` |
| Anti-pad | radius 0.50 default, at each via through ground | `sub6_final_rf_std.py:42,118-123`; `sub6_final_layout.py:26` |
| Top substrate | RT5880 εr 2.2, tanδ 0.0009; z = −3.175…0, thickness 3.175 | `sub6_final_rf_std.py:5-6,35,153`; `sub6_planA_anchor_build.py:120-130`; V_OFF build readback `stage1_work/ANCHOR_V_OFF_build.json` (`top_substrate_before`, `material_readback`) |
| Ground | solid PEC Cu, z = −3.210…−3.175, 35 µm; via anti-pads and feed opening radius 2.05 | `sub6_final_rf_std.py:5-8,41-42,118-123,149`; `sub6_planA_anchor_build.py:124-130` |
| Bias layer for requested comparison | vacuum, z = −4.010…−3.175, thickness 0.835 in imported STEP bounding box, including the 0.035 ground band; nominal material layer below ground is 0.8, z = −4.010…−3.210. Via pad bottom is −4.045. | `sub6_final_rf_std.py:5-8,35-41,151`; `sub6_planA_anchor_build.py:45-49,115-130`; `stage1_work/ANCHOR_V_OFF_build.json` (`bias_after`, `material_readback`) |

The selected calibrated design stores **directional β/k₀**. Its input mesh-96 row separately stores an uncalibrated mean β/k₀:

| Reference | OFF | ON | Method/source |
|---|---:|---:|---|
| Mesh-96 switched-row mean β/k₀ at 5.8 GHz | 1.0970076923 | 1.0167711398 | `input_island_rows_n96_d789d9a.json:136,166`; read into board model by `sub6_final_cal.py:40-52` and `sub6_board_model4.py:48-65`. These values are marked `extrapolated=true` in the row. The realized MoM dimensions are slot L 8.7, island L 4.35, island W 1.329167, diode gap 0.966667, cell grid 0.120833 mm (`input_island_rows_n96_d789d9a.json:123-132`). This is not the exact physical board geometry. |
| HFSS-calibrated along β/k₀ at 5.8 GHz | 1.0032195493 | 0.9901423123 | `calibration.json:9,95`; `sub6_final_cal.py:56-113`. Along = phase y after +90° rotation (`sub6_island_cell_packages.py:147`). |
| HFSS-calibrated across β/k₀ at 5.8 GHz | 1.1465229371 | 1.0071263681 | `calibration.json:51,128`; `sub6_final_cal.py:56-113`. Across = phase x after +90° rotation (`sub6_island_cell_packages.py:147`). |

The HFSS calibration chooses a top surface-wave branch and rejects low-Q via/island and FR4 underside modes (`sub6_final_cal.py:3-7`, `hfss_data.json:10-15`). Its OFF points are 85° and 95° in each direction; ON uses only 95° in each direction because the 80° branch was ambiguous (`sub6_final_cal.py:3-7`, `hfss_data.json:16-23`, `calibration.json:18-28,60-70,104-107,137-140`). The MoM prior comprises seven frequency points 5.2–6.4 GHz per state/direction (`sub6_final_cal.py:34,56-76`; `input_pred_table_mom_nodata.json:1`), not a fixed eigenmode phase sweep. The calibrated board β is mixed by cos²/sin² and the ON-along value below 1 is floored to 1 (`sub6_final_cal.py:171-190`). BAR50 sheet values in the archived eigenmode package are OFF 5.424/−79.320 and ON 1.05/+7.65 Ω/sq (`sub6_island_cell_packages.py:147,188`). The mesh-96 MoM table itself used a frequency-dependent diode circuit (`sub6_iter1_island_retune.py:91-94`), so it is not a direct four-impedance-sheet eigenmode reference.

## FLAGS

1. **Reference geometry conflict:** `hfss_data.json:2` says the directional calibration used an FR4 0.8 mm underside, bias stubs, 30 mm air, and radiation bottom. The requested comparison is vacuum underside with no stubs. The eigenmode package also builds stubs and a radiation bottom for its via case (`sub6_island_cell_packages.py:40-69,141-149`). The RF board exporter has pads/barrels but no stubs (`sub6_final_rf_std.py:102-123`).
2. **Table hierarchy:** the board input row stores a single mean β/k₀, while the selected calibrated design uses directional β/k₀ and interpolates them (`sub6_final_design.py:16-18`, `sub6_final_layout.py:368`, `sub6_final_cal.py:171-190`). The selected method is clear, but its reference geometry does not match the requested vacuum cell.
3. **Discretization:** the mesh-96 row's realized geometry differs numerically from the physical STEP geometry as listed above.
4. No model was edited or retuned and no BAR50/BAR90 island cell solve was launched. Continuing past Step 0 requires an explicit decision to compare the vacuum board cell against the FR4/stub calibration despite the geometry mismatch, or a new vacuum-cell BAR50 reference.

## FILES

`results/sub6_cell_param/bar90_island_A/RESULT.md`; `results/sub6_cell_param/bar90_island_A/table_audit.csv`.
