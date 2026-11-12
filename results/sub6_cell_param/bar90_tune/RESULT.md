# BAR90 island tuning — geometry hard stop

## STATUS

STOPPED before any AEDT build or solve. The Round 1 island-via rule is geometrically inconsistent with six of the seven requested configurations. Moving the via only to the island centre-line in y does not correct the x margin.

## GATES

| Gate | Result |
|---|---|
| G1, all per-row AEDT readbacks | NOT STARTED — geometry conflict before build |
| G2, clean SW branches | NOT STARTED |
| G3, T1–T4 ranking | NOT STARTED |

## DATA

The exported board's local island-via centre is **(x, y) = (−1.780, −0.530) mm** and pad radius is **0.300 mm** ([sub6_final_layout.py:29,35,74-77](../../../scripts/research/sub6_final_layout.py)). The island is centred in x, with end positions ±Lᵢ/2 ([sub6_final_rf_std.py:44,97-101](../../../scripts/research/sub6_final_rf_std.py)). Thus the via-pad margin to its nearest island end is `Lᵢ/2 − 1.780 − 0.300` mm. The requested minimum is 0.100 mm. The unrotated local frame is used below; rotating the whole cell does not change margins.

| Config | g mm | c mm | Lᵢ mm | Wᵢ mm | x margin mm | y margin at exported position mm | y centre-line mm | Valid after y move? |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Base | 0.30 | 0.60 | 4.348 | 2.00 | 0.094 | 0.020 | 0.150 | **No** |
| g0.60 | 0.60 | 0.60 | 4.348 | 1.70 | 0.094 | 0.020 | 0.000 | **No** |
| g1.00 | 1.00 | 0.60 | 4.348 | 1.30 | 0.094 | 0.020 | −0.200 | **No** |
| c0.40 | 0.30 | 0.40 | 4.348 | 2.20 | 0.094 | 0.220 | 0.050 | **No** |
| c0.80 | 0.30 | 0.80 | 4.348 | 1.80 | 0.094 | −0.180 | 0.250 | **No** |
| Lᵢ3.5 | 0.30 | 0.60 | 3.500 | 2.00 | −0.330 | 0.020 | 0.150 | **No** |
| Lᵢ5.2 | 0.30 | 0.60 | 5.200 | 2.00 | 0.520 | 0.020 | 0.150 | Yes |

The local transverse island bounds used for the y checks are `y_min = −1.45 + c` and `y_max = +1.45 − g` mm, from the fixed 2.9 mm slot width. The y margin at the exported via is `min(y_via−0.3−y_min, y_max−y_via−0.3)`; the y centre-line is `(y_min+y_max)/2`. These reproduce the exported base island bounds −0.85…+0.45 mm when g = 1.00 and c = 0.60 ([sub6_final_layout.py:30-31](../../../scripts/research/sub6_final_layout.py)).

## FLAGS

1. The base and every Lᵢ = 4.348 mm variant miss the x margin by **0.006 mm**. The Lᵢ = 3.5 mm variant's pad protrudes beyond the island by **0.330 mm**. Moving only y cannot repair either condition.
2. A possible common x coordinate is −1.34 mm: it leaves 0.11 mm at Lᵢ = 3.5 mm; every y can then be set to that island's centre-line. This changes the exported x position by 0.44 mm and has **not** been applied.
3. No AEDT project was opened, no geometry was built, and no solve was launched for this tuning job. No feed or board model was touched. No `.py` file is included in the output.

## FILES

`results/sub6_cell_param/bar90_tune/RESULT.md` and `results/sub6_cell_param/bar90_tune/geometry_gate.csv`.
