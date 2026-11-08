"""Writes results/sub6_board_model/BOARD_MODEL.md from points.csv and provenance_run.json (MODEL_ONLY)."""
import csv
import json
from pathlib import Path

D = Path(__file__).resolve().parents[2] / "results" / "sub6_board_model"
R = list(csv.DictReader(open(D / "points.csv")))
prov = json.load(open(D / "provenance_run.json"))
n_cells = {int(float(r["side_mm"])): int(float(r["N_cells"])) for r in R}


def f(r, k, d=2):
    return f"{float(r[k]):.{d}f}"


t1 = "| side mm | N_sw (f) | pl | OFF G_A | OFF G_B | balance | OFF AR A/B | ON G_A | B supp | ON AR_A | SLL off/on | gate (margin) |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
t1 += "".join(f"| {int(float(r['side_mm']))} | {r['N_sw']} ({f(r, 'f_sw')}) | {r['placement']} | {f(r, 'off_GA', 1)} | {f(r, 'off_GB', 1)} | {f(r, 'off_balance')} | {f(r, 'off_AR_A')}/{f(r, 'off_AR_B')} | {f(r, 'on_GA', 1)} | {f(r, 'on_suppression', 1)} | {f(r, 'on_AR_A')} | {f(r, 'sll_off', 1)}/{f(r, 'sll_on', 1)} | {r['gate']} ({f(r, 'worst_margin_db')}) |\n" for r in R)
t2 = "| side mm | N_sw | pl | phase error ON-OFF at static sites: rms / max / mean deg | switched d_beta (ON-OFF) | beam-A gain cost (dB) | ON AR_A with / without | ON SLL with / without |\n|---|---|---|---|---|---|---|---|\n"
t2 += "".join(f"| {int(float(r['side_mm']))} | {r['N_sw']} | {r['placement']} | {f(r, 'dphi_rms_deg', 1)} / {f(r, 'dphi_max_deg', 1)} / {f(r, 'dphi_mean_deg', 1)} | {f(r, 'sw_beta_on_minus_off_mean', 3)} | {f(r, 'A_loss_from_beta_db')} | {f(r, 'on_AR_A')} / {f(r, 'on_AR_A_no_beta_change')} | {f(r, 'sll_on', 1)} / {f(r, 'sll_on_no_beta_change', 1)} |\n" for r in R)
t3 = "| side mm | N_sw | alpha0 mean (Np/m) | cap1 | q90 demand | eta OFF / ON | alpha-saturated static / switched sites | switched amplitude error rms % |\n|---|---|---|---|---|---|---|---|\n"
t3 += "".join(f"| {int(float(r['side_mm']))} | {r['N_sw']} | {f(r, 'alpha0_mean')} | {f(r, 'cap1_Np_m')} | {f(r, 'q90_demand', 1)} | {f(r, 'eta_off')} / {f(r, 'eta_on')} | {r['st_sites_alpha_saturated']} / {r['sw_sites_alpha_saturated']} | {f(r, 'sw_amp_err_rms_pct', 2)} |\n" for r in R if r["placement"] == "a")
n_meet = sum(r["gate"] == "MEETS_ALL" for r in R)
best = max(R, key=lambda r: float(r["worst_margin_db"]))
same = 0
for side, ns in ((200, (150, 200, 250, 300)), (150, None)):
    for r in R:
        pass
pairs = {}
for r in R:
    pairs.setdefault((r["side_mm"], r["N_sw"]), {})[r["placement"]] = r
same = sum(1 for v in pairs.values() if abs(float(v["a"]["worst_margin_db"]) - float(v["b"]["worst_margin_db"])) < 1e-9 and v["a"]["off_GA"] == v["b"]["off_GA"])
md = f"""# SUB6 board model with real cell tables (stage: sub6-board-model)

**MODEL_ONLY / PENDING_HFSS_CELL_CHECK / NOT_HFSS_VALIDATED.** Every number is a model number for cells whose tensors come from `results/sub6_slotdisc/rows_RT5880.jsonl` (mesh 64, MoM); no cell has been checked in HFSS yet. Job `jobs/sub6-board-model.md`; branch `compute/sub6-board-model`; data `points.csv`, `layout_<side>mm_N<n>_<a|b>.csv` (one per point: site, x, y, STATIC/SWITCHED, row id, slot axis, alpha and beta OFF/ON, required alpha, ON-OFF phase), `provenance_run.json`; code `scripts/research/sub6_board_model.py` (+ `sub6_board_req.py` copied unmodified from compute/sub6-board-requirements 7694971).

## Setup
5.8 GHz, RT5880 3.175 mm on PEC, p = 11.6 mm, hexagon side 200 mm ({n_cells[200]} cells) and 150 mm ({n_cells[150]} cells), centre probe, **launch fraction {prov['launch_fraction']:.4f}** (`src/dispersion/pin_launch.py` at 5.8 GHz; the Phase 26 pin scaled geometrically by the substrate-height ratio 3.175/1.575, as in stage 1b), RHCP beams A (25, 0) and B (25, 180). Phase 26 `march_rays` / `cell_design` unmodified; no frozen code edited.

**Architecture.** OFF: static cells radiate A+B (one linear current, amplitude via the cell alpha, axis psi). ON: static cells unchanged, switched cells radiate -B (near-silent in OFF: alpha_OFF <= 0.1 Np/m). Group-normalised weights: static (uA+uB)/(1-f), switched -uB/f, f = N_sw/N.
**Cells.** Static = NOCUT rows, switched = notch rows (OFF_BAR50 / ON_BAR50 with their own alpha and beta). Of {prov['rows_total']} rows {prov['rows_resonance_flagged']} are resonance-flagged (proxy) and excluded; the static filter beta_NOCUT <= 1.2 leaves **{prov['static_rows']} static rows**; the switched filter alpha_OFF <= 0.1, alpha_ON >= 0.3, beta_ON <= 1.2 leaves **{prov['switched_rows']} switched rows** (the filters are mine). Each site takes the row nearest in sqrt(alpha) to the required alpha. Rows differ in gap and W, so neighbouring cells differ in disc size; coupling between different cells is not modelled.
**Real beta.** Per state, along every ray, the cell under the ray gives its own alpha and beta (static beta_NOCUT 1.05-1.2; switched OFF 1.04-1.06, ON 1.06-1.19). The slot axes are designed with the OFF phase; ON keeps the axes and takes the ON wave phase.
**Placement.** (a) low-discrepancy interleave; (b) the ideal-cell projection optimiser of stage 1b (rho = -1000, the limit of this architecture), then evaluated with real cells. **(b) returned the same result as (a) at {same} of {len(pairs)} points**: the optimiser found no improvement there.
**Sweep.** N_sw 150/200/250/300 at 200 mm (f 0.19-0.39); 87/117/146/175 at 150 mm (same fractions). {prov['seconds']:.0f} s for 16 points (the 16 tasks limit the pool to 16 workers; the machine has 64).

## Gates (targets: balance <= 3 dB, AR <= 3 dB (OFF A, B and ON A), B suppression >= 12 dB, SLL <= -13 dB, both states)
**MEETS_ALL: {n_meet} of {len(R)} points.** Balance and suppression are met everywhere. The misses are AR in ON (beam A built from the static lattice with the switched sites removed) and SLL in ON. Best margin {float(best['worst_margin_db']):.2f} dB ({int(float(best['side_mm']))} mm, N {best['N_sw']}, placement {best['placement']}, {best['gate']}).

{t1}
Gain is referred to the accepted power x the model conversion eta, so the absolute G_A are model numbers; ON G_A is higher than OFF because the switched cells add conversion in ON. B suppression is the coherent cancellation of discrete real cells (switched amplitude error rms <= 0.4 %, the table is dense) and is the most optimistic number here.

## Beam-A phase error from the switched cells' beta change (ON - OFF)
{t2}
The switched cells' beta rises by about 0.06 in ON (mean over the chosen rows). Downstream static cells see a wave phase shifted by up to 57 deg (rms 7-16 deg); the cost is **<= 0.14 dB of beam-A gain and <= 0.3 dB in AR / SLL** (columns "without" = the ON amplitudes with the OFF phase). The beta change is therefore not what makes the ON AR / SLL miss: the miss is present without it.

## Design bookkeeping (placement a)
{t3}
alpha0 = mean alpha of the base march (engine, static rows with alpha <= cap1 = min(static, switched max alpha) / q90 of the demand, the alpha cap of stage 1b). The switched cells need 1/f of the amplitude, so the demand q90 is 3.1-10.4 and the base alpha (hence the conversion) is small at small N_sw. No site is alpha-saturated.

## Limits
Cells are mesh-64 MoM, not HFSS-checked (the HFSS eigenmode package is `hfss_exports/sub6_slotdisc_cell/`, 2 of the 80 rows). beta/alpha of a row are those of the infinite periodic array of that row; the ON tensor is the BAR50 model at 5.8 GHz. The A+B current of a static cell is treated as one ideal linear current whose amplitude is the chosen row's alpha-amplitude. The switched OFF polarisation is assumed equal to the ON direction. The placement (b) objective is the ideal-cell one. Stopped here: no cell choice, no board design, no HFSS package.
"""
(D / "BOARD_MODEL.md").write_text(md, encoding="utf-8", newline="\n")
print(n_meet, same, len(pairs))
