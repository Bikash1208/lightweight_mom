"""Phase 19 items 2+3: through-hole vs blind-pocket hole family, and base-plate thickness sweep.

    .venv/bin/python scripts/research/phase19_plate_study.py

Total height fixed at 25 mm; plate in {0.6, 1, 2, 3.42} mm.  Per configuration a coarse library
(theta 0/12/25 deg, phi = 0, p-pol) is built; every lens cell is assigned the lowest-loss R that hits its
required phase (family seam at fill f_seam), and configurations are ranked by the illumination-weighted
aperture-averaged |S21|^2.  The theta interpolation here is for RANKING only; the final design (Part C)
re-solves every cell directly at its own (theta, phi).
Writes results/nylon_cells_28ghz/plate_study.json.
"""

import json
import math
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.dielectric_cell import OctagonPillarCell, solve_cell  # noqa: E402
from src.dielectric_cell.library import COS_FLAT  # noqa: E402
from src.transmitarray.design_inputs import dual_beam_phase_map  # noqa: E402

F, P, H_TOT = 28e9, 5.0e-3, 25.0e-3
EPS, TAND, MIN_F = 2.45, 0.0185, 0.6e-3
THETAS = (0.0, 12.0, 25.0)
N_R = 28
R_MIN, R_MAX = MIN_F / (2 * COS_FLAT), (P - MIN_F) / (2 * COS_FLAT)
SPAN_MARGIN = 366.0
PLATES_MM = (0.6, 1.0, 2.0, 3.42)


def _job(a):
    kind, r, plate, h, th = a
    c = OctagonPillarCell(P, r, h, EPS, TAND, plate, 22.5, kind)
    t = solve_cell(c, F, th, 0.0, 8).t[1, 1]
    return (kind, r, plate, th, c.fill_fraction, complex(t))


def family_curves(rows, kind, plate):
    """{theta: (R array, fill, unwrapped phase deg, loss dB)} for one family at one plate thickness."""
    out = {}
    for th in THETAS:
        sel = sorted((r for r in rows if r[0] == kind and r[2] == plate and r[3] == th), key=lambda r: r[1])
        R = np.array([r[1] for r in sel])
        fill = np.array([r[4] for r in sel])
        ph = np.degrees(np.unwrap(np.angle([r[5] for r in sel])))
        loss = np.array([20 * math.log10(abs(r[5])) for r in sel])
        out[th] = (R, fill, ph, loss)
    return out


def cap_fill_for_span(pillar0, hole0):
    """Smallest fill cap giving >= SPAN_MARGIN deg of phase span at normal incidence (families joined at the seam)."""
    Rp, fp, pp, _ = pillar0
    Rh, fh, ph, _ = hole0
    sel_p = np.ones(len(fp), bool)  # both families over their full printable range (as in the final design)
    sel_h = np.ones(len(fh), bool)
    fill = np.concatenate([fp[sel_p], fh[sel_h]])
    # unwrap along fill across the joined curve
    wrapped = np.radians(np.concatenate([pp[sel_p], ph[sel_h]]))
    order = np.argsort(fill)
    fill = fill[order]
    phase = np.degrees(np.unwrap(wrapped[order]))
    span = np.abs(phase - phase[0])
    idx = np.nonzero(span >= SPAN_MARGIN)[0]
    return (float(fill[idx[0]]), float(span.max())) if idx.size else (None, float(span.max()))


def theta_blend(th_cell):
    """Linear-in-theta weights (i0, i1, w1) of each cell between the THETAS samples."""
    ths = np.array(THETAS)
    i1 = np.clip(np.searchsorted(ths, th_cell), 1, len(ths) - 1)
    i0 = i1 - 1
    w1 = np.clip((th_cell - ths[i0]) / (ths[i1] - ths[i0]), 0.0, 1.0)
    return i0, i1, w1


def family_cell_tables(cur, fam, cap, th_cell):
    """(phase[Ncell, nR], loss[Ncell, nR]) on a dense R grid, restricted to the family's fill range."""
    R0 = cur[0.0][0]
    Rd = np.linspace(0.0 if fam == "pillar" else R0.min(), R0.max(), 500)
    fl = np.interp(Rd, R0, cur[0.0][1])
    ok = np.ones(len(fl), bool) if fam == "pillar" else (fl <= cap)
    ph = np.array([np.interp(Rd[ok], cur[t][0], cur[t][2]) for t in THETAS])
    ls = np.array([np.interp(Rd[ok], cur[t][0], cur[t][3]) for t in THETAS])
    i0, i1, w1 = theta_blend(th_cell)
    return (ph[i0] * (1 - w1)[:, None] + ph[i1] * w1[:, None], ls[i0] * (1 - w1)[:, None] + ls[i1] * w1[:, None])


def main():
    jobs = []
    radii = np.linspace(R_MIN, R_MAX, N_R)
    # through-hole family: solid 25 mm, pocket through (plate irrelevant -> plate 0)
    for th in THETAS:
        for r in [0.0, *radii]:
            jobs.append(("hole", float(r), 0.0, H_TOT, th))
    for tb in PLATES_MM:
        t = tb * 1e-3
        for th in THETAS:
            for kind in ("pillar", "hole"):
                for r in [0.0, *radii]:
                    jobs.append((kind, float(r), t, H_TOT - t, th))
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(_job, jobs, chunksize=8))

    x, y, th_cell, weight, target_base = dual_beam_phase_map()
    results = {}
    for tb in PLATES_MM:
        t = tb * 1e-3
        pill = family_curves(rows, "pillar", t)
        for mode in ("blind", "through"):
            hole = family_curves(rows, "hole", t if mode == "blind" else 0.0)
            cap, span = cap_fill_for_span(pill[0.0], hole[0.0])
            if cap is None:
                results[f"plate{tb}_{mode}"] = {"plate_mm": tb, "hole": mode, "max_span_deg": span, "feasible": False}
                continue
            tabs = [family_cell_tables(pill, "pillar", cap, th_cell), family_cell_tables(hole, "hole", cap, th_cell)]
            best = None
            gap_best = (1e9, None)
            for psi0 in np.arange(0.0, 360.0, 3.0):
                tgt = np.mod(target_base + psi0, 360.0)
                cell_loss = np.full(len(x), -99.0)
                miss = np.full(len(x), 1e9)
                for ph, ls in tabs:
                    d = np.abs((ph - tgt[:, None] + 180.0) % 360.0 - 180.0)
                    miss = np.minimum(miss, d.min(axis=1))
                    d = np.where(d < 3.0, d, np.inf)
                    cand = np.where(np.isfinite(d), ls, -99.0).max(axis=1)
                    cell_loss = np.maximum(cell_loss, cand)
                if miss.max() < gap_best[0]:
                    gap_best = (float(miss.max()), float(psi0))
                if np.all(cell_loss > -50):
                    wmean = float(np.sum(weight * 10 ** (cell_loss / 10)) / np.sum(weight))
                    if best is None or wmean > best[1]:
                        best = (float(psi0), wmean, float(np.mean(cell_loss)), float(cell_loss.min()))
            results[f"plate{tb}_{mode}"] = {
                "plate_mm": tb,
                "hole": mode,
                "fill_cap": cap,
                "max_span_deg": span,
                "feasible": best is not None,
                "best_worst_cell_phase_miss_deg": gap_best[0],
                "psi0_of_best_miss_deg": gap_best[1],
                **(
                    {
                        "psi0_deg": best[0],
                        "aperture_avg_power_transmission": best[1],
                        "aperture_avg_loss_db": 10 * math.log10(best[1]),
                        "mean_cell_loss_db": best[2],
                        "worst_cell_loss_db": best[3],
                    }
                    if best
                    else {}
                ),
            }
            print(tb, mode, results[f"plate{tb}_{mode}"], flush=True)
    out = {
        "notes": "ranking study, theta-interpolated p-pol phi=0 library; NOT the final per-cell solve",
        "family_rule": "both families over their full printable range; lowest-loss candidate within 3 deg of the target",
        "span_margin_deg": SPAN_MARGIN,
        "results": results,
    }
    (ROOT / "results" / "nylon_cells_28ghz" / "plate_study.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
