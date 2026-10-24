"""Cell assignment for the dual-beam lens: library used only as the initial guess, every cell then solved
directly at its own (theta, phi) with the full 2x2 Jones response (no table or angle interpolation).

Family rule: both families keep their full printable radius range (pillar R in {0} + [R_min, R_max]; hole
nylon fill in [1 - fill(R_max), cap]).  Their phase ranges overlap, so the family is chosen by the target's *delay
position* d (0 = the empty cell, growing with nylon): pillar for d <= ``d_seam_deg``, hole above, falling back to the
other family when the preferred one cannot reach the target.
"""

from __future__ import annotations

import csv
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from .cells import CellFamilyModel, jones_xy, solve_radius_for_phase, wrap_deg

SQRT8 = 2.0 * math.sqrt(2.0)


def radius_for_fill(fill_octagon: float, period_m: float) -> float:
    """Circumradius whose octagon covers ``fill_octagon`` of the cell (octagon area = 2 sqrt2 R^2)."""
    return period_m * math.sqrt(fill_octagon / SQRT8)


@dataclass(frozen=True)
class FamilyLimits:
    """Allowed radius intervals [m] per family (min feature, seam and fill cap applied)."""

    pillar_min_m: float
    pillar_max_m: float
    hole_min_m: float
    hole_max_m: float

    @staticmethod
    def build(period_m, min_feature_m, cap_fill):
        c = math.cos(math.radians(22.5))
        r_min = min_feature_m / (2 * c)
        r_max = (period_m - min_feature_m) / (2 * c)
        return FamilyLimits(
            pillar_min_m=r_min,
            pillar_max_m=r_max,
            hole_min_m=radius_for_fill(1.0 - cap_fill, period_m),
            hole_max_m=r_max,
        )


def load_guess_library(path):
    """Rows of the committed library CSV as dicts with float fields and a complex 2x2 ``t`` (s/p)."""
    rows = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            t = np.array(
                [
                    [
                        complex(float(r["t_ss_re"]), float(r["t_ss_im"])),
                        complex(float(r["t_sp_re"]), float(r["t_sp_im"])),
                    ],
                    [
                        complex(float(r["t_ps_re"]), float(r["t_ps_im"])),
                        complex(float(r["t_pp_re"]), float(r["t_pp_im"])),
                    ],
                ]
            )
            rows.append(
                dict(
                    kind=r["kind"],
                    R=float(r["radius_mm"]) * 1e-3,
                    fill=float(r["fill"]),
                    theta=float(r["theta_deg"]),
                    phi=float(r["phi_deg"]),
                    t=t,
                )
            )
    return rows


def fold_phi(phi_deg):
    """Fold an azimuth into [0, 45] using the cell's mirror / C4 symmetry."""
    f = abs(phi_deg) % 90.0
    return 90.0 - f if f > 45.0 else f


def guess_table(rows, cap_fill):
    """Initial-guess candidates allowed by the family rule: {phi_lib: (kind[], R[], theta[], phase[nc, nth], loss[nc, nth])}.

    Phase is the co-polar Jones phase J_xx (deg) at the library's own (theta, phi); unwrapped along theta.
    """
    by = {}
    for r in rows:
        if r["kind"] == "hole" and r["fill"] > cap_fill + 1e-9:
            continue
        j = jones_xy(r["t"], r["theta"], r["phi"])
        by.setdefault(r["phi"], {}).setdefault((r["kind"], r["R"]), {})[r["theta"]] = (
            math.degrees(np.angle(j[0, 0])),
            20 * math.log10(abs(j[0, 0])),
        )
    out = {}
    for phi, cands in by.items():
        keys = sorted(cands)
        thetas = sorted({t for c in cands.values() for t in c})
        keys = [k for k in keys if all(t in cands[k] for t in thetas)]
        ph = np.array([[cands[k][t][0] for t in thetas] for k in keys])
        ls = np.array([[cands[k][t][1] for t in thetas] for k in keys])
        ph = np.degrees(np.unwrap(np.radians(ph), axis=1))
        out[phi] = (
            np.array([k[0] for k in keys], dtype=object),
            np.array([k[1] for k in keys]),
            np.array(thetas),
            ph,
            ls,
        )
    return out


def cell_guess(table, theta_deg, phi_deg):
    """Guess arrays (kind[], R[], phase_deg[], loss_db[]) for one cell: nearest library azimuth, theta-interpolated.

    The interpolation exists only to start the direct solve close to its answer and to score psi0; the final
    radius of every cell is a direct solve at its own (theta, phi).
    """
    phis = np.array(sorted(table))
    key = float(phis[np.argmin(np.abs(phis - fold_phi(phi_deg)))])
    kd, rr, ths, ph, ls = table[key]
    t = float(np.clip(theta_deg, ths[0], ths[-1]))
    i1 = int(np.clip(np.searchsorted(ths, t), 1, len(ths) - 1))
    w = (t - ths[i1 - 1]) / (ths[i1] - ths[i1 - 1])
    return kd, rr, ph[:, i1 - 1] * (1 - w) + ph[:, i1] * w, ls[:, i1 - 1] * (1 - w) + ls[:, i1] * w


def pick_candidate(guess, target_deg, d_seam_deg, tol_deg=10.0):
    """Index of the guess candidate for one cell under the family rule (see module docstring)."""
    kd, _rr, ph, ls = guess
    empty = np.flatnonzero((kd == "pillar") & (_rr == 0.0))
    ph_empty = ph[empty[0]] if empty.size else float(ph.max())
    d_target = (ph_empty - target_deg) % 360.0
    err = np.abs(wrap_deg(ph - target_deg))
    ok = err <= tol_deg
    pref = "pillar" if d_target <= d_seam_deg else "hole"
    for fam in (pref, "hole" if pref == "pillar" else "pillar"):
        m = ok & (kd == fam)
        if m.any():
            return int(np.argmax(np.where(m, ls, -np.inf))), True
    return int(np.argmin(err)), False


def choose_psi0(guesses, base_phase_deg, weight, d_seam_deg, step=2.0):
    """Global phase offset maximising the illumination-weighted mean |T|^2 of the rule-selected guesses."""
    best = None
    for psi0 in np.arange(0.0, 360.0, step):
        tgt = base_phase_deg + psi0
        loss = np.empty(len(tgt))
        feasible = True
        for i, g in enumerate(guesses):
            j, ok = pick_candidate(g, tgt[i], d_seam_deg)
            if not ok:
                feasible = False
                break
            loss[i] = g[3][j]
        if not feasible:
            continue
        w = float(np.sum(weight * 10 ** (loss / 10)) / np.sum(weight))
        if best is None or w > best[1]:
            best = (float(psi0), w)
    return best


def _job(a):
    model, kind, target, th, ph, r0, lo, hi = a
    return solve_radius_for_phase(model, kind, target, th, ph, r0, lo, hi)


def assign_cells(model: CellFamilyModel, limits: FamilyLimits, guesses, th, phi, target_deg, d_seam_deg, workers=4):
    """Direct solve for every cell; returns list of dicts (kind, R, residual_deg, J, n_solves)."""
    jobs, meta = [], []
    for i in range(len(th)):
        kd, rr, _ph, _ls = guesses[i]
        j, _ok = pick_candidate(guesses[i], target_deg[i], d_seam_deg)
        kind, r0 = str(kd[j]), float(rr[j])
        lo, hi = (
            (limits.pillar_min_m, limits.pillar_max_m) if kind == "pillar" else (limits.hole_min_m, limits.hole_max_m)
        )
        jobs.append((model, kind, float(target_deg[i]), float(th[i]), float(phi[i]), r0, lo, hi))
        meta.append((kind, r0))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(_job, jobs, chunksize=4))
    out = []
    for (kind, r0), (r, resid, _jxx, ns, jfull), job in zip(meta, res, jobs, strict=True):
        out.append(dict(kind=kind, R=r, residual_deg=resid, J=jfull, n_solves=ns, R_guess=r0, bounds=(job[6], job[7])))
    return out


def resolve_empty_cells(model, out, th, phi, target_deg, pillar_min_m):
    """Pillar cells whose best R sits at the pillar_min bound may be closer to the empty cell (R = 0): take the better."""
    for i, c in enumerate(out):
        if c["kind"] == "pillar" and c["R"] <= pillar_min_m * (1 + 1e-6):
            j0 = model.jones("pillar", 0.0, th[i], phi[i])
            e0 = float(wrap_deg(math.degrees(np.angle(j0[0, 0])) - target_deg[i]))
            if abs(e0) < abs(c["residual_deg"]):
                c.update(R=0.0, residual_deg=e0, J=j0)
    return out
