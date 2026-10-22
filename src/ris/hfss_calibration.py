"""HFSS comparison and (approval-gated) calibration hook for the CP-2 state library (Phase 20 R4 item 6). MODEL_ONLY until the user's HFSS run exists.

The user's HFSS package `hfss_exports/phase20_r3_cp_cell/` asks for Floquet-port S11 (TE/TM, co and cross) of the four one-hot states at theta 0 and 30 deg, phi 0, 8-12 GHz. This module
* `read_hfss_rows` parses such an export in the same long format as the package's `mom_reference.csv` (columns `state, f_ghz, theta_deg, phi_deg, quantity, hfss_mag_dB, hfss_phase_deg`;
  `quantity` in S_TE_TE, S_TM_TM, S_TE_TM, S_TM_TE),
* `compare` evaluates the library at the same points and reports, per state and point, the r_RR phase / level difference and the opposite-handed difference, with a named verdict,
* `proposed_offsets` summarises the per-state differences as a CANDIDATE correction (phase offset, level offset) and does NOT apply it,
* `apply_calibration` returns a calibrated copy of the library arrays ONLY with `approved_by_user=True` and then labels the result `HFSS_CALIBRATED` with the data provenance; without approval it
  raises. No calibration exists until the user supplies data and approves it.

HFSS TM sign convention: the package asks the user to run the bare grounded slab first and to apply any 180 degree TM shift to every TM entry before comparing; `read_hfss_rows(tm_shift_180=True)`
applies it literally (S_TM_TM, S_TE_TM and S_TM_TE are multiplied by -1).
"""

import csv
import math
from dataclasses import dataclass

import numpy as np

from .polarization import circular_reflection, to_structure_frame

QUANTITIES = ("S_TE_TE", "S_TM_TM", "S_TE_TM", "S_TM_TE")
VERDICT_AGREE = "HFSS_AGREES_WITHIN_LIBRARY_UNCERTAINTY"
VERDICT_DIFFER = "HFSS_DIFFERS_FROM_LIBRARY"


def read_hfss_rows(path, tm_shift_180=False):
    """{(state, f_ghz, theta_deg, phi_deg): 2x2 complex (TE, TM) reflection} from a long-format CSV (see module docstring). `tm_shift_180` applies the 180 degree TM shift to every TM entry."""
    entries = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["quantity"] not in QUANTITIES:
                continue
            key = (
                int(float(r["state"])),
                round(float(r["f_ghz"]), 6),
                float(r["theta_deg"]),
                float(r.get("phi_deg", 0.0)),
            )
            z = 10 ** (float(r["hfss_mag_dB"]) / 20) * np.exp(1j * math.radians(float(r["hfss_phase_deg"])))
            m = entries.setdefault(key, np.zeros((2, 2), dtype=complex))
            i, j = {"S_TE_TE": (0, 0), "S_TM_TM": (1, 1), "S_TE_TM": (0, 1), "S_TM_TE": (1, 0)}[r["quantity"]]
            sign = -1 if (tm_shift_180 and (i == 1 or j == 1)) else 1
            m[i, j] = sign * z
    return entries


@dataclass
class PointComparison:
    state: int
    f_ghz: float
    theta_deg: float
    phi_deg: float
    rr_phase_diff_deg: float
    rr_level_diff_db: float
    opp_diff_db: float
    within_uncertainty: bool


def compare(library, hfss, phase_floor_deg=1.0, db_floor=0.1):
    """Compare library (corrected) with HFSS at the HFSS points. A point `within_uncertainty` if the r_RR difference is covered by the library's stored uncertainty
    (converted to degrees / dB) with floors `phase_floor_deg`, `db_floor` so a tiny stored uncertainty is not demanded to be exact. Returns (list of PointComparison, summary dict)."""
    out = []
    for (state, f, th, ph), s_h in sorted(hfss.items()):
        c_h = to_structure_frame(circular_reflection(s_h), math.radians(ph))
        c_l = library.structure_circular(f, th, ph, state)[0]
        unc = library.uncertainty(f, th, ph, state)[0]
        r_h, r_l = c_h[0, 0], c_l[0, 0]
        dphase = math.degrees(np.angle(r_h / r_l))
        ddb = 20 * math.log10(abs(r_h) / abs(r_l))
        opp_h = 20 * math.log10(max(abs(c_h[1, 0]), abs(c_h[0, 1])))
        opp_l = 20 * math.log10(max(abs(c_l[1, 0]), abs(c_l[0, 1])))
        u_ph = max(phase_floor_deg, math.degrees(unc[0, 0] / abs(r_l)))
        u_db = max(db_floor, 8.686 * unc[0, 0] / abs(r_l))
        out.append(
            PointComparison(state, f, th, ph, dphase, ddb, opp_h - opp_l, abs(dphase) <= u_ph and abs(ddb) <= u_db)
        )
    ph = np.array([abs(p.rr_phase_diff_deg) for p in out])
    db = np.array([abs(p.rr_level_diff_db) for p in out])
    ok = np.array([p.within_uncertainty for p in out])
    summary = {
        "n_points": len(out),
        "phase_diff_deg": {"median": float(np.median(ph)), "p95": float(np.percentile(ph, 95)), "max": float(ph.max())}
        if len(out)
        else None,
        "level_diff_db": {"median": float(np.median(db)), "p95": float(np.percentile(db, 95)), "max": float(db.max())}
        if len(out)
        else None,
        "fraction_within_uncertainty": float(ok.mean()) if len(out) else None,
        "verdict": VERDICT_AGREE if len(out) and ok.all() else VERDICT_DIFFER,
    }
    return out, summary


def proposed_offsets(comparisons):
    """Per-state median r_RR phase and level differences (HFSS minus library): a CANDIDATE correction, not applied."""
    res = {}
    for state in sorted({p.state for p in comparisons}):
        sel = [p for p in comparisons if p.state == state]
        res[state] = {
            "phase_offset_deg": float(np.median([p.rr_phase_diff_deg for p in sel])),
            "level_offset_db": float(np.median([p.rr_level_diff_db for p in sel])),
            "n_points": len(sel),
        }
    return res


class CalibrationNotApproved(RuntimeError):
    """Raised when `apply_calibration` is called without the user's explicit approval."""


def apply_calibration(arrays, offsets, *, approved_by_user=False, source=""):
    """Return a calibrated copy of the library arrays (`c_corr` and `s_corr` only; `c_raw` untouched) with the per-state r_RR phase / level offsets applied, labelled
    `HFSS_CALIBRATED`. Refuses unless `approved_by_user=True` (the offsets are the user's decision, from data that must exist)."""
    if not approved_by_user:
        raise CalibrationNotApproved(
            "HFSS calibration is applied only with the user's explicit approval (approved_by_user=True)"
        )
    if not source:
        raise ValueError("a calibration needs its data provenance (`source`), e.g. the HFSS export file name and date")
    out = {k: np.array(v) for k, v in arrays.items()}
    phi = out["phi_deg"]
    for state, off in offsets.items():
        g = 10 ** (off["level_offset_db"] / 20) * np.exp(1j * math.radians(off["phase_offset_deg"]))
        out["c_corr"][:, :, :, state - 1, 0, 0] *= g
        out["c_corr"][:, :, :, state - 1, 1, 1] *= np.conj(
            g
        )  # the L-handed term mirrors the R-handed one (a rotation reversal): same level, opposite phase
    from .state_library import te_tm_from_structure_circular

    out["s_corr"] = te_tm_from_structure_circular(out["c_corr"], phi)
    out["calibration"] = np.array([f"HFSS_CALIBRATED: {source}"])
    return out
