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


C0 = 299792458.0


def _reim(row, i, j):
    return float(row[f"re(S(FloquetPort1:{i},FloquetPort1:{j})) []"]) + 1j * float(
        row[f"im(S(FloquetPort1:{i},FloquetPort1:{j})) []"]
    )


def read_hfss_reim(path, plane_offset_mm=0.0, phi_deg=0.0):
    """Same dict as `read_hfss_rows` from an HFSS 'scan_theta, state, Freq, re/im S(FloquetPort1:i,FloquetPort1:j)' export (mode 1 = TE, mode 2 = TM; S(i,j) = received i, incident j).
    `plane_offset_mm` de-embeds the Floquet-port reference plane: every S is multiplied by exp(+j 2 k0 d cos(theta)) (see `fit_reference_plane_mm`). No TM sign flip is applied here."""
    entries = {}
    with open(path, newline="") as fh:
        for raw in csv.DictReader(fh):
            r = {k.strip(): v for k, v in raw.items()}
            th, st, f = float(r["scan_theta [deg]"]), int(float(r["state []"])), float(r["Freq [GHz]"])
            g = np.exp(1j * 2 * (2 * np.pi * f * 1e9 / C0) * (plane_offset_mm * 1e-3) * math.cos(math.radians(th)))
            m = np.array([[_reim(r, 1, 1), _reim(r, 1, 2)], [_reim(r, 2, 1), _reim(r, 2, 2)]]) * g
            entries[(st, round(f, 6), th, float(phi_deg))] = m
    return entries


def fit_reference_plane_mm(bare_csv, reference_csv):
    """Least-squares reference-plane offset d [mm] with HFSS phase = analytic - 2 k0 d cos(theta) from the bare grounded slab export (state 1) against the analytic
    `bare_grounded_slab_reference.csv`; returns (d_mm, max residual deg)."""
    ref = {}
    with open(reference_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            ref[(round(float(r["f_ghz"]), 6), float(r["theta_deg"]), r["floquet_mode"][0])] = float(r["S11_phase_deg"])
    a, b = [], []
    with open(bare_csv, newline="") as fh:
        for raw in csv.DictReader(fh):
            r = {k.strip(): v for k, v in raw.items()}
            if int(float(r["state []"])) != 1:
                continue
            th, f = float(r["scan_theta [deg]"]), float(r["Freq [GHz]"])
            for mode in (1, 2):
                phase_ref = ref[(round(f, 6), th, str(mode))]
                dphi = np.angle(np.exp(1j * math.radians(phase_ref)) / _reim(r, mode, mode))
                a.append(2 * (2 * np.pi * f * 1e9 / C0) * math.cos(math.radians(th)))
                b.append(dphi)
    a, b = np.array(a), np.array(b)
    d = float(np.sum(a * b) / np.sum(a * a))
    return d * 1e3, float(np.abs(np.degrees(b - a * d)).max())


def run_shift_mm(rows_a, rows_b):
    """Reference-plane shift [mm] of the HFSS data set `rows_b` relative to `rows_a` (both from `read_hfss_reim` with offset 0): median (and 16 / 84 percentiles)
    over the TE / TM co-polar entries with |S| > 0.9 at the common (state, f, theta) of phase(b / a) / (2 k0 cos theta). A pure plane move is a common phase slope,
    so this is independent of the cell physics as long as both runs model the same cell."""
    ds = []
    for k, m1 in rows_a.items():
        if k not in rows_b:
            continue
        st, f, th, _ = k
        for i in (0, 1):
            if abs(m1[i, i]) > 0.9:
                ds.append(np.angle(rows_b[k][i, i] / m1[i, i]) / (2 * (2 * np.pi * f * 1e9 / C0) * math.cos(math.radians(th))) * 1e3)
    return float(np.median(ds)), float(np.percentile(ds, 16)), float(np.percentile(ds, 84))


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


# --- additive HFSS delta table (user-approved calibration of the CP-2 cell, 2026-10-23) -----------------------------------------------------------------------------------------------
class DeltaTable:
    """Additive structure-frame correction `dC[f, theta, state, 2, 2] = C_hfss - C_library` at the HFSS points (phi = 0, f 0.5 GHz grid, theta 0 and 30 deg). `apply` adds it to a
    structure-frame circular response computed for the SAME cell, linearly interpolated in f and theta (theta beyond the HFSS range is held at the last value and flagged). For a MODIFIED
    cell the delta is a TRANSFER assumption (the HFSS-vs-MoM offsets of the old cell are taken to carry over), not a measurement: label `HFSS_CALIBRATED_TRANSFERRED`. phi-independent by assumption
    (HFSS was run at phi = 0 only)."""

    LABEL = "HFSS_CALIBRATED"

    def __init__(self, f_ghz, theta_deg, delta, source):
        self.f, self.theta, self.delta, self.source = (
            np.asarray(f_ghz),
            np.asarray(theta_deg),
            np.asarray(delta),
            str(source),
        )

    @classmethod
    def build(cls, library, hfss, source, *, approved_by_user=False):
        """Delta table from `library` (corrected, structure frame) and HFSS rows `{(state, f, theta, phi): 2x2 TE/TM}` (phi = 0 rows used). Needs `approved_by_user=True` and a provenance string."""
        if not approved_by_user:
            raise CalibrationNotApproved(
                "building an HFSS delta table needs the user's explicit approval (approved_by_user=True)"
            )
        if not source:
            raise ValueError("a calibration needs its data provenance (`source`)")
        rows = {k: v for k, v in hfss.items() if k[3] == 0.0}
        fs = sorted({k[1] for k in rows})
        ths = sorted({k[2] for k in rows})
        delta = np.zeros((len(fs), len(ths), 4, 2, 2), dtype=complex)
        for i, f in enumerate(fs):
            for j, th in enumerate(ths):
                for s in range(4):
                    c_h = to_structure_frame(circular_reflection(rows[(s + 1, f, th, 0.0)]), 0.0)
                    delta[i, j, s] = c_h - library.structure_circular(f, th, 0.0, s + 1)[0]
        return cls(fs, ths, delta, source)

    def at(self, f_ghz, theta_deg, state):
        """Interpolated delta (2x2) for state 1..4; also returns the label."""
        fi = np.interp(f_ghz, self.f, np.arange(len(self.f)))
        ti = np.interp(min(theta_deg, self.theta[-1]), self.theta, np.arange(len(self.theta)))
        i0, j0 = int(np.floor(fi)), int(np.floor(ti))
        i1, j1 = min(i0 + 1, len(self.f) - 1), min(j0 + 1, len(self.theta) - 1)
        wf, wt = fi - i0, ti - j0
        d = self.delta[:, :, state - 1]
        val = (
            (1 - wf) * (1 - wt) * d[i0, j0]
            + wf * (1 - wt) * d[i1, j0]
            + (1 - wf) * wt * d[i0, j1]
            + wf * wt * d[i1, j1]
        )
        label = (
            self.LABEL
            if theta_deg <= self.theta[-1] and self.f[0] <= f_ghz <= self.f[-1]
            else "HFSS_CALIBRATED_EXTRAPOLATED"
        )
        return val, label

    def apply(self, c_structure, f_ghz, theta_deg, state):
        val, label = self.at(f_ghz, theta_deg, state)
        return np.asarray(c_structure) + val, label

    def save(self, path):
        np.savez_compressed(
            path,
            f_ghz=self.f,
            theta_deg=self.theta,
            delta=self.delta,
            source=np.array([self.source]),
            label=np.array([self.LABEL]),
        )

    @classmethod
    def load(cls, path):
        z = np.load(path, allow_pickle=False)
        return cls(z["f_ghz"], z["theta_deg"], z["delta"], str(z["source"][0]))
