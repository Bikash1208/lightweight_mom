"""Phase 20 R5-5: service layer of the "RIS designer" (plane-wave reflecting RIS on the CP-2 family). No heavy compute: one request is a few seconds (calibrated library lookup,
coding synthesis, forward-model pattern, tolerance band). Everything outside the registry domain ends in the named state `NOT_AVAILABLE_FROM_CURRENT_DATA` with the reasons; nothing is
extrapolated. Array-level (supercell) verification is a server job; this layer only reads its results from `results/phase20_r5/supercell_*/` when they exist.
Labels: MODEL_ONLY, NOT_HFSS_VALIDATED (HFSS_CALIBRATED only where the cell responses are inside the HFSS data), LOCAL_RESPONSE_AT_INCIDENT_ANGLE_ASSUMED.
"""

import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from ..hologram import cell_family_registry as registry
from . import coding as cd
from .calibrated_library import (
    LIBRARY_STATUS,
    CalibratedLibrary,
    ToleranceBand,
    read_hfss_reim_tolerance,
    tolerance_plane_offset_mm,
)
from .state_library import OutOfDomain

ROOT = Path(__file__).resolve().parents[2]
FAMILY = "ris_cp2_octagon_madp_rt5880_15"
VARIANTS = ("plain", "dither", "tilt", "refined", "auto")
AUTO_DITHER_FROM_DEG = 35.0  # R5-3: the plain 2-bit gradient's quantisation lobe is worse than dither from about here; array-level dither check pending (R5-4)
HF = ROOT / "hfss_results" / "phase20_r3_cp_cell"
TOL_CSV = HF / "hfss_cp2_tolerance_Coff_0225_025_030.csv"
COFFS = (0.0225, 0.025, 0.030)
MAX_CELLS = 40


def default_request():
    return {
        "frequency_ghz": 10.5, "rows": 20, "cols": 20, "theta_s_deg": 30.0, "phi_s_deg": 0.0, "theta_i_deg": 0.0, "phi_i_deg": 0.0,
        "handedness": "rhcp", "variant": "plain", "seed": 0, "sensitivity": True, "library": "v2",
    }  # fmt: skip


def normalise_request(payload):
    r = {**default_request(), **{k: v for k, v in (payload or {}).items() if v is not None}}
    if r["variant"] not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}")
    if r["library"] not in ("v1", "v2"):
        raise ValueError("library must be 'v1' or 'v2'")
    for k in ("rows", "cols"):
        r[k] = int(r[k])
        if not 2 <= r[k] <= MAX_CELLS:
            raise ValueError(f"{k} must be between 2 and {MAX_CELLS}")
    for k in ("frequency_ghz", "theta_s_deg", "phi_s_deg", "theta_i_deg", "phi_i_deg"):
        r[k] = float(r[k])
    r["seed"] = int(r["seed"])
    r["sensitivity"] = bool(r["sensitivity"])
    return r


def _domain_reasons(r):
    why = list(
        registry.check_request(
            FAMILY, f_ghz=r["frequency_ghz"], theta_deg=r["theta_i_deg"], phi_deg=r["phi_i_deg"]
        ).reasons
    )
    if r["handedness"] != "rhcp":
        why.append(
            "only the right-handed (propagation-relative RHCP in, RHCP out) design is characterised; LHCP is not available"
        )
    if not 0.0 <= r["theta_s_deg"] <= 50.0:
        why.append(
            f"steering angle {r['theta_s_deg']:g} deg is outside the designed 0-50 deg range (60 deg exists in the R5 table only as a band-centre row)"
        )
    if r["phi_s_deg"] not in (0.0, 90.0):
        why.append("only the phi_s = 0 and 90 deg planes are designed in R5")
    if r["theta_s_deg"] > 0 and cd.grating_onset_deg(r["frequency_ghz"]) <= r["theta_s_deg"]:
        why.append(
            f"a grating lobe of the 15 mm lattice is in the visible region at {r['frequency_ghz']:g} GHz for {r['theta_s_deg']:g} deg (onset {cd.grating_onset_deg(r['frequency_ghz']):.1f} deg)"
        )
    if r["theta_i_deg"] > 30.0:
        why.append("plane-wave incidence beyond 30 deg is outside the HFSS-calibrated domain and is not offered")
    return why


@lru_cache(maxsize=4)
def _library(version, with_tolerance=True):
    cal = CalibratedLibrary.load(version=version)
    if with_tolerance and TOL_CSV.exists():
        plane, _, _ = tolerance_plane_offset_mm(
            TOL_CSV,
            HF / "hfss_cp2_states1-4_theta0_30.csv",
            HF / "hfss_bare_slab_theta0_30.csv",
            ROOT / "hfss_exports/phase20_r3_cp_cell/bare_grounded_slab_reference.csv",
        )
        cal.tolerance = ToleranceBand.from_hfss(read_hfss_reim_tolerance(TOL_CSV, plane_offset_mm=plane))
    return cal


def _metrics_dict(m, g_t):
    return {
        "gain_at_target_dbi": round(float(g_t), 2), "peak_gain_dbi": round(float(m.peak_gain_dbi), 2),
        "peak_direction_deg": [round(float(m.peak_theta_deg), 2), round(float(m.peak_phi_deg), 1)],
        "pointing_error_deg": round(float(m.pointing_error_deg), 2), "directivity_dbi": round(float(m.directivity_dbi), 2),
        "sll_db": round(float(m.sll_db), 1), "sll_direction_deg": [round(float(m.sll_theta_deg), 1), round(float(m.sll_phi_deg), 1)],
        "lobes_db_theta_phi": [[round(float(a), 1), round(float(b), 1), round(float(c), 1)] for a, b, c in m.lobes],
        "opposite_handed_peak_db": round(float(m.cross_peak_db), 1), "opposite_handed_at_specular_db": round(float(m.cross_at_specular_db), 1),
        "axial_ratio_db": round(float(m.axial_ratio_db), 2), "radiated_over_incident": round(float(m.radiated_over_incident), 3),
    }  # fmt: skip


def coupling_status():
    """Array-level verification results found in `results/phase20_r5/supercell_*` (merged from the server jobs): which reference tiles pass the 1 dB / 15 deg gate. The designed
    array itself is never called verified by this: it is `COUPLING_NOT_VERIFIED_FOR_THIS_CODING` unless it is one of the verified tiles."""
    out = {
        "state": "COUPLING_NOT_AVAILABLE_FROM_CURRENT_DATA",
        "tasks": [],
        "note": "no supercell results merged yet (jobs r5-supercell-gradients / r5-supercell-random)",
    }
    for sub in ("supercell_gradients", "supercell_random"):
        for f in sorted((ROOT / "results" / "phase20_r5" / sub).glob("*.json")):
            try:
                d = json.loads(f.read_text())
            except ValueError:
                continue
            for res in d.get("results", []):
                out["tasks"].append(
                    {
                        "task": d.get("task", {}).get("name", f.stem),
                        "coding": res.get("label"),
                        "gate_1db_15deg": res.get("gate_1db_15deg"),
                        "max_abs_err_db_above_floor": res.get("max_abs_err_db_above_floor"),
                    }
                )
    if out["tasks"]:
        out["state"] = "COUPLING_RESULTS_FOR_REFERENCE_TILES"
        out["note"] = (
            "reference tiles / windows only; a designed 20 x 20 array is not called verified by them (R5-4 report)"
        )
    return out


def run_design(req, progress=None):
    r = normalise_request(req)
    prog = progress or (lambda f, m: None)
    why = _domain_reasons(r)
    fam = registry.get(FAMILY)
    base = {
        "family": FAMILY,
        "validation_label": fam.validation_label,
        "request": r,
        "known_biases": list(fam.known_biases),
        "labels": list(cd.ASSUMPTIONS),
    }
    if why:
        return {**base, "state": "NOT_AVAILABLE_FROM_CURRENT_DATA", "reasons": why}
    prog(0.1, "loading the calibrated library")
    cal = _library(r["library"], r["sensitivity"])
    variant = r["variant"]
    notes = []
    if variant == "auto":
        variant = "dither" if r["theta_s_deg"] >= AUTO_DITHER_FROM_DEG else "plain"
        notes.append(
            f"auto -> {variant} (dither from {AUTO_DITHER_FROM_DEG:g} deg, R5-3); array-level behaviour of dither is pending the R5-4 tile runs"
        )
    try:
        prog(0.3, "synthesising the coding")
        c = cd.design(
            cal,
            variant,
            r["frequency_ghz"],
            r["theta_s_deg"],
            r["phi_s_deg"],
            r["theta_i_deg"],
            r["phi_i_deg"],
            rows=r["rows"],
            cols=r["cols"],
            seed=r["seed"],
        )
        prog(0.6, "forward model")
        m, g_t = cd.evaluate(cal, c, step=0.005)
        refs = {
            k: _metrics_dict(v, v.peak_gain_dbi)
            for k, v in cd.ideal_references(
                r["theta_s_deg"],
                r["phi_s_deg"],
                r["frequency_ghz"],
                r["theta_i_deg"],
                r["phi_i_deg"],
                r["rows"],
                r["cols"],
                step=0.005,
            ).items()
        }
        sens = None
        if r["sensitivity"] and cal.tolerance is not None and r["theta_i_deg"] <= 30.0:
            prog(0.8, "diode-tolerance band")
            sens = {}
            for coff in COFFS:
                mm, gg = cd.evaluate(cal, c, coff=coff, step=0.005)
                sens[str(coff)] = _metrics_dict(mm, gg)
            sens = {
                "label": "MODEL_ONLY with the user's HFSS C_off differences (theta 0/30, phi 0) added to the calibrated states",
                "by_coff_pf": sens,
                "gain_spread_db": round(
                    max(v["gain_at_target_dbi"] for v in sens.values())
                    - min(v["gain_at_target_dbi"] for v in sens.values()),
                    2,
                ),
            }
        resp = cal.response(r["frequency_ghz"], r["theta_i_deg"], r["phi_i_deg"], 1)
    except OutOfDomain as exc:
        return {**base, "state": "NOT_AVAILABLE_FROM_CURRENT_DATA", "reasons": [str(exc)]}
    prog(1.0, "done")
    counts = {int(s): int(np.sum(c.states == s)) for s in (1, 2, 3, 4)}
    return {
        **base, "state": "OK", "library_label": LIBRARY_STATUS[r["library"]], "reasons": [], "variant_used": variant, "seed": c.seed, "notes": notes + list(resp.notes), "cell_response_labels": list(resp.labels),
        "states": c.states.tolist(), "state_counts": counts, "metrics": _metrics_dict(m, g_t), "ideal_references": refs, "tolerance_sensitivity": sens,
        "grating_onset_deg": round(cd.grating_onset_deg(r["frequency_ghz"]), 1), "coupling": coupling_status(),
        "coupling_for_this_coding": "COUPLING_NOT_VERIFIED_FOR_THIS_CODING: isolated-cell (local-response) model; array-level check = server jobs r5-supercell-gradients / r5-supercell-random",
        "server_job_hint": "To verify this coding at array level, send its state grid (below) to a supercell job: tiles of at most 4 x 4 cells at F = 2n (hours per task); the page does not run it.",
    }  # fmt: skip
