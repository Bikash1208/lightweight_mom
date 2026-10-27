"""Phase 28 (HFSS-first path): coded hologram of the small hexagon from the cell G, and its MODEL_ONLY prediction.

The design is the Phase 26 / Phase 28 item-2 pipeline unchanged (`coded_feasibility`): continuous tensor-CP design on the lattice, slot axis quantised to the four axes of G
(0, 45, 90, 135 deg = states 0 .. 3), pattern by the Phase 26 aperture sum. Two table variants are reported:

* `table_uncapped`: the 10.0 GHz supercell table as it is (leakage up to 4.1 Np/m): what the continuous design asks for. G cannot deliver it.
* `table_capped_at_G`: every leakage constant limited to the mean alpha_beam of the four G states at 10.0 GHz (the fixed anisotropy of G; `capped_table`). THIS is the
  prediction for the board built from G. The conversion efficiency `eta_conv` of the march follows from it.

Held from the earlier phases (stated assumptions, not results): launch fraction 0.592 (Phase 26 bare slab, 10.54 GHz); cells respond locally; the four G states have the
same alpha (they differ by up to 25 %: not modelled); state axes are the ideal 0/45/90/135 (the as-built axes differ by < 0.3 deg).
"""

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.hologram import coded_feasibility as cf  # noqa: E402
from src.hologram.coded_cell_candidates import RationalSectorCell  # noqa: E402
from src.hologram.coded_cell_eval import principal_metrics  # noqa: E402,F401
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402

FREQ_GHZ = 10.0
PITCH_MM = 6.4
ACROSS_FLATS_MM = 70.0
BEAMS = [(0.0, 0.0), (30.0, 225.0), (45.0, 225.0)]
TABLE_GLOB = str(ROOT / "results/phase28/f10/supercell_table_*.json")
G_PARAMS = {
    "hub_half": 0.0012,
    "neck_gap": 0.0004,
    "neck_width": 0.0004,
    "gap_t": 0.0002,
    "arm_outer": 0.0030,
    "diag_outer": 0.0028,
    "clip": 0.0030,
    "mode": "bowtie",
}
cf.set_frequency(FREQ_GHZ * 1e9)


def g_cell():
    return RationalSectorCell(period=PITCH_MM * 1e-3, **G_PARAMS)


def circumradius_mm(across_flats_mm=ACROSS_FLATS_MM):
    return across_flats_mm / 2.0 / math.cos(math.radians(30.0))


def layout(across_flats_mm=ACROSS_FLATS_MM):
    """Lattice cells of the board: `(x_mm, y_mm)` arrays (node at the feed, outside the 6.5 mm feed keep-out, inside the hexagon with a vertex at +-y)."""
    x, y = cf.lattice(cf.hexagon_m(circumradius_mm(across_flats_mm) * 1e-3), PITCH_MM * 1e-3)
    order = np.lexsort((x, y))
    return x[order] * 1e3, y[order] * 1e3


def great_circle_cut(t_pk, p_pk, azimuth_deg, half_deg=70.0, step=0.5):
    """Points `(theta, phi, angle)` on the great circle through the peak, `azimuth_deg` = direction of the cut in the local tangent plane."""
    t, p = np.radians(t_pk), np.radians(p_pk)
    n = np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
    e_t = np.array([np.cos(t) * np.cos(p), np.cos(t) * np.sin(p), -np.sin(t)])
    e_p = np.array([-np.sin(p), np.cos(p), 0.0])
    a = np.radians(azimuth_deg)
    e = np.cos(a) * e_t + np.sin(a) * e_p
    ang = np.arange(-half_deg, half_deg + 1e-9, step)
    v = np.cos(np.radians(ang))[:, None] * n + np.sin(np.radians(ang))[:, None] * e
    ok = v[:, 2] >= 0.0
    th = np.degrees(np.arccos(np.clip(v[:, 2], -1, 1)))
    ph = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
    return th[ok], ph[ok], ang[ok]


def hpbw_deg(jx, jy, x, y, green, pitch, t_pk, p_pk):
    """Half-power width of the RHCP main beam on two orthogonal great-circle cuts through the peak: the contiguous region around the peak with power >= half the peak."""
    out = []
    for az in (0.0, 90.0):
        th, ph, ang = great_circle_cut(t_pk, p_pk, az)
        et, ep = cf.far_field(jx, jy, x, y, green, th, ph, pitch)
        pw = cf._main_power(et, ep, "cp")
        i0 = int(np.argmin(np.abs(ang)))
        pk = pw[i0]
        lo = i0
        while lo > 0 and pw[lo - 1] >= 0.5 * pk:
            lo -= 1
        hi = i0
        while hi < len(pw) - 1 and pw[hi + 1] >= 0.5 * pk:
            hi += 1
        out.append(float(ang[hi] - ang[lo]) if (lo > 0 and hi < len(pw) - 1) else None)
    return out


def sll_outside_db(jx, jy, x, y, green, pitch, t_pk, p_pk, hpbw):
    """Strongest RHCP level [dB re peak] farther than TWICE the mean half-power width from the peak (`analyse` excludes only 15 deg, which for a 23 deg beam is inside the main lobe skirt)."""
    if any(h is None for h in hpbw):
        return None
    th, ph = np.meshgrid(np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 2.0), indexing="ij")
    et, ep = cf.far_field(jx, jy, x, y, green, th, ph, pitch)
    pw = cf._main_power(et, ep, "cp")
    pk, _, _, _, _ = cf._refined_peak(jx, jy, x, y, green, pitch, "cp", t_pk, p_pk)
    far = cf._angle_deg(th, ph, t_pk, p_pk) > 2.0 * float(np.mean(hpbw))
    return float(10 * math.log10(np.max(np.where(far, pw, 0.0)) / pk)) if far.any() else None


def predict(beam, alpha_g_np_per_m, across_flats_mm=ACROSS_FLATS_MM):
    table = SupercellTable.from_files(TABLE_GLOB)
    poly = cf.hexagon_m(circumradius_mm(across_flats_mm) * 1e-3)
    pitch = PITCH_MM * 1e-3
    green = cf.green_function()
    res = {"beam_theta_phi_deg": list(beam)}
    for name, tab in (("table_uncapped", table), ("table_capped_at_G", cf.capped_table(table, alpha_g_np_per_m))):
        x, y, d, march = cf.design_cells(tab, poly, pitch, beam[0], beam[1])
        eta = float(march["eta_conv"])
        entry = {"n_cells": int(x.size), "eta_conv": eta}
        for label, bits in (("continuous", None), ("quantised_4_axes", 2)):
            idx, axis = cf.cp_states(d, bits)
            jx, jy = cf.cp_currents(d, axis)
            r = cf.analyse(jx, jy, x, y, green, pitch, "cp", beam)
            r["gain_dbi_launch_bare_0p592"] = cf.launch_gain_dbi(r["d_ap_dbi"], eta)
            r["hpbw_deg_two_cuts"] = hpbw_deg(jx, jy, x, y, green, pitch, r["peak_theta_deg"], r["peak_phi_deg"])
            r["sll_db_outside_2hpbw"] = sll_outside_db(jx, jy, x, y, green, pitch, r["peak_theta_deg"], r["peak_phi_deg"], r["hpbw_deg_two_cuts"])
            if bits is not None:
                r["state_counts"] = [int(np.sum(idx == k)) for k in range(4)]
                entry["state_index"] = [int(v) for v in idx]
                entry["x_mm"] = [float(v) for v in x * 1e3]
                entry["y_mm"] = [float(v) for v in y * 1e3]
            entry[label] = r
        entry["fixed_source_quantisation_loss_db"] = cf.fixed_source_loss_db(entry["continuous"], entry["quantised_4_axes"])
        res[name] = entry
    return res
