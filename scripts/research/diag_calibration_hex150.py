"""FOLLOW-UP DIAGNOSTIC item 1: the same variant test as diag_hfss_vs_model.py on the HFSS-validated Phase 26 static tensor CP board (hex150_tensor_cp_single, 10.54 / 10.5 GHz, RHCP beam at (30, 225)).
    python diag_calibration_hex150.py <out.json>
Model: the Phase 26 engine as used for the HFSS export (SupercellTable, march_rays, cell_design hand=+1, LayeredGreenFunction RT5880 1.575 mm, `_vector_far_field` = the same far-field routine and exp(+j k.r) sign that
sub6_board_model4 uses).  HFSS data: realized-gain RHCP / LHCP in dB (power only, theta 0..180, phi 0..359, 1 deg), so the test is on power patterns (RHCP and LHCP separately, linear power correlation over theta <= 60).
Variants: dihedral transform of the board x conj of the radiated field x conj of the cell excitation phase (the model fields are transformed analytically: mirror flips the handedness, conj of the field swaps RHCP / LHCP).
"""
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np

from src.green_function import LayeredGreenFunction
from src.hologram import tensor_cp_design as TD
from src.hologram.tensor_leaky_cp import SupercellTable, cell_design, march_rays
from src.pattern.vector_aperture_pattern import _vector_far_field

SQ2 = math.sqrt(2.0)


def load_hfss(path):
    th, ph, lh, rh = [], [], [], []
    for r in list(csv.reader(open(path)))[1:]:
        ph.append(float(r[1])); th.append(float(r[2])); lh.append(float(r[3])); rh.append(float(r[4]))
    th, ph = np.array(th), np.array(ph)
    T, P = np.unique(th), np.unique(ph)
    R = np.full((T.size, P.size), np.nan)
    L = np.full_like(R, np.nan)
    ti, pi = np.searchsorted(T, th), np.searchsorted(P, ph)
    R[ti, pi] = 10 ** (np.array(rh) / 10)
    L[ti, pi] = 10 ** (np.array(lh) / 10)
    return T, P[:360], R[:, :360], L[:, :360]    # phi 0 ... 360 is exported (361 columns): drop the duplicate 360


def pidx(phi_deg, P):
    return np.rint(np.mod(phi_deg, 360.0)).astype(int) % P.size


def main():
    outp = sys.argv[1]
    T, P, HR, HL = load_hfss(ROOT / "hfss_results/hex150_tensor_cp_single/realized_gain_RHCP_LHCP_10p5GHz_1deg.csv")
    assert P.size == 360 and T[0] == 0
    table = SupercellTable.from_files(TD.TABLE_GLOB)
    poly = TD.hexagon()
    x, y = TD.lattice_cells(poly)
    th0, ph0 = np.radians(30.0), np.radians(225.0)
    m = march_rays(table, poly, (0.0, 0.0), TD.K0, th0, ph0, 0.9, taper_c=0.5, n_psi=720)
    d = cell_design(m, table, x, y, (0.0, 0.0), TD.K0, th0, ph0, hand=+1)
    # sanity: the design reproduces the cell table that was sent to HFSS
    lay = list(csv.DictReader(open(ROOT / "hfss_exports/hex150_tensor_cp_single/layout.csv")))
    lx = np.array([float(r["x_mm"]) for r in lay]); ly = np.array([float(r["y_mm"]) for r in lay]); la = np.array([float(r["slot_angle_deg"]) for r in lay])
    k = [int(np.argmin(np.hypot(lx - a * 1e3, ly - b * 1e3))) for a, b in zip(x, y)]
    dang = ((np.degrees(d["slot_angle_rad"]) - la[k] + 90) % 180) - 90
    sanity = dict(n_cells_design=int(x.size), n_cells_layout=int(lx.size), max_pos_err_mm=float(np.max(np.hypot(lx[k] - x * 1e3, ly[k] - y * 1e3))), max_axis_err_deg=float(np.max(np.abs(dang))))
    print("layout.csv sanity", sanity, flush=True)
    green = LayeredGreenFunction(TD.SUB, TD.F, backing="pec")
    th_g = np.arange(0, 91, 1.0)
    ph_g = np.arange(0, 360, 1.0)
    TG, PG = np.meshgrid(th_g, ph_g, indexing="ij")

    def field(jx, jy):
        et = np.zeros(TG.size, complex); ep = np.zeros(TG.size, complex)
        a, b = TG.ravel(), PG.ravel()
        for lo in range(0, a.size, 4000):
            sl = slice(lo, lo + 4000)
            et[sl], ep[sl] = _vector_far_field(a[sl], b[sl], x, y, jx, jy, TD.F, green, TD.PITCH, TD.PITCH)
        return et.reshape(TG.shape), ep.reshape(TG.shape)

    F0 = field(d["jx"], d["jy"])
    F1 = field(np.conj(d["jx"]), np.conj(d["jy"]))
    nt = 91
    # HFSS power on the model grid (theta 0..90)
    HRm, HLm = HR[:nt], HL[:nt]

    def trans(e, kind, arg):
        et, ep = e
        if kind is None:
            return et, ep
        if kind == "rot":
            s = pidx(np.arange(360) - arg, P)
            return et[:, s], ep[:, s]
        s = pidx(2 * arg - np.arange(360), P)
        return et[:, s], -ep[:, s]

    def corr(a, b, tmax=60):
        x_, y_ = a[:tmax + 1].ravel(), b[:tmax + 1].ravel()
        x_, y_ = x_ - x_.mean(), y_ - y_.mean()
        return float(np.dot(x_, y_) / math.sqrt(np.dot(x_, x_) * np.dot(y_, y_) + 1e-300))

    def peak(p):
        i = np.unravel_index(np.argmax(p[:61]), p[:61].shape)
        return (int(th_g[i[0]]), int(ph_g[i[1]]))

    rows = []
    geoms = [("id", None, None), ("rot90", "rot", 90), ("rot180", "rot", 180), ("rot270", "rot", 270), ("mirror x->-x", "mir", 90), ("mirror y->-y", "mir", 0), ("swap x<->y", "mir", 45), ("swap x<->-y", "mir", 135)]
    for gname, kind, arg in geoms:
        for conjE in (0, 1):
            for conjJ in (0, 1):
                et, ep = trans(F1 if conjJ else F0, kind, arg)
                if conjE:
                    et, ep = np.conj(et), np.conj(ep)
                rh, lh = np.abs((et + 1j * ep) / SQ2) ** 2, np.abs((et - 1j * ep) / SQ2) ** 2
                rows.append(dict(variant=gname + (" + conj(E)" if conjE else "") + (" + conj(cell phase)" if conjJ else ""), corr_RHCP=corr(rh, HRm), corr_LHCP=corr(lh, HLm), peak_RHCP=peak(rh), peak_LHCP=peak(lh)))
    rows.sort(key=lambda r: -(r["corr_RHCP"] + r["corr_LHCP"]))
    hp = peak(HRm)
    tot = HRm + HLm
    res = dict(label="MODEL_ONLY calibration on the HFSS-validated hex150 tensor CP board (10.5 GHz)", layout_sanity=sanity, hfss_peak_RHCP_theta_phi=hp, hfss_peak_RHCP_dBi=float(10 * np.log10(HR[:91].max())),
               model_peak_RHCP_theta_phi=peak(np.abs((F0[0] + 1j * F0[1]) / SQ2) ** 2), variants_sorted=rows)
    Path(outp).write_text(json.dumps(res, indent=1, default=str), encoding="utf-8", newline="\n")
    print("HFSS RHCP peak", hp, round(res["hfss_peak_RHCP_dBi"], 2), "dBi; model peak", res["model_peak_RHCP_theta_phi"])
    for r in rows[:8]:
        print("%-44s corr RHCP %.3f LHCP %.3f  peak RHCP %s LHCP %s" % (r["variant"], r["corr_RHCP"], r["corr_LHCP"], r["peak_RHCP"], r["peak_LHCP"]))
    print("--- requested")
    for r in rows:
        if r["variant"] in ("id", "rot180", "mirror x->-x", "mirror y->-y", "id + conj(E)", "id + conj(cell phase)", "swap x<->y"):
            print("%-44s corr RHCP %.3f LHCP %.3f  peak RHCP %s LHCP %s" % (r["variant"], r["corr_RHCP"], r["corr_LHCP"], r["peak_RHCP"], r["peak_LHCP"]))


if __name__ == "__main__":
    main()
