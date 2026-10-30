"""Phase 29 option A: segmented aperture for +-45 deg coverage (MODEL_ONLY). 150 mm hexagon, 10.54 GHz, 7 regions of equal area (central hexagon of circumradius R/sqrt(7) + 6 outer 60-deg sectors), one feed per region
(centre: broadside; sectors: feed at radius 90 mm on the sector axis, beam theta 30 deg, phi 0/60/..., 300). Each region's cells take the slot/orientation of the EXISTING single-feed tensor synthesis (`design_tensor`, one feed,
the region's own target) with the target strength tapered to zero over `BAND` mm inside the region boundary (cosine; `cell_weight`). The board = union of the regions; every feed's current is then re-evaluated on the combined board
(its surface wave crosses the other regions' cells). Outputs hfss_exports/phase29_segmented/ (STEP, AEDT script, port table, steering LUT) and results/phase29/segmented_summary.json.
    python scripts/research/phase29_segmented.py"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/research"))
RES = ROOT / "results/phase29"
OUT = ROOT / "hfss_exports/phase29_segmented"
SYNC = Path(r"C:\hmsa_sync\data\phase29_segmented")
R_MM, BAND, FEED_R, TILT = 150.0, 12.8, 90.0, 30.0
RC = R_MM / np.sqrt(7.0)
NORMALS = np.array([[np.cos(np.radians(60 * k)), np.sin(np.radians(60 * k))] for k in range(6)])
A_IN = RC * np.cos(np.radians(30))


def region_of(x, y):
    """Region index 0 (central hexagon) or 1..6 (sector with axis 60 (j-1) deg), and the distance (mm) to the region's inner boundaries."""
    ap = (NORMALS @ np.vstack([x, y])).max(axis=0)  # hexagon 'apothem' coordinate along the edge normals
    inner = ap <= A_IN
    phi = np.degrees(np.arctan2(y, x)) % 360.0
    j = (np.rint(phi / 60.0).astype(int) % 6) + 1
    dphi = np.abs(((phi - 60.0 * (j - 1) + 180.0) % 360.0) - 180.0)
    d_rad = np.hypot(x, y) * np.sin(np.radians(np.maximum(30.0 - dphi, 0.0)))
    reg = np.where(inner, 0, j)
    dist = np.where(inner, A_IN - ap, np.minimum(d_rad, ap - A_IN))
    return reg, dist


def taper(dist):
    return np.where(dist >= BAND, 1.0, 0.5 * (1 - np.cos(np.pi * np.clip(dist, 0, BAND) / BAND)))


def main():
    import phase29_calibration as cal
    import phase29_tensor_ladder as L
    import phase29_twofeed_step as ts
    from phase29_pointB_sweep import port_metrics

    from src.export import shared7_export as ex
    from src.hologram import shared_aperture_multiport as sa
    from src.hologram import shared_aperture_tensor as st
    from src.hologram import tensor_cp_design as t

    L.init()
    table, lut = L._S["table"], L._S["lut"]
    ports = [(0.0, 0.0, 0.0, 0.0)] + [
        (FEED_R * np.cos(np.radians(60 * k)), FEED_R * np.sin(np.radians(60 * k)), TILT, 60.0 * k) for k in range(6)
    ]
    cx, cy = st.aperture_cells(R_MM, ports, (0.0, 0.0))  # board cells (mm): every feed's keep-out removed
    reg, _ = region_of(cx, cy)
    slot, ia = np.zeros(cx.size, int), np.zeros(cx.size, int)
    key = {(round(x, 2), round(y, 2)): i for i, (x, y) in enumerate(zip(cx, cy, strict=True))}
    for j, p in enumerate(ports):
        px, py = st.aperture_cells(R_MM, [p], (0.0, 0.0))
        r_, d_ = region_of(px, py)
        w = np.where(r_ == j, taper(d_), 0.0)
        dj = st.design_tensor(table, lut, R_MM, [p], [0.0], n_psi=360, cell_weight=w)
        ia_j = np.rint(dj["axis_rad"] / lut.step).astype(int)
        for x, y, s_, a_, rr in zip(px, py, dj["slot"], ia_j, r_, strict=True):
            if rr == j:
                i = key.get((round(x, 2), round(y, 2)))
                if i is not None:
                    slot[i], ia[i] = s_, a_
        print("region", j, "cells", int((reg == j).sum()), flush=True)
    xm, ym = cx * 1e-3, cy * 1e-3
    poly = sa.hexagon_m(R_MM)
    feeds = [(p[0] * 1e-3, p[1] * 1e-3) for p in ports]
    beams = [(np.radians(p[2]), np.radians(p[3])) for p in ports]
    cxy = np.c_[xm, ym]
    ana = [
        st.march_analysis_tensor(table, lut, poly, feeds[j], t.K0, beams[j][0], beams[j][1], cxy, slot, ia, 240)
        for j in range(7)
    ]
    cur = [st.cell_current(ana[j], lut, table, feeds[j], xm, ym, slot, ia, beams[j][0], beams[j][1]) for j in range(7)]
    d = {"x_m": xm, "y_m": ym, "J": cur, "eta_conv": [a["eta_conv"] for a in ana], "beams": beams, "feeds_m": feeds}
    th, ph = np.arange(0.0, 90.0 + 1e-9, 1.0), np.arange(0.0, 360.0, 1.0)
    F = sa.port_fields(d, L._S["green"], th, ph)
    pm = port_metrics(F, th, ph, beams)
    beam_deg = [(float(np.degrees(b[0])), float(np.degrees(b[1]))) for b in beams]
    rows = []
    for j, p in enumerate(pm):
        gr, sl = cal.calibrate(p["peak_gain_dbi"], p["sll_db"])
        rows.append(
            {
                "port": j + 1,
                "region": "centre" if j == 0 else f"sector {j} (axis {60 * (j - 1)} deg)",
                "feed_x_mm": round(ports[j][0], 3),
                "feed_y_mm": round(ports[j][1], 3),
                "design_beam_theta_deg": beam_deg[j][0],
                "design_beam_phi_deg": beam_deg[j][1],
                "pred_peak_theta_phi_deg": [round(x, 1) for x in p["peak_theta_phi"]],
                "hpbw_deg": round(p["hpbw_deg"], 1),
                "eta_conv": round(d["eta_conv"][j], 3),
                "gain_raw_dbi": round(p["peak_gain_dbi"], 1),
                "sll_raw_db": round(p["sll_db"], 1),
                "gain_calibrated_dbi": round(gr, 1),
                "sll_calibrated_db": round(sl, 1),
            }
        )
    # coverage: switched (best single port) and blended (MRC over all 7, ports taken as uncoupled)
    dirs = sa.scan_directions(45.0, 3.0, 15.0)
    G = np.array([sa.gain_db(*f) for f in F])  # (7, th, ph)
    tab = sa.scan_table(F, th, ph, dirs, n_active=7)
    sw, bl, best = [], [], []
    for r_ in tab:
        it, ip = int(round(r_["theta_deg"])), int(round(r_["phi_deg"])) % 360
        gj = G[:, it, ip]
        sw.append(float(gj.max()))
        best.append(int(np.argmax(gj)) + 1)
        bl.append(float(10 * np.log10(np.sum(10 ** (gj / 10)))))
    sw, bl = np.array(sw), np.array(bl)
    th_d = np.array([r_["theta_deg"] for r_ in tab])
    bands = {}
    for lo, hi in ((0, 15), (15, 30), (30, 45)):
        m = (th_d >= lo) & (th_d <= hi)
        bands[f"{lo}-{hi}"] = {"switched_min": round(float(sw[m].min()), 1), "blended_min": round(float(bl[m].min()), 1)}
    cov = {
        "switched_min_dbi_theta_le_45": round(float(sw.min()), 1),
        "blended_min_dbi_theta_le_45": round(float(bl.min()), 1),
        "switched_mean_dbi": round(float(sw.mean()), 1),
        "blended_mean_dbi": round(float(bl.mean()), 1),
        "calibrated_min_dbi_switched_blended": [
            round(float(sw.min()) + cal.GAIN_OFFSET_DB, 1),
            round(float(bl.min()) + cal.GAIN_OFFSET_DB, 1),
        ],
        "worst_switched_direction": [float(th_d[int(np.argmin(sw))]), float(tab[int(np.argmin(sw))]["phi_deg"])],
        "worst_blended_direction": [float(th_d[int(np.argmin(bl))]), float(tab[int(np.argmin(bl))]["phi_deg"])],
        "by_theta_band": bands,
        "blended_mrc_sll_median_db": round(float(np.nanmedian([r_["sll_db"] for r_ in tab])), 1),
        "blended_mrc_pointing_err_max_deg": round(float(max(r_["pointing_error_deg"] for r_ in tab)), 1),
    }
    S = {
        "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED", "INTER_PORT_COUPLING_NOT_MODELLED", cal.LABEL],
        "n_cells": int(cx.size),
        "region_cells": [int((reg == j).sum()) for j in range(7)],
        "ports": rows,
        "coverage": cov,
    }
    (RES / "segmented_summary.json").write_text(json.dumps(S, indent=1, default=float))
    # ---- package
    OUT.mkdir(parents=True, exist_ok=True)
    cells = [
        (float(x), float(y), float(table.slots[s_][0]), float(table.slots[s_][1]), float(np.degrees(a_ * lut.step)))
        for x, y, s_, a_ in zip(cx, cy, slot, ia, strict=True)
    ]
    fd = [(p[0], p[1]) for p in ports]
    cd = [dict(zip(("x_mm", "y_mm", "slot_L_mm", "slot_W_mm", "slot_angle_deg"), r_, strict=True)) for r_ in cells]
    (OUT / "build_board.py").write_text(ex.aedt_build_script(cd, fd, name="phase29_segmented"), encoding="utf-8")
    with open(OUT / "ports_and_excitations.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(list(rows[0]) + [f"excite_P{k}" for k in range(1, 8)])
        for r_ in rows:
            w_.writerow(list(r_.values()) + [("1.0 W / 0 deg" if k == r_["port"] else "0 (matched)") for k in range(1, 8)])
    with open(OUT / "steering_lut.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(
            ["theta_deg", "phi_deg", "best_single_port", "switched_gain_raw_dbi", "blended_gain_raw_dbi", "blended_gain_calibrated_dbi"]
            + [f"P{k}_power_W" for k in range(1, 8)]
            + [f"P{k}_phase_deg" for k in range(1, 8)]
        )
        for r_, s1, b1, bp in zip(tab, sw, bl, best, strict=True):
            wv = np.array([complex(a, b) for a, b in r_["weights"]])
            w_.writerow(
                [r_["theta_deg"], r_["phi_deg"], bp, f"{s1:.1f}", f"{b1:.1f}", f"{b1 + cal.GAIN_OFFSET_DB:.1f}"]
                + [f"{abs(v) ** 2:.4f}" for v in wv]
                + [f"{np.degrees(np.angle(v)):.1f}" for v in wv]
            )
    stepf = OUT / "segmented_board.step"
    ts.build(cells, fd, 270.0, 310.0, stepf)
    names, ok = ts.check(stepf, len(cells), 7)
    SYNC.mkdir(parents=True, exist_ok=True)
    for f_ in (stepf, OUT / "ports_and_excitations.csv", OUT / "steering_lut.csv"):
        shutil.copy2(f_, SYNC / f_.name)
    (OUT / "README.md").write_text(
        "# phase29_segmented (option A, MODEL_ONLY)\n\n150 mm hexagon, 10.54 GHz, same cell/stack-up/settings as `phase29_twofeed` and `phase29_sevenport`: 7 equal-area regions (central hexagon + 6 sectors), "
        "one port per region (P1 centre, broadside; P2..P7 sector axes 0, 60, ... 300 deg, beam theta 30 deg). Cells of a region come from the single-feed synthesis; target strength tapered (cosine) over 12.8 mm inside each region boundary.\n\n"
        "- `segmented_board.step` (gitignored; copy in `C:/hmsa_sync/data/phase29_segmented/`), `build_board.py` (AEDT), `ports_and_excitations.csv` (port table; predictions raw and calibrated: gain -2.5 dB, SLL +3.5 dB from the two-feed HFSS).\n"
        "- `steering_lut.csv`: per direction (theta 0..45 step 3, phi step 15): best single port, switched and blended gain, and the 7-port MRC weights (power W and phase deg, sum of power = 1 W) for Results > Edit Sources. Ports taken as uncoupled.\n"
        "- Manual STEP setup: STEP section of `../phase29_twofeed/README.md` (ports 1..7).\n",
        encoding="utf-8",
    )
    print("check_ok", ok, "cells", len(cells))
    print(json.dumps(rows))
    print(json.dumps(cov))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
