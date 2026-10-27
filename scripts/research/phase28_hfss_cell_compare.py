"""Phase 28: HFSS check of the loop cell R4 (user's run 2026-10-27, hfss_results/phase28_loop_R4/ON.csv, OFF.csv): sheet admittances per Floquet mode from the complex S, compared with the MoM prediction.

    python scripts/research/phase28_hfss_cell_compare.py <out.json>

Method (the advisor's, reproduced): per mode Y_in = Y0 (1 - S)/(1 + S) (normal incidence, port de-embedded to the substrate surface), Y_s = Y_in - Y_slab (grounded RT/duroid 5880 3.175 mm, tan d 0.0009). The two modes are diagonal
(|S12| ~ 1e-5). The MoM states compared are the HFSS-representable ones (R only 2.4 ohm, C only 28.2 fF; tags in measured_rows.jsonl). Labels: HFSS_CHECKED (one run, one board, ROT 0), MODEL_ONLY.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts" / "research"))
from src.dispersion.sheet_model import _y_slab  # noqa: E402
from src.frequency import angular_frequency, free_space_wavenumber  # noqa: E402
from src.materials import Substrate  # noqa: E402

Y0 = 1.0 / 376.730313668
SUB = Substrate("rt5880_3p175", 2.2, 0.0009, 3.175e-3)
HF = ROOT / "hfss_results/phase28_loop_R4"


def read_hfss(name):
    rows = list(csv.reader((HF / name).open(encoding="utf-8")))[1:]
    a = np.array([[float(v) for v in r] for r in rows])
    f = a[:, 0]
    s = np.empty((len(f), 2, 2), complex)
    s[:, 0, 0] = a[:, 1] + 1j * a[:, 2]; s[:, 1, 0] = a[:, 3] + 1j * a[:, 4]; s[:, 0, 1] = a[:, 5] + 1j * a[:, 6]; s[:, 1, 1] = a[:, 7] + 1j * a[:, 8]
    return f, s


def sheet_y(f_ghz, s):
    out = []
    for f, m in zip(f_ghz, s, strict=True):
        yin = np.array([Y0 * (1 - m[i, i]) / (1 + m[i, i]) for i in range(2)])
        ys = _y_slab(0.0, angular_frequency(f * 1e9), SUB.eps_complex, free_space_wavenumber(f * 1e9), SUB.thickness)
        out.append(yin - ys)
    return np.array(out)  # (n, 2 modes) complex


def mom_states():
    """{(label, f_ghz, mesh): (Yxx, Yyy)} for the MoM R-only / C-only states."""
    import phase28_loop_measured as lm

    rows = []
    for n in ("measured_rows.jsonl", "measured_extra_rows.jsonl"):
        rows += lm.rows_of(ROOT / "results/phase28/f10/loop_R4" / n)
    out = {}
    for r in rows:
        if r["scale"] is None and r["state"] in ("ON_10mA_Ronly", "OFF_0V_Conly", "ON_10mA", "OFF_0V"):
            y = np.array([[complex(*c) for c in row] for row in r["y_s"]])
            out[(r["state"], round(r["f_hz"] / 1e9, 4), r["mesh"])] = (y[0, 0], y[1, 1], y[0, 1])
    return out


def main(dst):
    fo, so = read_hfss("ON.csv"); ff, sf = read_hfss("OFF.csv")
    assert np.allclose(fo, ff)
    yh = {"ON.csv": sheet_y(fo, so), "OFF.csv": sheet_y(ff, sf)}
    mom = mom_states()
    res = {"labels": ["HFSS_CHECKED", "MODEL_ONLY"], "hfss_files": ["ON.csv", "OFF.csv"], "frequencies_ghz": [float(v) for v in fo], "rows": []}
    mm = lambda lab, f, mesh=60: mom[(lab, f, mesh)]  # noqa: E731
    for f in (9.5, 9.75, 10.0, 10.25, 10.5):
        i = int(np.argmin(np.abs(fo - f)))
        row = {"f_ghz": f}
        for name, y in yh.items():
            row[name] = {"B_mode1_ms": y[i, 0].imag * 1e3, "B_mode2_ms": y[i, 1].imag * 1e3, "G_mode1_ms": y[i, 0].real * 1e3, "G_mode2_ms": y[i, 1].real * 1e3}
        for lab in ("ON_10mA_Ronly", "OFF_0V_Conly"):
            xx, yy, xy = mm(lab, f)
            row["MoM_" + lab] = {"B_xx_ms": xx.imag * 1e3, "B_yy_ms": yy.imag * 1e3, "G_xx_ms": xx.real * 1e3, "G_yy_ms": yy.real * 1e3}
        res["rows"].append(row)
    # band ΔB along the diode axis (mode 2 = x) in HFSS, either pairing
    res["dB_band_ms"] = {"hfss_|B2(OFF.csv) - B2(ON.csv)|": [float(abs(yh["OFF.csv"][i, 1].imag - yh["ON.csv"][i, 1].imag) * 1e3) for i in range(len(fo))],
                         "mom_R_only_minus_C_only_xx": {f"{f:g}": float(abs(mm("ON_10mA_Ronly", f)[0].imag - mm("OFF_0V_Conly", f)[0].imag) * 1e3) for f in (9.5, 9.75, 10.0, 10.25, 10.5)}}
    Path(dst).write_text(json.dumps(res, indent=1), encoding="utf-8")
    r = next(x for x in res["rows"] if x["f_ghz"] == 10.0)
    print(json.dumps(r, indent=1))
    print("dB band hfss", [round(v, 2) for v in res["dB_band_ms"]["hfss_|B2(OFF.csv) - B2(ON.csv)|"][::5]])


if __name__ == "__main__":
    main(sys.argv[1])


def calibration(dst_cal, dst_cmp):
    """HFSS-checked sheet admittances of the two diode states of the loop cell R4 (mode 2 = x = along the diode, mode 1 = y), the Delta to the MoM and the leakage with each. The state NAMES follow the physics, not the file
    names: the lossy state (R = 2.4 ohm) is the one in OFF.csv and the open one (C = 28.2 fF) is in ON.csv (labels swapped in the HFSS project)."""
    import phase28_scalar_ring_cell as rc

    fo, so = read_hfss("ON.csv"); ff, sf = read_hfss("OFF.csv")
    y_c = sheet_y(fo, so)   # ON.csv  -> the C (open diode) state
    y_r = sheet_y(ff, sf)   # OFF.csv -> the R (closed diode) state
    mom = mom_states()
    cal = {"labels": ["HFSS_CHECKED", "MODEL_ONLY"], "scope": "loop cell R4, RT/duroid 5880 3.175 mm, 9 mm period, ROT 0, normal incidence, diode as R only 2.4 ohm (closed) / C only 28.2 fF (open); one HFSS run (2026-10-27); valid 9 - 11 GHz",
           "mode_to_axis": {"mode1": "y (across the diode)", "mode2": "x (along the diode)"}, "state_names": {"R": "OFF.csv (labelled OFF in the HFSS project, carries the 2.4 ohm boundary: lossy, isotropic)", "C": "ON.csv (labelled ON in the HFSS project, carries 28.2 fF: open, anisotropic)"},
           "f_ghz": [float(v) for v in fo],
           "Y_R": {"yy_re_im_S": [[float(v.real), float(v.imag)] for v in y_r[:, 0]], "xx_re_im_S": [[float(v.real), float(v.imag)] for v in y_r[:, 1]]},
           "Y_C": {"yy_re_im_S": [[float(v.real), float(v.imag)] for v in y_c[:, 0]], "xx_re_im_S": [[float(v.real), float(v.imag)] for v in y_c[:, 1]]},
           "vs_mom": {}}
    for f in (9.5, 9.75, 10.0, 10.25, 10.5):
        i = int(np.argmin(np.abs(fo - f)))
        r_m, c_m = mom[("ON_10mA_Ronly", f, 60)], mom[("OFF_0V_Conly", f, 60)]
        def pm(r_xx, r_yy, c_xx, c_yy):
            mk = lambda xx, yy: {"substrate": "rt5880_3p175", "f_hz": f * 1e9, "y_s": [[[xx.real, xx.imag], [0, 0]], [[0, 0], [yy.real, yy.imag]]]}  # noqa: E731
            return rc.pair_metrics(mk(r_xx, r_yy), mk(c_xx, c_yy))
        h = pm(y_r[i, 1], y_r[i, 0], y_c[i, 1], y_c[i, 0])
        m = pm(r_m[0], r_m[1], c_m[0], c_m[1])
        cal["vs_mom"][f"{f:g}"] = {"hfss": {k: float(h[k]) for k in ("d_b_ms", "b_bar_ms", "alpha_exact", "alpha_first_order", "g_on_ms") if k in h},
                                    "mom": {k: float(m[k]) for k in ("d_b_ms", "b_bar_ms", "alpha_exact", "alpha_first_order", "g_on_ms") if k in m},
                                    "dB_ratio_hfss_over_mom": float(h["d_b_ms"] / m["d_b_ms"]), "b_bar_ratio": float(h["b_bar_ms"] / m["b_bar_ms"]), "alpha_ratio": float(h["alpha_exact"] / m["alpha_exact"])}
    cal["delta_to_mom"] = {}
    for f in (9.5, 9.75, 10.0, 10.25, 10.5):
        i = int(np.argmin(np.abs(fo - f)))
        r_m, c_m = mom[("ON_10mA_Ronly", f, 60)], mom[("OFF_0V_Conly", f, 60)]
        d = lambda h, m: [float((h - m).real), float((h - m).imag)]  # noqa: E731
        cal["delta_to_mom"][f"{f:g}"] = {"R": {"xx": d(y_r[i, 1], r_m[0]), "yy": d(y_r[i, 0], r_m[1])}, "C": {"xx": d(y_c[i, 1], c_m[0]), "yy": d(y_c[i, 0], c_m[1])}}
    Path(dst_cal).write_text(json.dumps(cal, indent=1), encoding="utf-8")
    for f, v in cal["vs_mom"].items():
        print(f, {k: round(x, 3) for k, x in v.items() if not isinstance(x, dict)}, "alpha hfss/mom", round(v["hfss"]["alpha_exact"], 2), round(v["mom"]["alpha_exact"], 2))
    return cal


if __name__ == "__main__" and len(sys.argv) > 2:
    calibration(sys.argv[2], sys.argv[1])
