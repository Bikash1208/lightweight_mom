"""Analysis of the HFSS far fields of board_K2_default (RT5880, 150 mm): state 1 all OFF, state 2 all ON (rE [V], theta 0..180 step 2, phi -180..180 step 2, theta-outer, 10.5 GHz).  Not a model.
Output: results/sidequest_hybrid/hfss_ffd_analysis.json  (power split, peak, sidelobes with directions, common / difference decomposition)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

UP = Path(r"C:\Users\Bikash\.claude\uploads\e5f5ea6a-da86-4c57-8877-0e8d73bc3087")
FILES = {1: UP / "d04320ce-state1_allOFF_exportfields.ffd", 2: UP / "3a400938-state2_allON_exportfields.ffd"}
NT, NP = 91, 181


def load(p):
    rows = []
    with open(p, encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            if i < 4:
                continue
            v = line.split()
            if len(v) >= 4:
                rows.append([float(x) for x in v[:4]])
    a = np.array(rows)
    assert a.shape[0] == NT * NP, a.shape
    et = (a[:, 0] + 1j * a[:, 1]).reshape(NT, NP)
    ep = (a[:, 2] + 1j * a[:, 3]).reshape(NT, NP)
    return et, ep


def main():
    th = np.radians(np.arange(0, 181, 2.0))
    ph = np.radians(np.arange(-180, 181, 2.0))
    TG, PG = np.meshgrid(th, ph, indexing="ij")
    w = np.sin(TG) * np.radians(2) * np.radians(2)
    data = {k: load(p) for k, p in FILES.items()}
    out = {}
    for k, (et, ep) in data.items():
        u = np.abs(et) ** 2 + np.abs(ep) ** 2
        prad = float(np.sum(u[:, :-1] * w[:, :-1]))  # phi -180 and 180 duplicate
        d = 4 * math.pi * u / prad
        ddb = 10 * np.log10(d + 1e-30)
        ipk = np.unravel_index(np.argmax(u), u.shape)
        band = lambda a, b, m=None: float(np.sum((u * w)[(np.degrees(TG) >= a) & (np.degrees(TG) < b)][:, None] if False else 0))  # noqa: E731
        tdeg = np.degrees(TG)
        frac = {}
        for nm, (a, b) in {"theta<30": (0, 30), "30-60": (30, 60), "60-90": (60, 90.001), "back(90-180)": (90.001, 181)}.items():
            m = (tdeg >= a) & (tdeg < b)
            frac[nm] = float(np.sum((u * w)[:, :-1][m[:, :-1]]) / prad)
        # local maxima on the grid (8-neighbourhood, periodic in phi) above -15 dB of the peak
        peaks = []
        for i in range(1, NT - 1):
            for j in range(NP - 1):
                c = u[i, j]
                nb = [u[i + a, (j + b) % (NP - 1)] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
                if c >= max(nb) and c > 10 ** (-1.5) * u.max():
                    peaks.append((float(np.degrees(th[i])), float(np.degrees(ph[j])), float(ddb[i, j]), float(10 * np.log10(c / u.max()))))
        peaks.sort(key=lambda t: -t[2])
        # co / cross at the peak: LP phi-hat design
        xpd = float(10 * np.log10(abs(ep[ipk]) ** 2 / (abs(et[ipk]) ** 2 + 1e-30)))
        out[k] = {"D_peak_dbi": float(ddb[ipk]), "peak_dir": [float(np.degrees(th[ipk[0]])), float(np.degrees(ph[ipk[1]]))], "power_fraction": frac, "peaks_theta_phi_Ddbi_reldB": peaks[:12], "Ephi_over_Etheta_at_peak_db": xpd,
                  "D_at_(25,0)": float(ddb[12, 90]), "D_at_(25,180)": float(ddb[12, 0]), "D_at_(24,176)": float(ddb[12, 2]), "D_at_(20,154)": float(ddb[10, 167 - 90 if False else (154 + 180) // 2])}
    # common / difference decomposition on the upper hemisphere
    (e1t, e1p), (e2t, e2p) = data[1], data[2]
    u1, u2 = np.abs(e1t) ** 2 + np.abs(e1p) ** 2, np.abs(e2t) ** 2 + np.abs(e2p) ** 2
    cm = ((e1t + e2t) / 2, (e1p + e2p) / 2)
    df = ((e1t - e2t) / 2, (e1p - e2p) / 2)
    ucm, udf = np.abs(cm[0]) ** 2 + np.abs(cm[1]) ** 2, np.abs(df[0]) ** 2 + np.abs(df[1]) ** 2
    prad1 = float(np.sum(u1[:, :-1] * w[:, :-1]))
    out["common_mode"] = {"peak_dir": [float(np.degrees(th[i])) for i in [np.unravel_index(np.argmax(ucm), ucm.shape)[0]]] + [float(np.degrees(ph[np.unravel_index(np.argmax(ucm), ucm.shape)[1]]))], "peak_D_dbi_rel_state1_power": float(10 * np.log10(4 * math.pi * ucm.max() / prad1)),
                          "power_ratio_common_to_total": float(np.sum(ucm[:, :-1] * w[:, :-1]) / (np.sum(ucm[:, :-1] * w[:, :-1]) + np.sum(udf[:, :-1] * w[:, :-1])))}
    out["difference_mode"] = {"peak_dir": [float(np.degrees(th[np.unravel_index(np.argmax(udf), udf.shape)[0]])), float(np.degrees(ph[np.unravel_index(np.argmax(udf), udf.shape)[1]]))], "peak_D_dbi_rel_state1_power": float(10 * np.log10(4 * math.pi * udf.max() / prad1))}
    # azimuth structure of the 60-90 ring and of the back hemisphere: power per 30 deg phi sector
    for k, (et, ep) in data.items():
        u = np.abs(et) ** 2 + np.abs(ep) ** 2
        prad = float(np.sum(u[:, :-1] * w[:, :-1]))
        tdeg = np.degrees(TG)
        for nm, (a, b) in {"ring60_90": (60, 90.001), "back": (90.001, 181)}.items():
            m = (tdeg >= a) & (tdeg < b)
            sect = []
            for s0 in range(-180, 180, 30):
                pm = (np.degrees(PG) >= s0) & (np.degrees(PG) < s0 + 30)
                sect.append(round(float(np.sum((u * w)[m & pm]) / prad * 100), 1))
            out[k][f"{nm}_pct_per_30deg_phi_sector_from_-180"] = sect
    Path("results/sidequest_hybrid/hfss_ffd_analysis.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline=chr(10))
    for k in (1, 2):
        o = out[k]
        print("state", k, "D_peak %.2f dBi at (%.0f, %.0f)" % (o["D_peak_dbi"], *o["peak_dir"]), "power", {a: round(100 * b, 1) for a, b in o["power_fraction"].items()}, "Ephi/Etheta %.1f dB" % o["Ephi_over_Etheta_at_peak_db"])
        print("  D at (25,0) %.1f, (25,180) %.1f, (24,176) %.1f" % (o["D_at_(25,0)"], o["D_at_(25,180)"], o["D_at_(24,176)"]))
        print("  peaks (theta, phi, D dBi, rel dB):", [(round(a), round(b), round(c, 1), round(d, 1)) for a, b, c, d in o["peaks_theta_phi_Ddbi_reldB"]])
        print("  ring60_90 pct/30deg sectors from -180:", o["ring60_90_pct_per_30deg_phi_sector_from_-180"])
        print("  back pct/30deg sectors:", o["back_pct_per_30deg_phi_sector_from_-180"])
    print("common", out["common_mode"], "difference", out["difference_mode"])


if __name__ == "__main__":
    main()
