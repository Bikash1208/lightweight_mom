"""Phase 20 R4 item 5: library metrics from the assembled CP-2 library (corrected table): per (f, theta, phi) the R3 target quantities, aggregated worst-case over phi, usable
bands, LP 1-bit codes, plots. MODEL_ONLY; the ~10 deg model uncertainty is not in these numbers.

    python scripts/research/phase20_r4_metrics.py

Quantities (handedness-preserving r_RR of the four states, which differ only by the common frame phase so the ladder is read directly):
  lattice deviation   best-offset max deviation of the four phases from the 0/90/180/270 ladder IN STATE ORDER 1-2-3-4 (target <= 15 deg)
  step error          max |cyclic step - 90| over the four steps 1->2->3->4->1 (target <= 20 deg)
  min amplitude, spread   (>= -1.5 dB, <= 1 dB);   opposite-handed   max over states of max(|r_LR|, |r_RL|) (<= -15 dB)
LP 1-bit codes (phi = 0, 90, 180, 270 only, where the incident linear polarisation lies along a structure axis): co-LP code = states 1 vs 3 (wanted: co-polar sign
reversal), cross-LP code = states 2 vs 4 (wanted: cross-polar), for TM and TE incidence: phase difference, amplitude, unwanted-term level.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LIB = ROOT / "results" / "phase20_r4" / "cp2_library_v1.npz"
OUT = ROOT / "results" / "phase20_r4"
TARGETS = {"dev": 15.0, "step": 20.0, "amp_min": -1.5, "spread": 1.0, "opp": -15.0}


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def point_metrics(c):
    """`c[state, 2, 2]` structure-frame circular reflection of the four states at one (f, theta, phi)."""
    rr = c[:, 0, 0]
    ph = np.degrees(np.angle(rr))
    d = np.exp(1j * np.radians(ph - 90.0 * np.arange(4)))
    off = np.angle(d.sum())
    dev = float(np.max(np.abs(np.degrees(np.angle(d * np.exp(-1j * off))))))
    steps = np.array([(ph[(k + 1) % 4] - ph[k]) % 360.0 for k in range(4)])
    step_err = float(np.max(np.abs(steps - 90.0)))
    amp = 20 * np.log10(np.abs(rr))
    opp = 20 * np.log10(np.max(np.abs(np.concatenate([c[:, 1, 0], c[:, 0, 1]]))))
    return {
        "dev": dev,
        "step": step_err,
        "amp_min": float(amp.min()),
        "spread": float(amp.max() - amp.min()),
        "opp": float(opp),
    }


def met(m, with_opp=True):
    ok = (
        m["dev"] <= TARGETS["dev"]
        and m["step"] <= TARGETS["step"]
        and m["amp_min"] >= TARGETS["amp_min"]
        and m["spread"] <= TARGETS["spread"]
    )
    return ok and (m["opp"] <= TARGETS["opp"] if with_opp else True)


def band(freqs, flags, f0=10.0):
    i0 = int(np.argmin(np.abs(freqs - f0)))
    if not flags[i0]:
        return 0.0, None
    lo = hi = i0
    while lo > 0 and flags[lo - 1]:
        lo -= 1
    while hi < len(freqs) - 1 and flags[hi + 1]:
        hi += 1
    return float(freqs[hi] - freqs[lo]), (float(freqs[lo]), float(freqs[hi]))


def lp_codes(lib, ip):
    """Phase difference [deg], amplitudes [dB] and unwanted-term level (re wanted) of the two 1-bit LP codes for TM and TE incidence at phi index ip."""
    s = lib["s_corr"][:, :, ip]  # (f, theta, state, 2, 2) [rx (TE,TM), inc (TE,TM)]
    out = {}
    for inc_name, inc in (("TE", 0), ("TM", 1)):
        co_i, cr_i = inc, 1 - inc  # co-polar row = same index as incident; cross = the other
        for code, (a, b), wanted_cross in (("co", (0, 2), False), ("cross", (1, 3), True)):
            wi, ui = (cr_i, co_i) if wanted_cross else (co_i, cr_i)
            wa, wb = s[:, :, a, wi, inc], s[:, :, b, wi, inc]
            ua, ub = s[:, :, a, ui, inc], s[:, :, b, ui, inc]
            dphi = np.abs(wrap(np.degrees(np.angle(wa / wb))))
            amp = 20 * np.log10(np.stack([np.abs(wa), np.abs(wb)]))
            unwanted = 20 * np.log10(np.maximum(np.abs(ua) / np.abs(wa), np.abs(ub) / np.abs(wb)))
            out[f"{inc_name}_{code}"] = {"dphase": dphi, "amp_min": amp.min(axis=0), "unwanted": unwanted}
    return out


def main():
    z = np.load(LIB)
    lib = {k: z[k] for k in z.files}
    f, th, phi = lib["f_ghz"], lib["theta_deg"], lib["phi_deg"]
    c = lib["c_corr"]
    LP_PHI_IDX = [int(np.argmin(np.abs(phi - p))) for p in (0.0, 90.0, 180.0, 270.0)]  # incident LP along a structure axis
    margin = lib[
        "higher_order_margin"
    ]  # (f, theta, phi): < 1 means a higher diffraction order propagates -> point excluded (never "met")
    M = {k: np.zeros((len(f), len(th), len(phi))) for k in TARGETS}
    for i in range(len(f)):
        for j in range(len(th)):
            for p in range(len(phi)):
                m = point_metrics(c[i, j, p])
                for k in TARGETS:
                    M[k][i, j, p] = (
                        m[k]
                        if margin[i, j, p] >= 1.0
                        else {"dev": 360.0, "step": 360.0, "amp_min": -99.0, "spread": 99.0, "opp": 0.0}[k]
                    )
    worst = {
        "dev": M["dev"].max(2),
        "step": M["step"].max(2),
        "amp_min": M["amp_min"].min(2),
        "spread": M["spread"].max(2),
        "opp": M["opp"].max(2),
    }
    rows = []
    for j, t in enumerate(th):
        flags_all = [met({k: worst[k][i, j] for k in TARGETS}) for i in range(len(f))]
        flags_noopp = [met({k: worst[k][i, j] for k in TARGETS}, with_opp=False) for i in range(len(f))]
        b_all, r_all = band(f, np.array(flags_all))
        b_no, r_no = band(f, np.array(flags_noopp))
        i10 = int(np.argmin(np.abs(f - 10.0)))
        rows.append(
            {
                "theta_deg": float(t),
                "label": "EXTRAPOLATED_BEYOND_R3_TARGETS" if t > 30 else "within_R3_targets",
                "at_10GHz_worst_over_phi": {k: float(worst[k][i10, j]) for k in TARGETS},
                "band_all_targets_GHz": b_all,
                "band_all_targets_range": r_all,
                "band_without_opposite_handed_GHz": b_no,
                "band_without_opposite_handed_range": r_no,
            }
        )
    lp = {}
    for ip in LP_PHI_IDX:
        lp[f"phi{phi[ip]:g}"] = {
            name: {k: v[[int(np.argmin(np.abs(f - 10.0)))]].tolist() for k, v in d.items()}
            for name, d in lp_codes(lib, ip).items()
        }
    lp_at = {}
    for ip in LP_PHI_IDX:
        i10 = int(np.argmin(np.abs(f - 10.0)))
        for name, d in lp_codes(lib, ip).items():
            for j, t in enumerate(th):
                if t in (0.0, 30.0):
                    lp_at[f"phi{phi[ip]:g}_theta{t:g}_{name}"] = {k: float(v[i10, j]) for k, v in d.items()}
    (OUT / "library_metrics.json").write_text(
        json.dumps({"targets": TARGETS, "by_theta": rows, "lp_codes_10GHz": lp_at}, indent=1)
    )
    plots(lib, M, worst)
    for r in rows:
        print(
            r["theta_deg"],
            {k: round(v, 2) for k, v in r["at_10GHz_worst_over_phi"].items()},
            "band(all)",
            r["band_all_targets_GHz"],
            "band(no opp)",
            r["band_without_opposite_handed_GHz"],
        )


def plots(lib, M, worst):
    f, th, phi = lib["f_ghz"], lib["theta_deg"], lib["phi_deg"]
    c = lib["c_corr"]
    (OUT / "plots").mkdir(exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for col, t in enumerate((0.0, 30.0)):
        j = int(np.argmin(np.abs(th - t)))
        for k in range(4):
            r = c[:, j, 0, k, 0, 0]
            ax[0, col].plot(f, np.degrees(np.unwrap(np.angle(r))), label=f"state {k + 1}")
            ax[1, col].plot(f, 20 * np.log10(np.abs(r)), label=f"state {k + 1}")
        ax[0, col].set_title(f"r_RR phase (structure frame), theta {t:g} deg, phi 0")
        ax[1, col].set_title("r_RR level [dB]")
        ax[1, col].axhline(-1.5, color="k", ls=":")
        ax[1, col].set_xlabel("f [GHz]")
    ax[0, 0].legend()
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "state_response_vs_f.png", dpi=110)
    plt.close(fig)
    i10 = int(np.argmin(np.abs(f - 10.0)))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for a, key, title in (
        (ax[0], "opp", "opposite-handed level [dB], worst state"),
        (ax[1], "dev", "ladder deviation [deg]"),
    ):
        im = a.pcolormesh(phi, th, M[key][i10], shading="auto")
        a.set_xlabel("phi [deg]")
        a.set_ylabel("theta [deg]")
        a.set_title(f"{title}, 10 GHz")
        fig.colorbar(im, ax=a)
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "maps_theta_phi_10GHz.png", dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4))
    for t in (0.0, 15.0, 30.0, 45.0):
        j = int(np.argmin(np.abs(th - t)))
        ax.plot(f, worst["dev"][:, j], label=f"theta {t:g}")
    ax.axhline(15, color="k", ls=":")
    ax.set_xlabel("f [GHz]")
    ax.set_ylabel("ladder deviation, worst over phi [deg]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "ladder_deviation_vs_f.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
