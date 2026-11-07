"""Report for the sub-6 GHz split-disc cell scan: plots + results/sub6_splitdisc/CELL_SCAN.md from results/sub6_splitdisc/cell_scan.json (run `sub6_splitdisc_cell_scan.py collect` first).
MODEL_ONLY, IDEAL_BIAS, NOT_HFSS_VALIDATED."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sub6_splitdisc"
BLUE = "#2a78d6"
ORANGE_RAMP = {0.6: "#f6b79c", 0.7: "#f08e66", 0.8: "#eb6834", 0.9: "#a8431c"}  # T: q = r2/r1, one hue light -> dark
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
FAM_ROWS = [("BAR50", "BAR50-02V (primary)"), ("MA", "MA4AGP907 (reference)"), ("SMP", "SMP1345 (fallback)")]
STATE_PAIR = {"BAR50": ("OFF_BAR50", "ON_BAR50"), "MA": ("OFF_MA", "ON_MA"), "SMP": ("OFF_SMP", "ON_SMP")}
HS = (1.575, 3.175)


def cx(v):
    return complex(v[0], v[1])


def series(d, fam, variant, q, h, fn):
    pts = []
    for g in d["geometries"]:
        if g["tag"] == "grid" and g["variant"] == variant and abs(g["q"] - q) < 1e-9 and abs(g["h_mm"] - h) < 1e-9:
            v = fn(g, fam)
            if v is not None:
                pts.append((g["r1_mm"], v))
    return sorted(pts)


def panel_fig(d, fname, title, fn, ylabel, hline=None):
    fig, axs = plt.subplots(3, 2, figsize=(10, 9), sharex=True)
    for i, (fam, famname) in enumerate(FAM_ROWS):
        for j, h in enumerate(HS):
            ax = axs[i, j]
            ax.set_facecolor("#fcfcfb")
            ax.grid(color=GRID, lw=0.6)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if hline is not None:
                ax.axhline(hline, color=MUTED, lw=0.8, ls=":")
            pts = series(d, fam, "S", 1.0, h, fn)
            if pts:
                ax.plot(*zip(*pts), color=BLUE, lw=2, marker="o", ms=4, label="S")
            if fam == "BAR50":
                for q, c in ORANGE_RAMP.items():
                    pts = series(d, fam, "T", q, h, fn)
                    if pts:
                        ax.plot(*zip(*pts), color=c, lw=1.4, marker="o", ms=3, label=f"T q={q}")
            ax.set_title(f"{famname}, h = {h} mm", fontsize=9, color=INK, loc="left")
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9, color=MUTED)
            if i == 2:
                ax.set_xlabel("r (r1) [mm]", fontsize=9, color=MUTED)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle(f"{title}   [MODEL_ONLY, IDEAL_BIAS, NOT_HFSS_VALIDATED]", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / fname, dpi=130)
    plt.close(fig)


def delta_of(state_idx):
    def f(g, fam):
        st = STATE_PAIR[fam][state_idx]
        e = g["states"].get(st)
        return cx(e["Delta_mS"]).imag if e else None

    return f


def rho_re(g, fam):
    r = g["derived"].get(fam)
    return r["rho_re"] if r else None


def beta_mean(state_idx):
    def f(g, fam):
        e = g["states"].get(STATE_PAIR[fam][state_idx])
        return 0.5 * (e["beta_par_over_k0"] + e["beta_perp_over_k0"]) if e else None

    return f


def candidates(d):
    rows = []
    for g in d["geometries"]:
        if g["tag"] != "grid":
            continue
        for fam in ("BAR50", "MA", "SMP"):
            r = g["derived"].get(fam)
            if not r:
                continue
            kind = []
            if r["rho_abs"] <= 0.2:
                kind.append("B-cell |rho|<=0.2")
            if 0.8 <= r["rho_re"] <= 1.2 and abs(r["rho_im"]) <= 0.2:
                kind.append("A-cell 0.8<=Re rho<=1.2")
            if r["rho_re"] <= -0.5:
                kind.append("sign-flip Re rho<=-0.5")
            if kind:
                rows.append((fam, g["variant"], g["h_mm"], g["r1_mm"], g["q"], kind, r))
    return rows


def span(d):
    out = {}
    for g in d["geometries"]:
        if g["tag"] != "grid":
            continue
        for fam, r in g["derived"].items():
            out.setdefault((fam, g["variant"], g["h_mm"]), []).append((r["rho_re"], g["r1_mm"], g["q"]))
    return out


def main():
    d = json.loads((OUT / "cell_scan.json").read_text(encoding="utf-8"))
    panel_fig(d, "plot_delta_off.png", "Im Delta_OFF = Im(Y_par - Y_perp), OFF state", delta_of(0), "Im Delta [mS]")
    panel_fig(d, "plot_delta_on.png", "Im Delta_ON = Im(Y_par - Y_perp), ON state", delta_of(1), "Im Delta [mS]")
    panel_fig(d, "plot_rho.png", "Re rho = Re(Delta_ON / Delta_OFF)", rho_re, "Re rho", hline=0.0)
    panel_fig(d, "plot_beta_off.png", "beta/k0 (mean of par and perp), OFF uniform lattice", beta_mean(0), "beta/k0")
    panel_fig(d, "plot_beta_on.png", "beta/k0 (mean of par and perp), ON uniform lattice", beta_mean(1), "beta/k0")
    g = d["gates"]
    L = ["# CELL_SCAN: sub-6 GHz split-disc cell scan (stage 1)", "", "Labels: **MODEL_ONLY, IDEAL_BIAS** (diode = lumped load only), **NOT_HFSS_VALIDATED**. 5.8 GHz, p = 10 mm, RT5880, grounded, mesh 64, F = 2n.", ""]
    L += ["## Gates", "", f"- `SYMMETRY`: **{g['SYMMETRY']['state']}**, by variant {g['SYMMETRY']['by_variant']}; worst |Y_off|/|Delta| S {g['SYMMETRY']['worst_by_variant']['S']:.2e}, T {g['SYMMETRY']['worst_by_variant']['T']:.2e} (tolerance {g['SYMMETRY']['tolerance']}); {g['SYMMETRY']['n_violations']} of {g['SYMMETRY']['geometry_states_checked']} geometry-states violate."]
    rot = g["ROTATION"]
    L.append(f"- `ROTATION`: **{rot['state']}**, max gap-frame deviation {rot['max_error']} (tolerance {rot['tolerance']})")
    cv = g["CONVERGENCE"]
    L.append(f"- `CONVERGENCE`: **{cv['state']}**, {cv.get('n_failing')} failing rows of {len(cv['rows'])} (beta {cv.get('tolerances', {}).get('beta')}, Delta {cv.get('tolerances', {}).get('Delta')})")
    rs = g["RESONANCE"]
    L.append(f"- `RESONANCE`: **{rs['state']}**; scope: {rs['scope']}; {len(rs['flags'])} flags.")
    L += ["- `RHO_SPAN`: below (report only).", ""]
    L += ["## RHO_SPAN: Re(rho) over the r grid (min / max)", "", "| diode | variant | h [mm] | min Re rho | at | max Re rho | at |", "|---|---|---|---|---|---|---|"]
    for (fam, var, h), v in sorted(span(d).items()):
        lo, hi = min(v), max(v)
        L.append(f"| {fam} | {var} | {h} | {lo[0]:.2f} | r={lo[1]}, q={lo[2]} | {hi[0]:.2f} | r={hi[1]}, q={hi[2]} |")
    c = candidates(d)
    L += ["", f"## Candidates ({len(c)} geometry-diode entries)", "", "| diode | var | h | r1 | q | class | Re rho | Im rho | d_beta par/perp | alpha_ON par/perp [Np/m] |", "|---|---|---|---|---|---|---|---|---|---|"]
    for fam, var, h, r, q, kind, x in c:
        L.append(f"| {fam} | {var} | {h} | {r} | {q} | {'; '.join(kind)} | {x['rho_re']:.2f} | {x['rho_im']:.2f} | {x['delta_beta_par']:+.4f} / {x['delta_beta_perp']:+.4f} | {x['alpha_on_par']:.3f} / {x['alpha_on_perp']:.3f} |")
    L += ["", "## Plots", "", "`plot_delta_off.png`, `plot_delta_on.png`, `plot_rho.png`, `plot_beta_off.png`, `plot_beta_on.png` (one panel per h and diode; blue S, orange ramp T by q).", ""]
    (OUT / "CELL_SCAN_tables.md").write_text("\n".join(L), encoding="utf-8", newline="\n")
    print("\n".join(L[:30]))


if __name__ == "__main__":
    sys.exit(main())
