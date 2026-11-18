"""MODEL_ONLY 11.6-mm binary dither of two CP sectors and B-OFF state."""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import maximum_filter

from scripts.research.sub6_meander_aperture import C0, F, K0, THETA_BEAM, angular_distance

P = .0116
R = .150
ALPHA_OFF = .13
ALPHAS_ON = (1.16, 1.33, 1.48)
WIDTHS_MM = (30, 40)
SEEDS = tuple(range(8))
ANG_STEP = 2
FEED_FLOOR = P / 2


def lattice():
    ix, iy = np.meshgrid(np.arange(-12, 13), np.arange(-12, 13), indexing="xy")
    x, y = ix * P, iy * P
    keep = (x*x + y*y <= R*R) & ((ix != 0) | (iy != 0))  # reserve centre for the feed
    return ix[keep], iy[keep], x[keep], y[keep]


def dither(ix, iy, x, width_mm, seed):
    """Column-stratified stochastic dither of p_B(x); no hand mixing within a cell."""
    p_b = .5 * (1 - np.tanh(x / (width_mm * 1e-3)))
    is_b = np.zeros(len(x), bool)
    rng = np.random.default_rng(seed)
    for col in np.unique(ix):
        ids = np.flatnonzero(ix == col)
        exact = float(np.sum(p_b[ids]))
        n_b = int(math.floor(exact) + (rng.random() < exact % 1))
        is_b[rng.choice(ids, size=n_b, replace=False)] = True
    return p_b, is_b


def ray_integrals(ix, iy, x, y, is_b, alpha_a, alpha_b):
    """Midpoint integration through nearest assigned cells along each feed-to-cell ray."""
    n = int(math.ceil(R / (P/4)))
    rho = np.hypot(x, y)
    s = (np.arange(n) + .5) / n
    sx, sy = x[:, None] * s, y[:, None] * s
    nmax = 12
    grid = np.zeros((2*nmax+1, 2*nmax+1), bool)
    grid[iy+nmax, ix+nmax] = is_b
    jj = np.clip(np.rint(sy/P).astype(int) + nmax, 0, 2*nmax)
    ii = np.clip(np.rint(sx/P).astype(int) + nmax, 0, 2*nmax)
    local = np.where(grid[jj, ii], alpha_b, alpha_a)
    local[(jj == nmax) & (ii == nmax)] = 0.0
    return rho * np.mean(local, axis=1)


def conversion(is_b_grid, alpha_a, alpha_b):
    """Radial launched-power conversion assigned to A and B cells."""
    nmax = 12
    angles = np.radians(np.arange(360))
    s = np.arange(.0005, R, .001)
    sx = np.cos(angles)[:, None] * s
    sy = np.sin(angles)[:, None] * s
    jj = np.clip(np.rint(sy/P).astype(int) + nmax, 0, 2*nmax)
    ii = np.clip(np.rint(sx/P).astype(int) + nmax, 0, 2*nmax)
    is_b = is_b_grid[jj, ii]
    a = np.where(is_b, alpha_b, alpha_a)
    a[(jj == nmax) & (ii == nmax)] = 0.0
    optical_depth_before = np.cumsum(a, axis=1) * .001 - a * .001
    segment_conversion = np.exp(-2 * optical_depth_before) * (1 - np.exp(-2 * a * .001))
    out_a = float(np.mean(np.sum(np.where(is_b, 0, segment_conversion), axis=1)))
    out_b = float(np.mean(np.sum(np.where(is_b, segment_conversion, 0), axis=1)))
    return {"A": out_a, "B": out_b, "total": out_a + out_b}


def directions(step):
    theta_d = np.arange(0, 90+step, step)
    phi_d = np.arange(0, 360, step)
    tt, pp = np.meshgrid(np.radians(theta_d), np.radians(phi_d), indexing="ij")
    t, p = tt.ravel(), pp.ravel()
    uv = np.column_stack([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p)])
    return tt.shape, t, p, uv


def pattern_metrics(fields_a, fields_b, shape, t, p, step):
    ia, ib = abs(fields_a)**2, abs(fields_b)**2
    domega = math.radians(step)**2 * np.sin(t)
    power = float(np.sum((ia+ib)*domega))
    first_null = math.asin(min(1, 1.22*(C0/F)/(2*R/math.sqrt(2))))
    results = {}
    for hand, target_phi, co, cross in (("A_RHCP", 0.0, ia, ib), ("B_LHCP", math.pi, ib, ia)):
        distance = angular_distance(t, p, target_phi)
        in_beam = distance <= math.radians(20)
        target = np.argmin(abs(t-math.radians(30))+abs((p-target_phi+math.pi)%(2*math.pi)-math.pi))
        if not np.any(co):
            results[hand] = {"peak_theta_deg": None, "peak_phi_deg": None, "directivity_dbi": None,
                             "sll_db": None, "sll_theta_deg": None, "sll_phi_deg": None,
                             "target_intensity": 0.0, "target_cross_intensity": float(cross[target]),
                             "target_xpd_db": None}
            continue
        peak = np.flatnonzero(in_beam)[int(np.argmax(co[in_beam]))]
        local_max = co.reshape(shape) >= maximum_filter(co.reshape(shape), size=(3,3), mode=("nearest","wrap"))*(1-1e-12)
        outer = (distance > first_null) & local_max.ravel()
        sidelobe = np.flatnonzero(outer)[int(np.argmax(co[outer]))]
        results[hand] = {"peak_theta_deg": float(np.degrees(t[peak])), "peak_phi_deg": float(np.degrees(p[peak])),
                         "directivity_dbi": 10*math.log10(4*math.pi*co[peak]/power),
                         "sll_db": 10*math.log10(co[sidelobe]/co[peak]),
                         "sll_theta_deg": float(np.degrees(t[sidelobe])), "sll_phi_deg": float(np.degrees(p[sidelobe])),
                         "target_intensity": float(co[target]), "target_cross_intensity": float(cross[target]),
                         "target_xpd_db": 10*math.log10(co[target]/cross[target]) if cross[target]>0 else math.inf}
    return results


def summarize(rows):
    keys = ("two_beam_A_directivity_dbi", "two_beam_B_directivity_dbi", "two_beam_worst_sll_db",
            "two_beam_B_minus_A_db", "one_beam_A_change_db", "one_beam_B_suppression_db",
            "one_beam_B_to_A_db", "two_beam_conversion", "one_beam_conversion",
            "two_beam_A_conversion", "one_beam_A_conversion", "two_beam_B_conversion", "one_beam_B_conversion",
            "two_beam_A_gain_dbi_ideal", "one_beam_A_gain_dbi_ideal")
    out = {}
    for key in keys:
        a = sorted(float(r[key]) for r in rows if r[key] is not None)
        out[key] = {"min": a[0], "median": float(np.median(a)), "max": a[-1]} if a else {"min": None, "median": None, "max": None}
    return out


def run(outdir, off_anisotropy_scale=1.0):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ix, iy, x, y = lattice()
    shape, t, phi, uv = directions(ANG_STEP)
    phase = np.exp(1j*K0*(uv @ np.column_stack([x,y]).T))*np.cos(t)[:,None]
    launch_a = np.exp(-1j*K0*math.sin(THETA_BEAM)*x)
    launch_b = np.exp(+1j*K0*math.sin(THETA_BEAM)*x)
    rows = []
    for width in WIDTHS_MM:
        for seed in SEEDS:
            p_b, is_b = dither(ix, iy, x, width, seed)
            grid = np.zeros((25,25), bool)
            grid[iy+12,ix+12] = is_b
            for alpha_on in ALPHAS_ON:
                states = {}
                for state, alpha_b in (("two_beam",alpha_on),("one_beam",ALPHA_OFF)):
                    rho = np.hypot(x,y)
                    integral = ray_integrals(ix,iy,x,y,is_b,alpha_on,alpha_b)
                    local_alpha = np.where(is_b,alpha_b,alpha_on)
                    amp = np.sqrt(2*local_alpha/np.maximum(rho,FEED_FLOOR))*np.exp(-integral)
                    fa = phase @ (amp*(~is_b)*launch_a)
                    b_coherent_scale = 1.0 if state == "two_beam" else off_anisotropy_scale
                    fb = phase @ (amp*is_b*launch_b*b_coherent_scale)
                    conv = conversion(grid,alpha_on,alpha_b)
                    states[state] = {"metrics":pattern_metrics(fa,fb,shape,t,phi,ANG_STEP),
                                     "conversion":conv["total"], "conversion_A":conv["A"], "conversion_B":conv["B"],
                                     "radiated_conversion":conv["A"]+b_coherent_scale**2*conv["B"]}
                two, one = states["two_beam"], states["one_beam"]
                aa2, bb2 = two["metrics"]["A_RHCP"], two["metrics"]["B_LHCP"]
                aa1, bb1 = one["metrics"]["A_RHCP"], one["metrics"]["B_LHCP"]
                rows.append({"width_mm":width,"seed":seed,"alpha_on_np_per_m":alpha_on,
                             "alpha_b_one_beam_np_per_m":ALPHA_OFF,"n_cells":len(x),"n_A":int(sum(~is_b)),"n_B":int(sum(is_b)),
                             "column_probability_formula":"p_B(x)=(1-tanh(x/w))/2; count-stratified seeded stochastic assignment",
                             "two_beam":two,"one_beam":one,
                             "two_beam_A_directivity_dbi":aa2["directivity_dbi"],
                             "two_beam_B_directivity_dbi":bb2["directivity_dbi"],
                             "two_beam_worst_sll_db":max(aa2["sll_db"],bb2["sll_db"]),
                             "two_beam_B_minus_A_db":10*math.log10(bb2["target_intensity"]/aa2["target_intensity"]),
                             "one_beam_A_change_db":10*math.log10(aa1["target_intensity"]/aa2["target_intensity"]),
                             "one_beam_B_suppression_db":(10*math.log10(bb2["target_intensity"]/bb1["target_intensity"])
                                                           if bb1["target_intensity"]>0 else None),
                             "one_beam_B_to_A_db":(10*math.log10(bb1["target_intensity"]/aa1["target_intensity"])
                                                   if bb1["target_intensity"]>0 else None),
                             "one_beam_B_coherent_zero":bb1["target_intensity"]==0,
                             "two_beam_conversion":two["conversion"],"one_beam_conversion":one["conversion"],
                             "two_beam_A_conversion":two["conversion_A"],"one_beam_A_conversion":one["conversion_A"],
                             "two_beam_B_conversion":two["conversion_B"],"one_beam_B_conversion":one["conversion_B"],
                             "two_beam_A_gain_dbi_ideal":aa2["directivity_dbi"]+10*math.log10(two["radiated_conversion"]),
                             "one_beam_A_gain_dbi_ideal":aa1["directivity_dbi"]+10*math.log10(one["radiated_conversion"])})
            if seed == 0:
                fig,ax=plt.subplots(figsize=(5,5),layout="constrained")
                ax.scatter(x[~is_b]*1e3,y[~is_b]*1e3,s=30,c="#d95f02",label="A RHCP")
                ax.scatter(x[is_b]*1e3,y[is_b]*1e3,s=30,c="#1b9e77",label="B LHCP")
                ax.set(xlabel="x (mm)",ylabel="y (mm)",title=f"{width} mm stochastic dither; 11.6 mm cell grid")
                ax.set_aspect("equal");ax.legend()
                fig.savefig(outdir/f"layout_w{width}_seed0.png",dpi=170);plt.close(fig)
    summary = {f"w{w}_alpha{a:.2f}":summarize([r for r in rows if r["width_mm"]==w and r["alpha_on_np_per_m"]==a])
               for w in WIDTHS_MM for a in ALPHAS_ON}
    result={"label":"MODEL_ONLY / NOT_HFSS_VALIDATED","frequency_hz":F,"pitch_mm":P*1e3,"radius_mm":R*1e3,
            "centre_cell":"reserved for feed","soft_widths_mm":WIDTHS_MM,"dither_seeds":SEEDS,
            "alpha_on_range_np_per_m":ALPHAS_ON,"alpha_B_off_np_per_m":ALPHA_OFF,
            "phase_assumption":("B-OFF coherent LHCP B beam amplitude set to zero; OFF attenuation retained with unspecified destination"
                                if off_anisotropy_scale == 0 else
                                "B-OFF residual radiation retains coherent LHCP B phase ramp (conservative residual beam)"),
            "off_anisotropy_coherent_scale":off_anisotropy_scale,
            "radial_model":"rho^-1/2 sqrt(2alpha_local) exp(-integral along nearest-cell radial path alpha ds)",
            "angle_step_deg":ANG_STEP,"summary":summary,"rows":rows}
    (outdir/"dither_grid.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), layout="constrained")
    for alpha_on in ALPHAS_ON:
        ss = [summary[f"w{w}_alpha{alpha_on:.2f}"] for w in WIDTHS_MM]
        for ax, key, label in ((axes[0], "two_beam_worst_sll_db", "Two-beam SLL (dB)"),
                               (axes[1], "one_beam_B_suppression_db", "B suppression (dB)")):
            if any(s[key]["median"] is None for s in ss):
                continue
            med = np.array([s[key]["median"] for s in ss])
            lo = np.array([s[key]["min"] for s in ss])
            hi = np.array([s[key]["max"] for s in ss])
            ax.errorbar(WIDTHS_MM, med, yerr=[med-lo, hi-med], marker="o", capsize=3,
                        label=f"alpha_ON={alpha_on:.2f}")
            ax.set(xlabel="Soft boundary width (mm)", ylabel=label, xticks=WIDTHS_MM)
            ax.grid(alpha=.25)
    axes[0].axhline(-15, color="black", linestyle="--", linewidth=1)
    axes[1].axhline(10, color="black", linestyle="--", linewidth=1)
    axes[0].legend(fontsize=8)
    fig.suptitle("R150, 11.6 mm stochastic dither: median and 8-seed range")
    fig.savefig(outdir/"dither_summary.png",dpi=180)
    plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == "__main__":
    import sys
    run(sys.argv[1])
