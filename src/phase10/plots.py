"""Sec. 46: the required plots, following ``phase9.plots``'s own style (matplotlib Agg backend,
raise ImportError with an install hint rather than silently skip)."""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 10 plots; install with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


# ------------------------------------------------------------------------------- neighborhood


def plot_neighbor_distance_vs_perturbation(study_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for name, r in study_results.items():
        for label in r.contributions:
            cardinal_mean = sum(c.delta_s_abs for c in r.contributions[label].values()) / 4
            diagonal_mean = sum(c.delta_s_abs for c in r.diagonal_contributions[label].values()) / 4
            long_mean = sum(c.delta_s_abs for c in r.long_range_cases.values()) / len(r.long_range_cases)
            ax.plot(
                [1, 1.41, 2],
                [cardinal_mean, diagonal_mean, long_mean],
                "o-",
                alpha=0.5,
                label=f"{name}/{label}" if label == "strong_high" else None,
            )
    ax.set_xlabel("neighbor distance [periods]")
    ax.set_ylabel("mean |delta S| [A.m]")
    ax.set_xticks([1, 1.41, 2])
    ax.set_xticklabels(["cardinal (1)", "diagonal (1.41)", "long-range (2)"])
    ax.set_yscale("log")
    ax.set_title("neighbor distance vs. response perturbation")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_contrast_vs_complex_error(study_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    labels = ("weak", "moderate", "strong_low", "strong_high")
    for name, r in study_results.items():
        ys = [r.predictive_errors[l].relative_to_s0 for l in labels]
        ax.plot(labels, ys, "o-", label=name)
    ax.set_ylabel("5-cell-model predictive error / |S0|")
    ax.set_yscale("log")
    ax.set_title("neighbor contrast vs. 5-cell predictive (complex) error")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_cardinal_vs_diagonal_sensitivity(study_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    labels = ("weak", "moderate", "strong_low", "strong_high")
    names = list(study_results.keys())
    width = 0.1
    x = np.arange(len(labels))
    for k, name in enumerate(names):
        r = study_results[name]
        ratios = [r.cardinal_vs_diagonal[l].ratio_diagonal_to_cardinal for l in labels]
        ratios = [min(v, 3.0) if np.isfinite(v) else 3.0 for v in ratios]
        ax.bar(x + k * width, ratios, width, label=name)
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xticks(x + width * len(names) / 2)
    ax.set_xticklabels(labels)
    ax.set_ylabel("mean |diagonal delta| / mean |cardinal delta|")
    ax.set_title("cardinal vs. diagonal sensitivity")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_5cell_vs_9cell_error(study_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    names = list(study_results.keys())
    worst = [max(r.predictive_errors[l].relative_to_s0 for l in r.predictive_errors) for r in study_results.values()]
    colors = [
        "crimson" if d == "LARGER_NEIGHBORHOOD_REQUIRED" else ("orange" if d == "9_CELL_REQUIRED" else "seagreen")
        for d in (r.decision.final_decision for r in study_results.values())
    ]
    ax.bar(names, worst, color=colors)
    ax.set_ylabel("worst-case 5-cell predictive error / |S0|")
    ax.set_yscale("log")
    ax.set_title("5-cell vs. 9-cell: worst-case error missed by ignoring diagonals")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_long_range_decay(study_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for name, r in study_results.items():
        vals = (
            list(r.long_range_relative.values())
            if hasattr(r, "long_range_relative")
            else list(r.decision.long_range_relative.values())
        )
        ax.plot(list(r.decision.long_range_relative.keys()), vals, "o-", label=name)
    ax.set_ylabel("relative deviation at distance-2")
    ax.set_yscale("log")
    ax.set_title("long-range (2-period) coupling")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------------- surrogate


def plot_predicted_vs_reference(preds, refs, component, path):
    plt = _require_matplotlib()
    preds = np.asarray([p if isinstance(p, complex) else complex("nan") for p in preds])
    refs = np.asarray(refs, dtype=complex)
    x = getattr(refs, component)
    y = getattr(preds, component)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(x, y, s=10, alpha=0.6)
    lims = [min(np.nanmin(x), np.nanmin(y)), max(np.nanmax(x), np.nanmax(y))]
    ax.plot(lims, lims, "k--", lw=1)
    ax.set_xlabel(f"reference {component}(S)")
    ax.set_ylabel(f"predicted {component}(S)")
    ax.set_title(f"predicted vs. reference {component}(S)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_complex_plane_comparison(preds, refs, path):
    plt = _require_matplotlib()
    preds = [p if isinstance(p, complex) else complex("nan") for p in preds]
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter([p.real for p in refs], [p.imag for p in refs], s=14, label="reference (MoM)")
    ax.scatter([p.real for p in preds], [p.imag for p in preds], s=10, marker="x", label="surrogate")
    ax.set_xlabel("Re(S)")
    ax.set_ylabel("Im(S)")
    ax.set_aspect("equal")
    ax.set_title("complex-plane comparison (blind test)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_distribution(errors, path, xlabel="complex error [A.m]"):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.hist(errors, bins=30)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.set_title("blind-test error distribution")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_vs_geometry(samples, errors, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ls = [s.center_lw_mm[0] for s in samples]
    sc = ax.scatter(ls, errors, c=errors, cmap="viridis", s=14)
    ax.set_xlabel("center L [mm]")
    ax.set_ylabel("complex error [A.m]")
    ax.set_title("error vs. center geometry")
    fig.colorbar(sc, ax=ax, label="error")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_vs_neighbor_contrast(samples, errors, path):
    from .dataset import ALL_POSITIONS

    plt = _require_matplotlib()
    contrasts = []
    for s in samples:
        cl, cw = s.center_lw_mm
        contrasts.append(
            max(abs(s.neighbors_lw_mm[p][0] - cl) + abs(s.neighbors_lw_mm[p][1] - cw) for p in ALL_POSITIONS)
        )
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(contrasts, errors, s=14, alpha=0.6)
    ax.set_xlabel("max neighbor contrast [mm]")
    ax.set_ylabel("complex error [A.m]")
    ax.set_title("error vs. neighbor contrast")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_error_vs_frequency(samples, errors, path):
    """Sec. 46: required, but this project's dataset covers a single design frequency (documented
    scope, ``DomainSpec.frequency_hz``) -- the plot is produced but is necessarily a single point
    per sample at that one frequency, annotated as such rather than a fabricated sweep."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    freqs = [s.frequency_hz / 1e9 for s in samples]
    ax.scatter(freqs, errors, s=14, alpha=0.6)
    ax.text(
        0.5,
        0.9,
        "single design frequency only (dataset scope)",
        ha="center",
        transform=ax.transAxes,
        fontsize=9,
        color="crimson",
        style="italic",
    )
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("complex error [A.m]")
    ax.set_title("error vs. frequency")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_worst_case_error(reports_by_model, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    names = list(reports_by_model.keys())
    means = [reports_by_model[n].mean_complex_error for n in names]
    p95s = [reports_by_model[n].p95_complex_error for n in names]
    maxs = [reports_by_model[n].max_complex_error for n in names]
    x = np.arange(len(names))
    width = 0.25
    ax.bar(x - width, means, width, label="mean")
    ax.bar(x, p95s, width, label="p95")
    ax.bar(x + width, maxs, width, label="max")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("complex error [A.m]")
    ax.set_yscale("log")
    ax.set_title("mean / p95 / worst-case error by model")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------- inverse_validation


def plot_inverse_target_vs_achieved(results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    reached = [r for r in results if r.status == "REACHED"]
    unreachable = [r for r in results if r.status != "REACHED"]
    ax.scatter(
        [r.target.real for r in reached], [r.target.imag for r in reached], s=40, marker="o", label="target (REACHED)"
    )
    ax.scatter(
        [r.mom_achieved_at_candidate.real for r in reached if isinstance(r.mom_achieved_at_candidate, complex)],
        [r.mom_achieved_at_candidate.imag for r in reached if isinstance(r.mom_achieved_at_candidate, complex)],
        s=20,
        marker="x",
        label="MoM-achieved",
    )
    if unreachable:
        ax.scatter(
            [r.target.real for r in unreachable],
            [r.target.imag for r in unreachable],
            s=40,
            marker="s",
            color="crimson",
            label="target (UNREACHABLE)",
        )
    ax.set_xlabel("Re(S)")
    ax.set_ylabel("Im(S)")
    ax.set_aspect("equal")
    ax.set_title("inverse validation: target vs. MoM-achieved response")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------- holographic


def plot_spatial_complex_response_error(cases, grid_n, path):
    plt = _require_matplotlib()
    grid = np.full((grid_n, grid_n), np.nan)
    for c in cases:
        if isinstance(c.px_surrogate, complex):
            grid[c.i, c.j] = abs(c.px_mom - c.px_surrogate)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(grid, cmap="inferno")
    ax.set_title("spatial complex-response error (MoM vs. surrogate)")
    fig.colorbar(im, ax=ax, label="|error| [A.m]")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_aperture_field(cases, grid_n, component, path):
    plt = _require_matplotlib()
    grid = np.full((grid_n, grid_n), np.nan)
    for c in cases:
        grid[c.i, c.j] = abs(c.px_mom) if component == "magnitude" else np.degrees(np.angle(c.px_mom))
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(grid, cmap="viridis" if component == "magnitude" else "twilight")
    ax.set_title(f"aperture {component} (MoM reference)")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_far_field_pattern(far_field, theta_grid_deg, path, title="far-field pattern"):
    plt = _require_matplotlib()
    power_db = 10 * np.log10(np.maximum(np.abs(far_field) ** 2, 1e-30))
    power_db -= np.max(power_db)
    fig, ax = plt.subplots(figsize=(6.5, 5), subplot_kw={"projection": "polar"})
    theta_rad = np.radians(theta_grid_deg)
    c = ax.pcolormesh(
        np.radians(np.linspace(0, 359, power_db.shape[1])),
        theta_grid_deg[:, 0],
        power_db,
        shading="auto",
        vmin=-30,
        vmax=0,
    )
    ax.set_title(title)
    fig.colorbar(c, ax=ax, label="relative power [dB]")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_beam_pointing_comparison(metrics_mom, metrics_sur, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    labels = [f"beam {i}" for i in range(metrics_mom.n_intended_beams)]
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, metrics_mom.pointing_error_deg, width, label="MoM reference")
    if metrics_sur is not None:
        ax.bar(x + width / 2, metrics_sur.pointing_error_deg, width, label="surrogate")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("pointing error [deg]")
    ax.set_title("beam pointing comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_peak_power_comparison(metrics_mom, metrics_sur, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    labels = [f"beam {i}" for i in range(metrics_mom.n_intended_beams)]
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, metrics_mom.matched_peak_power_db, width, label="MoM reference")
    if metrics_sur is not None:
        ax.bar(x + width / 2, metrics_sur.matched_peak_power_db, width, label="surrogate")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("matched peak power [dB, relative]")
    ax.set_title("peak power comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_sidelobe_comparison(far_field_mom, far_field_sur, theta_grid_deg, phi_grid_deg, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    mid = theta_grid_deg.shape[1] // 2
    p_mom = 10 * np.log10(np.maximum(np.abs(far_field_mom[:, mid]) ** 2, 1e-30))
    p_mom -= np.max(p_mom)
    ax.plot(theta_grid_deg[:, mid], p_mom, label="MoM reference")
    if far_field_sur is not None:
        p_sur = 10 * np.log10(np.maximum(np.abs(far_field_sur[:, mid]) ** 2, 1e-30))
        p_sur -= np.max(p_sur)
        ax.plot(theta_grid_deg[:, mid], p_sur, "--", label="surrogate")
    ax.set_xlabel("theta [deg]")
    ax.set_ylabel("relative power [dB]")
    ax.set_title(f"pattern cut at phi={phi_grid_deg[0, mid]:.0f} deg (sidelobe comparison)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
