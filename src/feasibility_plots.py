"""Phase 8A Sec. 43: the 6 required plots. Plots 1-3 are library-only (Mode A, always produced);
plots 4-6 need an actual required-phase distribution (Mode B) and are produced only when one is
supplied -- Sec. 40/43's "only create/generate design-specific [files/plots] when design-specific
inputs were supplied".
"""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 8A plots (src/feasibility_plots.py); "
            "install it with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def _grid_from_table(table, l_values_mm, w_values_mm, value_key):
    l_axis = np.array(sorted(l_values_mm))
    w_axis = np.array(sorted(w_values_mm))
    grid = np.full((len(w_axis), len(l_axis)), np.nan)
    l_index = {round(v, 9): i for i, v in enumerate(l_axis)}
    w_index = {round(v, 9): i for i, v in enumerate(w_axis)}
    for row in table:
        li, wi = l_index.get(round(row["L_mm"], 9)), w_index.get(round(row["W_mm"], 9))
        if li is not None and wi is not None:
            grid[wi, li] = row[value_key]
    return l_axis, w_axis, grid


def plot_phase_vs_lw(table, l_values_mm, w_values_mm, frequency_ghz, path):
    """Plot 1: phase vs L/W."""
    plt = _require_matplotlib()
    l_axis, w_axis, grid = _grid_from_table(table, l_values_mm, w_values_mm, "phase_deg")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    mesh = ax.pcolormesh(l_axis, w_axis, grid, shading="nearest", cmap="twilight", vmin=-180, vmax=180)
    fig.colorbar(mesh, ax=ax, label="phase(S11xx) [deg]")
    ax.set_xlabel("L [mm]")
    ax.set_ylabel("W [mm]")
    ax.set_title(f"Library S11xx phase vs geometry at {frequency_ghz:g} GHz")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_circular_phase_coverage(table, coverage, frequency_ghz, path):
    """Plot 2: circular phase coverage."""
    plt = _require_matplotlib()
    phases_rad = np.radians([row["phase_deg"] for row in table])
    fig, ax = plt.subplots(figsize=(6.4, 7.2), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    pts = ax.scatter(phases_rad, np.ones_like(phases_rad), s=14, zorder=3, label="library phase (S11xx)")
    gap_start = np.radians(coverage.largest_gap_before_deg)
    gap_end = gap_start + np.radians(coverage.largest_gap_deg)
    theta_fill = np.linspace(gap_start, gap_end, 100)
    fill = ax.fill_between(
        theta_fill, 0, 1, color="red", alpha=0.2, zorder=1, label=f"largest gap: {coverage.largest_gap_deg:.1f} deg"
    )
    ax.set_yticklabels([])
    fig.suptitle(f"Library phase coverage at {frequency_ghz:g} GHz", y=0.97)
    fig.text(
        0.5,
        0.08,
        f"{coverage.n_points} cells | coverage = {coverage.coverage_deg:.1f} deg of 360 deg",
        ha="center",
        fontsize=9,
    )
    fig.legend(handles=[pts, fill], loc="lower center", bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.16, 1, 0.93))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_coverage_vs_frequency(freq_coverage, path):
    """Plot 3: phase coverage vs frequency."""
    plt = _require_matplotlib()
    freqs = [r["frequency_GHz"] for r in freq_coverage]
    coverage = [r["coverage_deg"] for r in freq_coverage]
    gap = [r["largest_gap_deg"] for r in freq_coverage]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(freqs, coverage, label="coverage [deg]")
    ax.plot(freqs, gap, label="largest gap [deg]", linestyle="--")
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("degrees")
    ax.set_ylim(0, 360)
    ax.set_title("Circular phase coverage vs frequency (S11xx)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_required_phase_map(x_mm, y_mm, phase_deg, path):
    """Plot 4: required holographic phase (design-specific only)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    sc = ax.scatter(x_mm, y_mm, c=phase_deg, cmap="twilight", vmin=-180, vmax=180, s=20)
    fig.colorbar(sc, ax=ax, label="required phase [deg]")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title("Required holographic phase")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_required_phase_histogram(phase_deg, path):
    """Plot 4b: required-phase histogram (design-specific only)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.hist(phase_deg, bins=36, range=(-180, 180))
    ax.set_xlabel("required phase [deg]")
    ax.set_ylabel("count")
    ax.set_title("Required holographic phase histogram")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_matching_error(match_results, path):
    """Plots 5-6: required vs nearest validated library phase, and the resulting phase error
    (design-specific only)."""
    plt = _require_matplotlib()
    errors = [m.phase_error_deg for m in match_results if m.phase_error_deg is not None]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    if errors:
        ax.hist(errors, bins=30)
    ax.set_xlabel("phase error (nearest library - required) [deg]")
    ax.set_ylabel("count")
    n_unmatched = sum(1 for m in match_results if not m.matched)
    ax.set_title(f"Required-vs-matched phase error ({n_unmatched} unmatched of {len(match_results)})")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
