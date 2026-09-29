"""Phase 12B's required plots: aligned vs unaligned convergence, phase vs L per period/frequency,
coverage vs frequency per period, |S11| vs L, duplicate-audit table.
"""

from pathlib import Path


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 12B plots (src/phase12b/plots.py); install it "
            "with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_aligned_vs_unaligned_convergence(unaligned_points, aligned_points, path):
    """``unaligned_points``/``aligned_points``: iterables with ``.cells``-ish and phase fields --
    accepts the Phase 12A ``ConvergencePoint`` (cells_per_side) and Phase 12B
    ``AlignedRefinementPoint`` (cells)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))

    def _n(p):
        return getattr(p, "cells", None) or getattr(p, "cells_per_side", None)

    ax.plot(
        [_n(p) for p in unaligned_points],
        [p.s11_xx_phase_deg for p in unaligned_points],
        marker="o",
        markersize=4,
        linestyle="--",
        label="unaligned (Phase 12A)",
    )
    ax.plot(
        [_n(p) for p in aligned_points],
        [p.s11_xx_phase_deg for p in aligned_points],
        marker="s",
        markersize=4,
        label="aligned (Phase 12B)",
    )
    ax.set_xlabel("cells per side")
    ax.set_ylabel("phase(S11_xx) [deg]")
    ax.set_title("Mesh refinement: aligned vs unaligned (L=W=2.8 mm, 12.5 GHz)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_phase_vs_length_per_period(period_results, path):
    plt = _require_matplotlib()
    periods = sorted(period_results)
    fig, axes = plt.subplots(1, len(periods), figsize=(6.4 * len(periods), 4.4), squeeze=False)
    for ax, period_mm in zip(axes[0], periods):
        points = period_results[period_mm]["points"]
        freqs = sorted({p.frequency_hz for p in points})
        for f in freqs:
            pts = sorted((p for p in points if p.frequency_hz == f), key=lambda p: p.length_width_mm)
            ax.plot(
                [p.length_width_mm for p in pts],
                [p.s11_xx_phase_deg for p in pts],
                marker="o",
                markersize=3,
                label=f"{f / 1e9:.2f} GHz",
            )
        ax.set_xlabel("L = W [mm]")
        ax.set_ylabel("phase(S11_xx) [deg]")
        ax.set_ylim(-180, 180)
        ax.set_title(f"P = {period_mm} mm")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_coverage_vs_frequency_per_period(period_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for period_mm, r in sorted(period_results.items()):
        freqs = sorted(r["coverage"])
        ax.plot(
            [f / 1e9 for f in freqs],
            [r["coverage"][f].coverage_deg for f in freqs],
            marker="o",
            markersize=4,
            label=f"P = {period_mm} mm",
        )
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("S11_xx circular phase coverage [deg]")
    ax.set_ylim(0, 360)
    ax.set_title("Phase coverage vs frequency, by period")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_magnitude_vs_length_per_period(period_results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for period_mm, r in sorted(period_results.items()):
        pts = sorted(
            (p for p in r["points"] if p.frequency_hz == sorted(r["coverage"])[len(r["coverage"]) // 2]),
            key=lambda p: p.length_width_mm,
        )
        ax.plot(
            [p.length_width_mm for p in pts],
            [p.s11_xx_mag for p in pts],
            marker="o",
            markersize=3,
            label=f"P = {period_mm} mm",
        )
    ax.set_xlabel("L = W [mm]")
    ax.set_ylabel("|S11_xx| (linear)")
    ax.set_title("Reflection magnitude (loss) vs size, by period")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_duplicate_audit_table(phase7_report, phase10_groups, path):
    plt = _require_matplotlib()
    rows = [
        (
            "Phase 7 library (324 cells, 32x32 mesh)",
            phase7_report.n_geometries,
            phase7_report.n_distinct,
            phase7_report.n_duplicate_groups,
        ),
        ("Phase 10 centers (5, 8-cell/side mesh)", 5, 5 - sum(len(g) - 1 for g in phase10_groups), len(phase10_groups)),
    ]
    fig, ax = plt.subplots(figsize=(7.5, 1.2 + 0.4 * len(rows)))
    ax.axis("off")
    table = ax.table(
        cellText=[[r[0], str(r[1]), str(r[2]), str(r[3])] for r in rows],
        colLabels=["Set", "n geometries", "n mesh-distinct", "n duplicate groups"],
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    ax.set_title("Duplicate-geometry audit summary", pad=20)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def generate_all_plots(
    unaligned_convergence_points,
    aligned_refinement_points,
    period_results,
    phase7_report,
    phase10_groups,
    plots_dir,
):
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    return {
        "aligned_vs_unaligned_convergence": plot_aligned_vs_unaligned_convergence(
            unaligned_convergence_points,
            aligned_refinement_points,
            plots_dir / "aligned_vs_unaligned_convergence.png",
        ),
        "phase_vs_length_per_period": plot_phase_vs_length_per_period(
            period_results, plots_dir / "phase_vs_length_per_period.png"
        ),
        "coverage_vs_frequency_per_period": plot_coverage_vs_frequency_per_period(
            period_results, plots_dir / "coverage_vs_frequency_per_period.png"
        ),
        "magnitude_vs_length_per_period": plot_magnitude_vs_length_per_period(
            period_results, plots_dir / "magnitude_vs_length_per_period.png"
        ),
        "duplicate_audit_table": plot_duplicate_audit_table(
            phase7_report, phase10_groups, plots_dir / "duplicate_audit_table.png"
        ),
    }
