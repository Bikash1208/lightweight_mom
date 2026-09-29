"""Phase 12A's 4 required plots (deliverables list): phase vs L at several frequencies,
|S11| vs L, grounded vs ungrounded comparison, convergence. matplotlib only (optional dependency,
same rule as ``cad_export.py``/``validation_export.py``); never imported from the solver core.
"""

from pathlib import Path


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 12A plots (src/phase12a/plots.py); install it "
            "with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_phase_vs_length(points, path):
    plt = _require_matplotlib()
    freqs = sorted({p.frequency_hz for p in points})
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for f in freqs:
        pts = sorted((p for p in points if p.frequency_hz == f), key=lambda p: p.length_width_mm)
        ax.plot(
            [p.length_width_mm for p in pts], [p.s11_xx_phase_deg for p in pts],
            marker="o", markersize=3, label=f"{f / 1e9:.2f} GHz",
        )
    ax.set_xlabel("L = W [mm]")
    ax.set_ylabel("phase(S11_xx) [deg]")
    ax.set_ylim(-180, 180)
    ax.set_title("Grounded rectangular-PEC-patch reflection phase vs size")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_magnitude_vs_length(points, path):
    plt = _require_matplotlib()
    freqs = sorted({p.frequency_hz for p in points})
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for f in freqs:
        pts = sorted((p for p in points if p.frequency_hz == f), key=lambda p: p.length_width_mm)
        ax.plot(
            [p.length_width_mm for p in pts], [p.s11_xx_mag for p in pts],
            marker="o", markersize=3, label=f"{f / 1e9:.2f} GHz",
        )
    ax.set_xlabel("L = W [mm]")
    ax.set_ylabel("|S11_xx| (linear)")
    ax.set_title("Grounded rectangular-PEC-patch reflection magnitude (loss) vs size")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_grounded_vs_ungrounded(grounded_points, ungrounded_table, design_frequency_hz, path):
    """``ungrounded_table``: list of {"L_mm", "phase_deg"} at ``design_frequency_hz`` (Phase 7's
    own library, read via ``ungrounded_phase7_coverage``'s sibling data if available; the plot
    just shows what is actually there, no fabrication if the Phase 7 checkpoint is missing)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    g = sorted(
        (p for p in grounded_points if p.frequency_hz == design_frequency_hz),
        key=lambda p: p.length_width_mm,
    )
    ax.plot([p.length_width_mm for p in g], [p.s11_xx_phase_deg for p in g],
            marker="o", markersize=3, label="grounded (backing=pec)")
    if ungrounded_table:
        u = sorted(ungrounded_table, key=lambda r: r["L_mm"])
        ax.plot([r["L_mm"] for r in u], [r["phase_deg"] for r in u],
                marker="s", markersize=3, linestyle="--", label="ungrounded (Phase 7 library)")
    else:
        ax.text(0.5, 0.5, "ungrounded Phase 7 data NOT_AVAILABLE", transform=ax.transAxes,
                ha="center", va="center", color="red")
    ax.set_xlabel("L = W [mm]")
    ax.set_ylabel("phase(S11_xx) [deg]")
    ax.set_ylim(-180, 180)
    ax.set_title(f"Grounded vs ungrounded reflection phase at {design_frequency_hz / 1e9:.2f} GHz")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_convergence(points, path):
    plt = _require_matplotlib()
    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
    n = [p.cells_per_side for p in points]
    ax_mag.plot(n, [p.s11_xx_mag for p in points], marker="o")
    ax_mag.set_ylabel("|S11_xx|")
    ax_mag.set_title("Mesh/Floquet convergence, one grounded patch (floquet_order = 2 x cells_per_side)")
    ax_mag.grid(True, alpha=0.3)
    ax_phase.plot(n, [p.s11_xx_phase_deg for p in points], marker="o")
    ax_phase.set_ylabel("phase(S11_xx) [deg]")
    ax_phase.set_xlabel("cells per side")
    ax_phase.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def generate_all_plots(sweep_points, ungrounded_table, design_frequency_hz, convergence_points, plots_dir):
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    return {
        "phase_vs_length": plot_phase_vs_length(sweep_points, plots_dir / "phase_vs_length.png"),
        "magnitude_vs_length": plot_magnitude_vs_length(sweep_points, plots_dir / "magnitude_vs_length.png"),
        "grounded_vs_ungrounded": plot_grounded_vs_ungrounded(
            sweep_points, ungrounded_table, design_frequency_hz, plots_dir / "grounded_vs_ungrounded.png"
        ),
        "convergence": plot_convergence(convergence_points, plots_dir / "convergence.png"),
    }
