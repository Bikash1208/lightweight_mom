"""Phase 12C's required plots: deviation vs phase step, embedded vs isolated per cell for the
phase-gradient supercell, coupling vs size (grounded vs ungrounded), fine-vs-coarse spot-check
table.
"""

from pathlib import Path


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 12C plots (src/phase12c/plots.py); install it "
            "with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_deviation_vs_phase_step(gradient_results, path):
    """Sec. B: relative complex LPA error vs the achieved phase step, one point per gradient case,
    with Phase 9's own noise-floor/substantial-deviation reference lines."""
    from ..phase9.local_periodicity import NOISE_FLOOR_RELATIVE, SUBSTANTIAL_DEVIATION_RELATIVE

    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    steps = [abs(r.achieved_step_deg) for r in gradient_results]
    errs = [r.relative_complex_error for r in gradient_results]
    ax.scatter(steps, errs, s=30, zorder=3)
    for r in gradient_results:
        ax.annotate(r.name, (abs(r.achieved_step_deg), r.relative_complex_error), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.axhline(NOISE_FLOOR_RELATIVE, color="green", linestyle=":", label=f"noise floor ({NOISE_FLOOR_RELATIVE:.0e})")
    ax.axhline(SUBSTANTIAL_DEVIATION_RELATIVE, color="red", linestyle="--",
               label=f"substantial ({SUBSTANTIAL_DEVIATION_RELATIVE})")
    ax.set_yscale("log")
    ax.set_xlabel("achieved neighbour phase step [deg]")
    ax.set_ylabel("relative complex LPA error")
    ax.set_title("Grounded 8 mm cell: local-periodicity deviation vs phase step")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_periodic_supercell_embedded_vs_isolated(periodic_results, path):
    """Sec. B.3: |response| and phase, embedded vs isolated, for every column of the periodic
    phase-gradient supercell."""
    import cmath

    plt = _require_matplotlib()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    js = [r.j for r in periodic_results]
    iso_mag = [abs(r.isolated_px) for r in periodic_results]
    emb_mag = [abs(r.embedded_px) for r in periodic_results]
    iso_phase = [cmath.phase(r.isolated_px) * 180 / 3.141592653589793 for r in periodic_results]
    emb_phase = [cmath.phase(r.embedded_px) * 180 / 3.141592653589793 for r in periodic_results]
    ax1.plot(js, iso_mag, marker="o", label="isolated")
    ax1.plot(js, emb_mag, marker="s", label="embedded")
    ax1.set_xlabel("column j")
    ax1.set_ylabel("|dipole moment|")
    ax1.set_title("Magnitude")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax2.plot(js, iso_phase, marker="o", label="isolated")
    ax2.plot(js, emb_phase, marker="s", label="embedded")
    ax2.set_xlabel("column j")
    ax2.set_ylabel("phase [deg]")
    ax2.set_title("Phase")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    fig.suptitle("Periodic phase-gradient supercell: embedded vs isolated, per column")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_coupling_vs_size_grounded_vs_ungrounded(grounded_reports, ungrounded_series, path):
    """Sec. C: |response| shell-increment (relative) vs n, grounded (this phase) vs ungrounded
    (Phase 11's own already-measured series, read-only)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for name, report in grounded_reports.items():
        ns = [inc.n_to for inc in report.increments]
        rels = [inc.delta_relative for inc in report.increments]
        ax.plot(ns, rels, marker="o", label=f"grounded: {name}")
    if ungrounded_series:
        for name, points in ungrounded_series.items():
            from ..phase11.coupling_radius import compute_shell_increments

            incs = compute_shell_increments(points)
            ax.plot([i.n_to for i in incs], [i.delta_relative for i in incs], marker="x",
                    linestyle="--", label=f"ungrounded: {name}")
    ax.set_yscale("log")
    ax.set_xlabel("supercell size n")
    ax.set_ylabel("relative shell increment")
    ax.set_title("Coupling extent: grounded vs ungrounded")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def plot_fine_vs_coarse_spot_check_table(spot_gate, path):
    """Sec. D: a table image (Phase 12B's own ``duplicate_audit_table.png`` precedent) of every
    spot-check case's verdict."""
    plt = _require_matplotlib()
    rows = [[name, verdict, detail[:70] + ("..." if len(detail) > 70 else "")]
            for name, verdict, detail in spot_gate.findings]
    fig, ax = plt.subplots(figsize=(11, 0.6 * (len(rows) + 1) + 0.5))
    ax.axis("off")
    table = ax.table(cellText=rows, colLabels=["case", "verdict", "detail"], loc="center",
                      cellLoc="left", colWidths=[0.2, 0.3, 0.5])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.6)
    ax.set_title("Phase 12C Task D: fine-vs-coarse mesh spot-check", pad=20)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return Path(path)


def generate_all_plots(gradient_results, periodic_results, grounded_coupling_reports,
                        ungrounded_coupling_series, spot_gate, plots_dir):
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    plot_deviation_vs_phase_step(gradient_results, plots_dir / "deviation_vs_phase_step.png")
    plot_periodic_supercell_embedded_vs_isolated(periodic_results, plots_dir / "periodic_supercell_embedded_vs_isolated.png")
    plot_coupling_vs_size_grounded_vs_ungrounded(grounded_coupling_reports, ungrounded_coupling_series,
                                                   plots_dir / "coupling_vs_size_grounded_vs_ungrounded.png")
    plot_fine_vs_coarse_spot_check_table(spot_gate, plots_dir / "fine_vs_coarse_spot_check_table.png")
