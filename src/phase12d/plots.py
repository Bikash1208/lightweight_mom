"""Phase 12D's required plots: error vs step size (Task A), per-order power for each N (Task B),
per-cell phase error (Task B), HFSS vs MoM for each available case (Task C)."""


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 12D plots (src/phase12d/plots.py); install it "
            "with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_error_vs_step_size(step_cases, path, design_tol_deg=10.0, design_tol_db=1.0,
                             strict_tol_deg=5.0, strict_tol_db=0.5):
    """Task A: the design-relevant reflection-phase error vs the achieved neighbour phase step,
    one point per step case, with both stated tolerances marked."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    steps = [c.achieved_step_deg for c in step_cases]
    errs = [c.reflection_phase_error_deg for c in step_cases]
    colors = ["tab:blue" if c.kind == "1d_row" else "tab:orange" for c in step_cases]
    ax.scatter(steps, errs, c=colors, s=30, zorder=3)
    ax.axhline(design_tol_deg, color="red", linestyle="--", label=f"design tolerance ({design_tol_deg} deg)")
    ax.axhline(strict_tol_deg, color="darkred", linestyle=":", label=f"strict tolerance ({strict_tol_deg} deg)")
    ax.set_xlabel("achieved neighbour phase step [deg]")
    ax.set_ylabel("centre-cell reflection phase error [deg]")
    ax.set_title("Phase 12D Task A: reflection-phase error vs neighbour step")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_per_order_power(aperture_results, path):
    """Task B: fraction of reflected power per propagating Floquet order, one subplot per tested
    (feasible) N -- shows how power concentrates into the intended order as N grows, or fails to."""
    plt = _require_matplotlib()
    n = len(aperture_results)
    fig, axes = plt.subplots(1, max(n, 1), figsize=(4.5 * max(n, 1), 4.0), squeeze=False)
    for ax, r in zip(axes[0], aperture_results):
        # The only orders worth labelling individually: specular, intended, and the largest spurious ones.
        labels = ["(0,0)", str(r.intended_order)] + [str(mn) for mn, _ in r.spurious_orders[:4]]
        values = [r.fraction_specular, r.fraction_intended] + [f for _, f in r.spurious_orders[:4]]
        ax.bar(range(len(labels)), values, color=["gray", "tab:green"] + ["tab:red"] * 4)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"N={r.n_cells} (span {r.achieved_span_deg:.0f} deg)", fontsize=9)
        ax.set_ylabel("fraction of reflected power")
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Phase 12D Task B: per-order reflected power")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_per_cell_phase_error(aperture_results, path):
    """Task B: per-cell (per-column) reflection-phase error, lookup-table prediction vs the real
    full-supercell MoM solve, one line per tested N."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for r in aperture_results:
        js = [pc.j for pc in r.per_cell_phase_errors]
        errs = [pc.phase_error_deg for pc in r.per_cell_phase_errors]
        ax.plot(js, errs, marker="o", label=f"N={r.n_cells}")
    ax.set_xlabel("cell index (column j)")
    ax.set_ylabel("reflection phase error vs lookup-table prediction [deg]")
    ax.set_title("Phase 12D Task B: per-cell phase-error profile")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_hfss_vs_mom(comparisons, path):
    """Task C: phase deviation (HFSS mode1 - MoM) vs frequency, one line per case with HFSS data."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for c in comparisons:
        ax.plot(c.frequency_ghz, c.dphase_mode1_deg, marker=".", label=f"{c.case_id} (L={c.l_mm} mm)")
    ax.axhline(0.0, color="black", linewidth=0.5)
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("HFSS - MoM S11 phase [deg]")
    ax.set_title("Phase 12D Task C: HFSS vs MoM, single-cell reflection phase")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def generate_all_plots(step_cases, aperture_results, comparisons, output_dir):
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    paths["error_vs_step_size"] = plot_error_vs_step_size(step_cases, output_dir / "error_vs_step_size.png")
    if aperture_results:
        paths["per_order_power"] = plot_per_order_power(aperture_results, output_dir / "per_order_power.png")
        paths["per_cell_phase_error"] = plot_per_cell_phase_error(aperture_results, output_dir / "per_cell_phase_error.png")
    if comparisons:
        paths["hfss_vs_mom"] = plot_hfss_vs_mom(comparisons, output_dir / "hfss_vs_mom.png")
    return paths
