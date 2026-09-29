"""Phase 8C Sec. 46: the required plots. Only a plot whose underlying quantity actually exists is
generated (Sec. 46's "only generate a plot if the corresponding quantity actually exists")."""

from pathlib import Path

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 8C plots (src/phase8c/plots.py); install with "
            "'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def _grid_from_rows(rows, l_values_mm, w_values_mm, value_fn):
    l_axis, w_axis = np.array(sorted(l_values_mm)), np.array(sorted(w_values_mm))
    grid = np.full((len(w_axis), len(l_axis)), np.nan)
    l_index = {round(v, 9): i for i, v in enumerate(l_axis)}
    w_index = {round(v, 9): i for i, v in enumerate(w_axis)}
    for r in rows:
        li, wi = l_index.get(round(r.l_mm, 9)), w_index.get(round(r.w_mm, 9))
        if li is not None and wi is not None:
            grid[wi, li] = value_fn(r)
    return l_axis, w_axis, grid


def plot_topology_summary(audit, path):
    """01_unit_cell_topology.png: a text-only summary figure (this project has no 3-D rendering
    dependency chain here; the STEP files themselves are the geometry artifact, Sec. 58)."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.axis("off")
    text = (
        f"Ground plane: {audit.ground_plane}\n\n"
        f"Metal: {audit.metal_layers}\n\n"
        f"Substrate: {audit.substrate_layers}\n\n"
        f"Architecture: {audit.architecture_classification}\n({audit.architecture_evidence})"
    )
    ax.text(0.02, 0.98, text, va="top", fontsize=10, family="monospace", wrap=True)
    ax.set_title("Unit-cell topology summary")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_representative_curve(database, cells, component_name, quantity, path):
    """02-05: representative-cell |S| or phase(S) vs frequency."""
    from ..validation_export import _db

    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for l_mm, w_mm in cells:
        cell = next(
            (c for c in database.cells if abs(c.parameters["length_m"] - l_mm * 1e-3) < 1e-15
             and abs(c.parameters["width_m"] - w_mm * 1e-3) < 1e-15), None,
        )
        if cell is None or cell.status != "ok":
            continue
        records = sorted(database.records_of(cell.index), key=lambda r: r.frequency_index)
        freq_ghz = np.array([r.frequency for r in records]) / 1e9
        label = f"L={l_mm:g}, W={w_mm:g} mm"
        if quantity == "mag":
            mag = np.array([r.magnitude(component_name) for r in records])
            ax.plot(freq_ghz, _db(mag), marker="o", markersize=2, label=label)
            ax.set_ylabel(f"|{component_name}| [dB]")
        else:
            phase = np.array([r.phase_deg(component_name) for r in records])
            ax.plot(freq_ghz, phase, marker="o", markersize=2, label=label)
            ax.set_ylabel(f"phase({component_name}) [deg] (wrapped)")
            ax.set_ylim(-180, 180)
    ax.set_xlabel("frequency [GHz]")
    ax.set_title(f"{component_name} {quantity}, representative cells")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_library_map(rows, l_values_mm, w_values_mm, value_fn, label, title, path, vmin=None, vmax=None, cmap="viridis"):
    """06-09: library-wide (L,W) map for any per-row value."""
    plt = _require_matplotlib()
    l_axis, w_axis, grid = _grid_from_rows(rows, l_values_mm, w_values_mm, value_fn)
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    mesh = ax.pcolormesh(l_axis, w_axis, grid, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(mesh, ax=ax, label=label)
    ax.set_xlabel("L [mm]"); ax.set_ylabel("W [mm]"); ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_complex_plane(rows, component_name, path):
    """10-11: complex-plane locus of the candidate response over the full library."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5.6, 5.6))
    re = [r.real for r in rows]; im = [r.imag for r in rows]
    ax.scatter(re, im, s=10)
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(theta), np.sin(theta), color="grey", linewidth=0.7, linestyle="--")
    ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.set_aspect("equal")
    ax.set_title(f"{component_name} complex-plane locus ({len(rows)} cells)")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_phase_vs_magnitude(rows, component_name, path):
    """12-13: phase vs magnitude scatter."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.scatter([r.phase_deg for r in rows], [r.mag for r in rows], s=10)
    ax.set_xlabel(f"phase({component_name}) [deg]"); ax.set_ylabel(f"|{component_name}|")
    ax.set_title(f"{component_name}: phase vs magnitude")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_phase_coverage(rows, coverage, component_name, path):
    """14-15: circular phase coverage (same style as Phase 7/8A/8B's coverage plot)."""
    plt = _require_matplotlib()
    phases_rad = np.radians([r.phase_deg for r in rows])
    fig, ax = plt.subplots(figsize=(6.4, 7.2), subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1)
    ax.scatter(phases_rad, np.ones_like(phases_rad), s=12, zorder=3)
    fig.suptitle(f"{component_name} phase coverage: {coverage.coverage_deg:.1f} deg of 360 deg", y=0.97)
    ax.set_yticklabels([])
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_power_balance(power_summary, path):
    """16_power_balance.png."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    labels = ["R (co)", "R (cross)", "T (co)", "T (cross)", "1-(R+T)"]
    values = [power_summary["mean_R_co"], power_summary["mean_R_cross"], power_summary["mean_T_co"],
              power_summary["mean_T_cross"], 1.0 - power_summary["mean_R_plus_T"]]
    ax.bar(labels, values)
    ax.set_ylabel("mean power fraction"); ax.set_title("Library-mean power balance (x-incidence)")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_phase_state_distribution(analysis, path, component_name):
    """17-18: distinct-complex-vs-phase-state bar chart."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.bar(["geometries", "complex states", "phase states"],
           [analysis.n_geometries, analysis.n_distinct_complex_states, analysis.n_distinct_phase_states])
    ax.set_title(f"{component_name}: geometric vs. complex vs. phase diversity")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
