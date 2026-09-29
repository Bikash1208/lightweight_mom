"""Phase 7 plots (req. 21), all built from already-computed library data (a `database.Database`,
`phase_extraction` rows/coverage/table, `library_classifier` results) -- no plot here fabricates
or resamples a curve; every plotted point is a value already stored by the solver or a documented
deterministic transform of it (dB via `validation_export._db`, the same floor/definition Phase 6
and 6A use, req. 23's "do not create a second incompatible convention").
"""

from pathlib import Path

import numpy as np

from .library_schema import RESONANCE_MIN_PROMINENCE
from .validation_export import _db


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 7 library plots (src/library_plots.py); "
            "install it with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def select_representative_cells(config, n_per_axis=3):
    """Deterministic representative-cell selection for the 4 per-frequency plots (req. 21.A-D):
    a `n_per_axis` x `n_per_axis` sub-grid of L/W index positions spread evenly across the sweep
    (`round(np.linspace(0, n-1, n_per_axis))`), i.e. for the locked 18-value default and
    n_per_axis=3: indices [0, 8 or 9, 17] -> L, W in {1.0, ~2.6, 4.4} mm, corners and centre of
    the L x W square -- 9 cells total, not all 324. Every plot built from this selection states
    exactly this rule and lists the (L, W) pairs used, so "representative" is reproducible, not
    a hand pick (req. 21.A's "document exactly which cells were selected")."""
    n_l, n_w = len(config.l_values_mm), len(config.w_values_mm)
    li = sorted(set(int(round(x)) for x in np.linspace(0, n_l - 1, n_per_axis)))
    wi = sorted(set(int(round(x)) for x in np.linspace(0, n_w - 1, n_per_axis)))
    return [(config.l_values_mm[i], config.w_values_mm[j]) for i in li for j in wi]


def _cell_for(database, l_mm, w_mm):
    for c in database.cells:
        if abs(c.parameters["length_m"] - l_mm * 1e-3) < 1e-15 and abs(c.parameters["width_m"] - w_mm * 1e-3) < 1e-15:
            return c
    raise ValueError(f"no cell for L={l_mm} mm, W={w_mm} mm")


def plot_s_parameter_curve(database, cells, family, quantity, path):
    """Req. 21.A-D: one figure, one curve per representative cell (`select_representative_cells`,
    documented and deterministic). ``family``: 'S11' or 'S21'; ``quantity``: 'mag' (dB) or
    'phase' (wrapped degrees, no unwrapping) of the co-polarized (`_xx`) component."""
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for l_mm, w_mm in cells:
        cell = _cell_for(database, l_mm, w_mm)
        if cell.status != "ok":
            continue
        records = sorted(database.records_of(cell.index), key=lambda r: r.frequency_index)
        freq_ghz = np.array([r.frequency for r in records]) / 1e9
        label = f"L={l_mm:g}, W={w_mm:g} mm"
        if quantity == "mag":
            mag = np.array([r.magnitude(f"{family}_xx") for r in records])
            ax.plot(freq_ghz, _db(mag), marker="o", markersize=2, label=label)
            ax.set_ylabel(f"|{family}xx| [dB]")
        else:
            phase = np.array([r.phase_deg(f"{family}_xx") for r in records])
            ax.plot(freq_ghz, phase, marker="o", markersize=2, label=label)
            ax.set_ylabel(f"phase({family}xx) [deg] (wrapped)")
            ax.set_ylim(-180, 180)
    ax.set_xlabel("frequency [GHz]")
    ax.set_title(f"{family}xx {quantity}, representative cells ({', '.join(f'L={l:g}/W={w:g}' for l, w in cells)})")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _grid_from_table(table, config, value_key):
    """(L_axis, W_axis, values[len(W), len(L)]) for pcolormesh, NaN where a geometry is absent
    from ``table`` (a failed/invalid cell)."""
    l_axis = np.array(config.l_values_mm)
    w_axis = np.array(config.w_values_mm)
    grid = np.full((len(w_axis), len(l_axis)), np.nan)
    l_index = {v: i for i, v in enumerate(l_axis)}
    w_index = {v: i for i, v in enumerate(w_axis)}
    for row in table:
        grid[w_index[row["W_mm"]], l_index[row["L_mm"]]] = row[value_key]
    return l_axis, w_axis, grid


def plot_phase_vs_geometry(table, config, path, design_frequency_ghz):
    """Req. 21.E: the design-frequency S11_xx phase over the full L/W geometry space, as a
    heatmap on the actual 18x18 sweep grid (every cell plotted, not a representative subset)."""
    plt = _require_matplotlib()
    l_axis, w_axis, phase_grid = _grid_from_table(table, config, "phase_deg")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    mesh = ax.pcolormesh(l_axis, w_axis, phase_grid, shading="nearest", cmap="twilight", vmin=-180, vmax=180)
    fig.colorbar(mesh, ax=ax, label="phase(S11xx) [deg]")
    ax.set_xlabel("L [mm]")
    ax.set_ylabel("W [mm]")
    ax.set_title(f"S11xx phase vs geometry at f_design = {design_frequency_ghz:g} GHz")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_phase_coverage_map(table, coverage, path, design_frequency_ghz):
    """Req. 21.F: every cell's design-frequency phase placed on the unit circle, with the single
    largest empty arc (`phase_extraction.phase_coverage`) shaded -- a direct visualization of
    what is and is not covered, not a claim of "full coverage" (req. 12's caution)."""
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
    fig.suptitle(f"Phase coverage at f_design = {design_frequency_ghz:g} GHz", y=0.97)
    fig.text(
        0.5,
        0.08,
        f"{coverage.n_points} cells  |  coverage = {coverage.coverage_deg:.1f} deg of 360 deg  |  "
        f"largest gap = {coverage.largest_gap_deg:.1f} deg",
        ha="center",
        fontsize=9,
    )
    fig.legend(handles=[pts, fill], loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=1, fontsize=8)
    fig.tight_layout(rect=(0, 0.16, 1, 0.93))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_resonance_map(classifications, config, path):
    """Req. 21.G: which (L, W) cells were classified RESONANT vs NON_RESONANT (req. 14's
    documented prominence criterion, `library_schema.RESONANCE_MIN_PROMINENCE`)."""
    plt = _require_matplotlib()
    l_axis = np.array(config.l_values_mm)
    w_axis = np.array(config.w_values_mm)
    grid = np.full((len(w_axis), len(l_axis)), np.nan)
    l_index = {v: i for i, v in enumerate(l_axis)}
    w_index = {v: i for i, v in enumerate(w_axis)}
    for c in classifications:
        value = np.nan if c.resonance is None else float(c.resonance.is_resonant)
        grid[w_index[c.w_mm], l_index[c.l_mm]] = value
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    mesh = ax.pcolormesh(l_axis, w_axis, grid, shading="nearest", cmap="RdYlGn_r", vmin=0, vmax=1)
    cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["NON_RESONANT", "RESONANT"])
    ax.set_xlabel("L [mm]")
    ax.set_ylabel("W [mm]")
    ax.set_title(f"Resonance classification (|S11xx| local-max prominence >= {RESONANCE_MIN_PROMINENCE:g})")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


_QUALITY_CODE = {
    ("INVALID", "FAILED"): 0,
    ("VALID", "FAILED"): 1,
    ("VALID", "OK"): 2,  # further split by phase_status via colour below is skipped for clarity
}
_QUALITY_LABELS = ["INVALID/FAILED", "SOLVE FAILED", "OK (see phase_status separately)"]


def plot_library_quality_map(classifications, config, path):
    """Req. 21.H: valid/invalid/unclassified cells over the L/W grid, from
    `library_classifier.classify_cell`'s independent geometry_status/solver_status flags."""
    plt = _require_matplotlib()
    l_axis = np.array(config.l_values_mm)
    w_axis = np.array(config.w_values_mm)
    grid = np.full((len(w_axis), len(l_axis)), np.nan)
    l_index = {v: i for i, v in enumerate(l_axis)}
    w_index = {v: i for i, v in enumerate(w_axis)}
    for c in classifications:
        code = _QUALITY_CODE.get((c.geometry_status, c.solver_status), np.nan)
        grid[w_index[c.w_mm], l_index[c.l_mm]] = code
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    mesh = ax.pcolormesh(l_axis, w_axis, grid, shading="nearest", cmap="RdYlGn", vmin=0, vmax=2)
    cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2])
    cbar.ax.set_yticklabels(_QUALITY_LABELS)
    ax.set_xlabel("L [mm]")
    ax.set_ylabel("W [mm]")
    ax.set_title("Library quality map (geometry_status / solver_status)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def generate_all_plots(database, config, table, coverage, classifications, plots_dir):
    """Write all 8 required plots (req. 21) into ``plots_dir``; returns {name: path}."""
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    cells = select_representative_cells(config)
    design_ghz = config.f_design_hz / 1e9

    paths = {
        "s11_magnitude_db": plots_dir / "s11_magnitude_db.png",
        "s11_phase_deg": plots_dir / "s11_phase_deg.png",
        "s21_magnitude_db": plots_dir / "s21_magnitude_db.png",
        "s21_phase_deg": plots_dir / "s21_phase_deg.png",
        "phase_vs_geometry": plots_dir / "phase_vs_geometry.png",
        "phase_coverage_map": plots_dir / "phase_coverage_map.png",
        "resonance_map": plots_dir / "resonance_map.png",
        "library_quality_map": plots_dir / "library_quality_map.png",
    }
    plot_s_parameter_curve(database, cells, "S11", "mag", paths["s11_magnitude_db"])
    plot_s_parameter_curve(database, cells, "S11", "phase", paths["s11_phase_deg"])
    plot_s_parameter_curve(database, cells, "S21", "mag", paths["s21_magnitude_db"])
    plot_s_parameter_curve(database, cells, "S21", "phase", paths["s21_phase_deg"])

    plot_phase_vs_geometry(table, config, paths["phase_vs_geometry"], design_ghz)
    plot_phase_coverage_map(table, coverage, paths["phase_coverage_map"], design_ghz)
    plot_resonance_map(classifications, config, paths["resonance_map"])
    plot_library_quality_map(classifications, config, paths["library_quality_map"])
    return paths
