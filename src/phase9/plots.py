"""Sec. 26: the required plots. Where HFSS data does not exist (every production run of this
project, as of Phase 9), the corresponding panel is explicitly annotated `HFSS data NOT AVAILABLE`
rather than left blank or silently omitted.
"""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError("matplotlib is required for the Phase 9 plots; install with "
                           "'pip install -r requirements-cad.txt'.") from e
    return plt


def _rows_for_cell(rows, cell_id):
    return sorted([r for r in rows if r.cell_id == cell_id], key=lambda r: r.frequency_ghz)


def _annotate_no_hfss(ax):
    ax.text(0.5, 0.5, "HFSS data NOT AVAILABLE", ha="center", va="center", transform=ax.transAxes,
            fontsize=11, color="crimson", style="italic")


# --------------------------------------------------------------------------------------- 9A: 1-8

def plot_s_magnitude(rows, hfss_records, cell_id, component, path):
    plt = _require_matplotlib()
    cell_rows = _rows_for_cell(rows, cell_id)
    freqs = [r.frequency_ghz for r in cell_rows]
    solver_mag = [abs(getattr(r, component)) for r in cell_rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(freqs, solver_mag, "o-", label="solver")
    if hfss_records:
        hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records}
        hfss_mag = [abs(getattr(hfss_by_key[(cell_id, round(f, 6))], component)) for f in freqs
                    if (cell_id, round(f, 6)) in hfss_by_key]
        if hfss_mag:
            ax.plot(freqs[: len(hfss_mag)], hfss_mag, "s--", label="HFSS")
    else:
        _annotate_no_hfss(ax)
    ax.set_xlabel("frequency [GHz]"); ax.set_ylabel(f"|{component}|")
    ax.set_title(f"{component} magnitude, {cell_id}")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_s_phase(rows, hfss_records, cell_id, component, path):
    plt = _require_matplotlib()
    cell_rows = _rows_for_cell(rows, cell_id)
    freqs = [r.frequency_ghz for r in cell_rows]
    solver_phase = [np.degrees(np.angle(getattr(r, component))) for r in cell_rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(freqs, solver_phase, "o-", label="solver")
    if not hfss_records:
        _annotate_no_hfss(ax)
    ax.set_xlabel("frequency [GHz]"); ax.set_ylabel(f"phase({component}) [deg]")
    ax.set_title(f"{component} phase, {cell_id}")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_s_complex_plane(rows, hfss_records, component, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    values = [getattr(r, component) for r in rows]
    ax.scatter([v.real for v in values], [v.imag for v in values], s=14, label="solver")
    if not hfss_records:
        _annotate_no_hfss(ax)
    ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_title(f"{component} complex-plane distribution (reference cells)")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_s_error_vs_frequency(rows, hfss_records, cell_id, component, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    if not hfss_records:
        _annotate_no_hfss(ax)
    else:
        from .error_metrics import complex_error as complex_error_fn
        hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records}
        cell_rows = _rows_for_cell(rows, cell_id)
        freqs, errs = [], []
        for r in cell_rows:
            key = (cell_id, round(r.frequency_ghz, 6))
            if key in hfss_by_key:
                freqs.append(r.frequency_ghz)
                errs.append(float(complex_error_fn(getattr(r, component), getattr(hfss_by_key[key], component))))
        ax.plot(freqs, errs, "o-")
    ax.set_xlabel("frequency [GHz]"); ax.set_ylabel(f"|{component} error|")
    ax.set_title(f"{component} complex error vs. frequency, {cell_id}")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_error_summary_across_geometries(rows, hfss_records, component, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if not hfss_records:
        _annotate_no_hfss(ax)
    else:
        from .error_metrics import complex_error as complex_error_fn
        hfss_by_key = {(r.cell_id, round(r.frequency_ghz, 6)): r for r in hfss_records}
        by_cell = {}
        for r in rows:
            key = (r.cell_id, round(r.frequency_ghz, 6))
            if key in hfss_by_key:
                by_cell.setdefault(r.cell_id, []).append(
                    float(complex_error_fn(getattr(r, component), getattr(hfss_by_key[key], component)))
                )
        cells = sorted(by_cell)
        ax.bar(cells, [np.mean(by_cell[c]) for c in cells])
        ax.tick_params(axis="x", rotation=45)
    ax.set_ylabel(f"mean |{component} error|")
    ax.set_title(f"{component} error summary across reference geometries")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


# --------------------------------------------------------------------------------------- 9B: 9-16

def plot_deviation_bar(results, names, value_fn, ylabel, title, path):
    plt = _require_matplotlib()
    subset = [r for r in results if r.name in names]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar([r.name for r in subset], [value_fn(r) for r in subset])
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_contrast_vs_error(results, value_fn, ylabel, title, path):
    plt = _require_matplotlib()
    order = ["3x3_weak", "3x3_moderate", "3x3_symmetric", "3x3_strong"]
    subset = [next(r for r in results if r.name == name) for name in order if any(r.name == name for r in results)]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot([s.name for s in subset], [value_fn(s) for s in subset], "o-")
    ax.set_ylabel(ylabel); ax.set_title(title)
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_locality_convergence(results, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    by_name = {r.name: r for r in results}
    for label, n3, n5 in (("moderate", "3x3_moderate", "5x5_moderate"), ("strong", "3x3_strong", "5x5_strong")):
        if n3 in by_name and n5 in by_name:
            ax.plot([3, 5], [by_name[n3].relative_complex_error, by_name[n5].relative_complex_error],
                    "o-", label=label)
    ax.set_xlabel("supercell size n (n x n)"); ax.set_ylabel("relative complex error")
    ax.set_xticks([3, 5])
    ax.set_title("3x3 vs 5x5 convergence (locality check)")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
