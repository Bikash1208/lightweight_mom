"""Sec. 27: the required plots."""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 11 plots; install with 'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_coupling_vs_distance(directional_points, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    by_dir = {}
    for p in directional_points:
        by_dir.setdefault(p.direction, []).append((p.distance, p.delta_abs))
    for direction, pts in by_dir.items():
        pts = sorted(pts)
        ax.plot([d for d, _ in pts], [v for _, v in pts], "o-", label=direction)
    ax.set_xlabel("distance [cells]")
    ax.set_ylabel("|delta R| [A.m]")
    ax.set_yscale("log")
    ax.set_title("directly computed interaction magnitude vs. distance (no fitted law)")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_shell_increment(increments_by_center, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for name, incs in increments_by_center.items():
        ax.plot([inc.n_to for inc in incs], [inc.delta_relative for inc in incs], "o-", label=name)
    ax.set_xlabel("grid size n (shell added going TO this n)")
    ax.set_ylabel("relative shell increment")
    ax.set_yscale("log")
    ax.set_title("incremental per-shell contribution")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_cardinal_vs_diagonal(cardinal_mean, diagonal_mean, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.bar(["cardinal", "diagonal"], [cardinal_mean, diagonal_mean])
    ax.set_ylabel("mean |delta R| [A.m]")
    ax.set_title(
        f"cardinal vs. diagonal (ratio={diagonal_mean / cardinal_mean:.2f})"
        if cardinal_mean > 0
        else "cardinal vs. diagonal"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_directional_map(directional_points, distance, path):
    plt = _require_matplotlib()
    angle_map = {"E": 0, "NE": 45, "N": 90, "NW": 135, "W": 180, "SW": 225, "S": 270, "SE": 315}
    pts = [p for p in directional_points if p.distance == distance]
    fig, ax = plt.subplots(figsize=(5.5, 5.5), subplot_kw={"projection": "polar"})
    thetas = [np.radians(angle_map[p.direction]) for p in pts]
    radii = [p.delta_abs for p in pts]
    ax.scatter(thetas, radii, s=60)
    for t, r, p in zip(thetas, radii, pts):
        ax.annotate(p.direction, (t, r), fontsize=8)
    ax.set_title(f"directional coupling map at distance={distance}")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_grid_geometry(grid_mm, path, title):
    plt = _require_matplotlib()
    n = len(grid_mm)
    l_vals = np.array([[grid_mm[i][j][0] for j in range(n)] for i in range(n)])
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(l_vals, cmap="viridis")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="L [mm]")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_smooth_gradient(grid_mm, path):
    _plot_grid_geometry(grid_mm, path, "smooth gradient pattern")


def plot_abrupt_boundary(grid_mm, path):
    _plot_grid_geometry(grid_mm, path, "abrupt boundary pattern")


def plot_alternating_pattern(grid_mm, path):
    _plot_grid_geometry(grid_mm, path, "alternating pattern")


def plot_checkerboard_pattern(grid_mm, path):
    _plot_grid_geometry(grid_mm, path, "checkerboard pattern")


def plot_frequency_dependence(frequency_points, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    freqs = [p.frequency_hz / 1e9 for p in frequency_points]
    max_incs = [p.max_shell_increment_relative for p in frequency_points]
    ax.plot(freqs, max_incs, "o-")
    ax.set_xlabel("frequency [GHz]")
    ax.set_ylabel("max shell-increment relative error")
    ax.set_title("coupling extent vs. frequency")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_convergence_vs_neighborhood(points, noise_floor, substantial, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ns = [p.n for p in points]
    mags = [abs(p.response) for p in points]
    ax.plot(ns, mags, "o-")
    ax.set_xlabel("grid size n")
    ax.set_ylabel("|center response| [A.m]")
    ax.set_title("center response vs. simulated extent")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_center_current_distribution(sg, n, path, n_samples=64):
    plt = _require_matplotlib()
    import numpy as np

    period = sg.supercell.sub_period_x * n
    lin = np.linspace(-period / 2, period / 2, n_samples)
    x, y = np.meshgrid(lin, lin)
    j = sg.basis.evaluate(x, y, sg.coeffs, 0.0, 0.0)
    mag = np.abs(j[..., 0])
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(mag, extent=[-period / 2 * 1e3, period / 2 * 1e3, -period / 2 * 1e3, period / 2 * 1e3])
    ax.set_title("surface current magnitude |Jx|")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_heterogeneous_far_field(far_field, theta_grid_deg, phi_grid_deg, path, title="heterogeneous far field"):
    plt = _require_matplotlib()
    power_db = 10 * np.log10(np.maximum(np.abs(far_field) ** 2, 1e-30))
    power_db -= np.max(power_db)
    fig, ax = plt.subplots(figsize=(6.5, 5), subplot_kw={"projection": "polar"})
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
