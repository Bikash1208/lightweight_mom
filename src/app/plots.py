"""PHASE13C: server-rendered PNG plots (matplotlib -- already an `app` extra dependency, see
`pyproject.toml`'s own comment "plots are rendered server-side with matplotlib"). No client-side
plotting library or build step; every plot is a small, self-contained PNG a job writes to its own
`job_output_dir` and serves through the existing `/api/files/{job_id}/{filename}` download route --
the same mechanism STEP export already uses.

Styled to (loosely) match the app's light/Apple-like theme (`style.css`'s `--accent`): white
figure background, one accent color, readable at ~800px card width. This is presentation only --
never a second source of truth for a number also shown as text.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

_ACCENT = "#0071e3"  # system blue -- matches style.css --accent
_FIGURE_DPI = 150


def _new_fig(figsize=(5.5, 4.5)):
    fig, ax = plt.subplots(figsize=figsize, dpi=_FIGURE_DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.grid(True, linewidth=0.4, color="#d2d2d7", alpha=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#d2d2d7")
    return fig, ax


def _save_png(fig, path):
    fig.savefig(path, format="png", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def reference_phase_map_png(x_m, y_m, phase_deg, path, title="Reference phase"):
    """Step 1's "a PNG of the reference phase map"."""
    fig, ax = _new_fig()
    im = ax.pcolormesh(np.asarray(x_m) * 1e3, np.asarray(y_m) * 1e3, phase_deg, shading="auto", cmap="twilight")
    fig.colorbar(im, ax=ax, label="phase (deg)")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)
    ax.set_aspect("equal")
    return _save_png(fig, path)


def layout_map_png(x_m, y_m, value, path, title="Layout", label="value", cmap="cividis"):
    """A per-cell scalar map (patch length, required/achieved phase, phase error, ...) -- reused
    for the hologram/reflectarray layout preview and the error-map plots."""
    fig, ax = _new_fig()
    im = ax.pcolormesh(np.asarray(x_m) * 1e3, np.asarray(y_m) * 1e3, value, shading="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=label)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)
    ax.set_aspect("equal")
    return _save_png(fig, path)


def pattern_uv_map_png(theta_deg, phi_deg, power_db, path, beams_theta_phi_deg=(), title="Radiation pattern"):
    """2-D u-v map in dB (13C Step 3), u=sin(theta)cos(phi), v=sin(theta)sin(phi) -- the standard
    direction-cosine projection, undistorted near broadside unlike a raw theta-phi plot. A
    theta/phi grid maps to a non-rectangular point cloud in (u,v), so this triangulates
    (`tripcolor`) rather than treating u/v as a regular mesh (`pcolormesh` would silently produce
    the wrong cell edges here)."""
    theta_rad = np.radians(theta_deg).ravel()
    phi_rad = np.radians(phi_deg).ravel()
    u = np.sin(theta_rad) * np.cos(phi_rad)
    v = np.sin(theta_rad) * np.sin(phi_rad)
    power_db = np.asarray(power_db).ravel()
    peak_db = float(np.nanmax(power_db))
    fig, ax = _new_fig()
    im = ax.tripcolor(u, v, power_db, shading="gouraud", cmap="viridis", vmin=peak_db - 40.0, vmax=peak_db)
    fig.colorbar(im, ax=ax, label="relative power (dB)")
    for theta_b, phi_b in beams_theta_phi_deg:
        tb, pb = np.radians(theta_b), np.radians(phi_b)
        ax.plot(
            np.sin(tb) * np.cos(pb),
            np.sin(tb) * np.sin(pb),
            marker="+",
            color="white",
            markersize=12,
            markeredgewidth=2,
        )
    ax.set_xlabel("u")
    ax.set_ylabel("v")
    ax.set_title(title)
    ax.set_aspect("equal")
    return _save_png(fig, path)


def pattern_cut_png(theta_deg, phi_deg, power_db, cut_phi_deg, path, title=None):
    """One principal cut (nearest available phi row) through the pattern: power_db vs theta."""
    phi_row = phi_deg[0, :]
    j = int(np.argmin(np.abs(phi_row - (cut_phi_deg % 360.0))))
    theta_col = theta_deg[:, j]
    cut_db = power_db[:, j]
    fig, ax = _new_fig(figsize=(5.5, 3.5))
    ax.plot(theta_col, cut_db, color=_ACCENT, linewidth=1.6)
    ax.set_xlabel("θ (°)")
    ax.set_ylabel("relative power (dB)")
    ax.set_title(title or f"φ = {phi_row[j]:.0f}° cut")
    finite = cut_db[np.isfinite(cut_db)]
    if finite.size:
        peak = float(np.max(finite))
        ax.set_ylim(peak - 50.0, peak + 3.0)  # a pattern null (e.g. 0 at 90 deg) must not squash the plot
    return _save_png(fig, path)


def vortex_phase_ring_png(alpha_deg, phase_deg, winding_number, path, title=None):
    """PHASE15 Part B: phase of the beam's own co-pol component around a small ring encircling its
    direction (`pattern.vector_aperture_pattern.measure_vortex_winding`'s own ring) -- a straight
    line means a pencil beam (no singularity inside), a `winding_number` full cycles of sawtooth
    means a vortex of that OAM charge. `alpha_deg`/`phase_deg`: same length 1-D arrays, the ring
    parameter and the unwrapped phase at each sample."""
    fig, ax = _new_fig(figsize=(5.5, 3.5))
    ax.plot(alpha_deg, phase_deg, color=_ACCENT, linewidth=1.6, marker=".", markersize=3)
    ax.set_xlabel("angle around the beam direction (°)")
    ax.set_ylabel("phase (°)")
    ax.set_title(title or f"Phase around beam axis -- winding {winding_number:+.1f}")
    return _save_png(fig, path)


_COPPER = "#E8893A"
_COPPER_EDGE = "#C8651B"
_SUBSTRATE = "#D9DADC"
_SUBSTRATE_EDGE = "#B8B9BC"


def patch_layout_png(cells_mm, period_mm, nx, ny, path, title=None, feed_xy_mm=(0.0, 0.0)):
    """True-geometry board drawing (PHASE13D item 5): every patch at its real L x W and position,
    bright copper on a light-grey substrate, feed point marked, mm axes, size legend.
    `cells_mm`: iterable of (x_mm, y_mm, l_mm, w_mm) patch centres/sizes."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    cells = [c for c in cells_mm if c[2] is not None and c[3] is not None]
    xs = np.array([c[0] for c in cells]) if cells else np.zeros(1)
    ys = np.array([c[1] for c in cells]) if cells else np.zeros(1)
    half_x = nx * period_mm / 2.0
    half_y = ny * period_mm / 2.0
    cx = (xs.min() + xs.max()) / 2.0
    cy = (ys.min() + ys.max()) / 2.0
    fig, ax = plt.subplots(figsize=(8, 8), dpi=260)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.add_patch(Rectangle((cx - half_x, cy - half_y), 2 * half_x, 2 * half_y,
                           facecolor=_SUBSTRATE, edgecolor=_SUBSTRATE_EDGE, linewidth=0.8, zorder=0))
    rects = [Rectangle((x - l / 2.0, y - w / 2.0), l, w) for x, y, l, w in cells]
    ax.add_collection(PatchCollection(rects, facecolor=_COPPER, edgecolor=_COPPER_EDGE, linewidth=0.15, zorder=1))
    if feed_xy_mm is not None:
        ax.plot(feed_xy_mm[0], feed_xy_mm[1], "o", markersize=4, color=_ACCENT, zorder=3)
    ax.set_xlim(cx - half_x - 2, cx + half_x + 2)
    ax.set_ylim(cy - half_y - 2, cy + half_y + 2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    for spine in ax.spines.values():
        spine.set_visible(False)
    if cells:
        sizes = np.array([c[2] for c in cells])
        ax.text(0.99, -0.07, f"Patch size {sizes.min():.2f}–{sizes.max():.2f} mm  ·  {len(cells)} patches  ·  pitch {period_mm:g} mm",
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555")
    ax.set_title(title or f"Patch layout — {nx} × {ny} cells, {period_mm:g} mm pitch", fontsize=10)
    return _save_png(fig, path)
