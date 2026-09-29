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
    ax.set_xlabel("theta (deg)")
    ax.set_ylabel("relative power (dB)")
    ax.set_title(title or f"Cut at phi={phi_row[j]:.0f}deg")
    return _save_png(fig, path)
