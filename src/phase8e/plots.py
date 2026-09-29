"""Sec. 23: the required plots, one function per file, all with physical units/axes labelled."""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 8E plots; install with "
            "'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def plot_reference_field_amplitude_phase(field, path):
    plt = _require_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    im0 = axes[0].pcolormesh(field.x_m * 1e3, field.y_m * 1e3, np.abs(field.Ex), shading="auto")
    axes[0].set_title("|E_ref| (Ex)"); axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    fig.colorbar(im0, ax=axes[0], label="V/m")
    im1 = axes[1].pcolormesh(field.x_m * 1e3, field.y_m * 1e3, np.degrees(np.angle(field.Ex)), shading="auto", cmap="twilight")
    axes[1].set_title("phase(E_ref) [deg]"); axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
    fig.colorbar(im1, ax=axes[1], label="deg")
    fig.suptitle(f"Reference field ({field.source})")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_reference_field_components(field, path):
    plt = _require_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, comp, name in ((axes[0], field.Ex, "Ex"), (axes[1], field.Ey, "Ey")):
        im = ax.pcolormesh(field.x_m * 1e3, field.y_m * 1e3, np.abs(comp), shading="auto")
        ax.set_title(f"|{name}|"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        fig.colorbar(im, ax=ax, label="V/m")
    fig.suptitle("Reference field components (Ez not modelled: scalar-feed synthetic field)")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_reference_field_spectrum(field, spectral_result, path):
    plt = _require_matplotlib()
    spectrum = np.fft.fftshift(np.fft.fft2(field.Ex))
    ny, nx = field.Ex.shape
    kx = np.fft.fftshift(np.fft.fftfreq(nx, d=field.dx_m)) * 2 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(ny, d=field.dy_m)) * 2 * np.pi
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.pcolormesh(kx, ky, 10 * np.log10(np.abs(spectrum) ** 2 + 1e-30), shading="auto")
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(spectral_result.k0_per_m * np.cos(theta), spectral_result.k0_per_m * np.sin(theta), "r--",
            label="k0 (propagating boundary)")
    ax.set_xlabel("kx [rad/m]"); ax.set_ylabel("ky [rad/m]")
    ax.set_title(f"Reference-field 2-D spectrum (propagating fraction = {spectral_result.propagating_fraction:.4g})")
    ax.legend(); fig.colorbar(im, ax=ax, label="dB")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_object_wave_amplitude_phase(x_m, y_m, e_obj, path):
    plt = _require_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    im0 = axes[0].pcolormesh(x_m * 1e3, y_m * 1e3, np.abs(e_obj), shading="auto")
    axes[0].set_title("|E_obj|"); axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].pcolormesh(x_m * 1e3, y_m * 1e3, np.degrees(np.angle(e_obj)), shading="auto", cmap="twilight")
    axes[1].set_title("phase(E_obj) [deg]"); axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
    fig.colorbar(im1, ax=axes[1])
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_object_wave_spectrum(object_spectrum, k0_per_m, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 6))
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.plot(k0_per_m * np.cos(theta), k0_per_m * np.sin(theta), "k--", label="k0 (visible-region boundary)")
    for p in object_spectrum.beams:
        ax.plot(p.kx_per_m, p.ky_per_m, "o", label=f"theta={p.theta_deg:g}, phi={p.phi_deg:g}")
    ax.set_xlabel("kx [rad/m]"); ax.set_ylabel("ky [rad/m]"); ax.set_aspect("equal")
    ax.set_title("Object-wave beam transverse wavevectors")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_control_response_amplitude_phase(x_m, y_m, complex_values, title, path):
    plt = _require_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    im0 = axes[0].pcolormesh(x_m * 1e3, y_m * 1e3, np.abs(complex_values), shading="auto")
    axes[0].set_title(f"|{title}|"); axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].pcolormesh(x_m * 1e3, y_m * 1e3, np.degrees(np.angle(complex_values)), shading="auto", cmap="twilight")
    axes[1].set_title(f"phase({title}) [deg]"); axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
    fig.colorbar(im1, ax=axes[1])
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_control_response_complex_plane(complex_values, title, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    values = np.asarray(complex_values).ravel()
    ax.scatter(values.real, values.imag, s=8)
    ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_title(f"C_req complex-plane distribution ({title})")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_library_complex_state_scatter(rows, component_name, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    reals = [r.real for r in rows]; imags = [r.imag for r in rows]
    ax.scatter(reals, imags, s=10)
    ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_title(f"Library complex-state scatter ({component_name})")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_library_magnitude_vs_geometry(rows, component_name, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter([r.l_mm for r in rows], [r.w_mm for r in rows], c=[r.mag for r in rows], cmap="viridis")
    ax.set_xlabel("L [mm]"); ax.set_ylabel("W [mm]")
    ax.set_title(f"|{component_name}| vs geometry")
    fig.colorbar(sc, ax=ax, label="magnitude")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_library_phase_vs_geometry(rows, component_name, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter([r.l_mm for r in rows], [r.w_mm for r in rows], c=[r.phase_deg for r in rows], cmap="twilight")
    ax.set_xlabel("L [mm]"); ax.set_ylabel("W [mm]")
    ax.set_title(f"phase({component_name}) vs geometry")
    fig.colorbar(sc, ax=ax, label="deg")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_amplitude_phase_relationship(rows, component_name, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter([r.phase_deg for r in rows], [r.mag for r in rows], s=10)
    ax.set_xlabel("phase [deg]"); ax.set_ylabel("magnitude")
    ax.set_title(f"Amplitude-phase coupling ({component_name})")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_aperture_field(x_m, y_m, aperture_field, title, path):
    plt = _require_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    im0 = axes[0].scatter(x_m * 1e3, y_m * 1e3, c=np.abs(aperture_field), s=40)
    axes[0].set_title(f"|{title}|"); axes[0].set_xlabel("x [mm]"); axes[0].set_ylabel("y [mm]")
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].scatter(x_m * 1e3, y_m * 1e3, c=np.degrees(np.angle(aperture_field)), s=40, cmap="twilight")
    axes[1].set_title(f"phase({title}) [deg]"); axes[1].set_xlabel("x [mm]"); axes[1].set_ylabel("y [mm]")
    fig.colorbar(im1, ax=axes[1])
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_far_field(theta_grid_deg, phi_grid_deg, far_field, title, path):
    plt = _require_matplotlib()
    power_db = 10 * np.log10(np.abs(far_field) ** 2 / np.max(np.abs(far_field) ** 2) + 1e-12)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.pcolormesh(phi_grid_deg, theta_grid_deg, power_db, shading="auto", vmin=-40, vmax=0)
    ax.set_xlabel("phi [deg]"); ax.set_ylabel("theta [deg]")
    ax.set_title(f"Far-field pattern, relative power [dB] ({title})")
    fig.colorbar(im, ax=ax, label="dB")
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_beam_cuts(theta_axis_deg, continuous_cut_db, discrete_cut_db, phi_fixed_deg, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(theta_axis_deg, continuous_cut_db, label="continuous")
    ax.plot(theta_axis_deg, discrete_cut_db, label="discrete", linestyle="--")
    ax.set_xlabel("theta [deg]"); ax.set_ylabel("relative power [dB]")
    ax.set_title(f"Beam cut at phi={phi_fixed_deg:g} deg")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def plot_beam_pointing_comparison(intended_beams, continuous_metrics, discrete_metrics, path):
    plt = _require_matplotlib()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    idx = np.arange(len(intended_beams))
    width = 0.35
    ax.bar(idx - width / 2, continuous_metrics.pointing_error_deg, width, label="continuous")
    ax.bar(idx + width / 2, discrete_metrics.pointing_error_deg, width, label="discrete")
    ax.set_xticks(idx); ax.set_xticklabels([f"beam {i+1}\n({b.theta_deg:g},{b.phi_deg:g})" for i, b in enumerate(intended_beams)])
    ax.set_ylabel("pointing error [deg]")
    ax.set_title("Beam-pointing error: continuous vs. discrete reconstruction")
    ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
