"""Phase 8B Sec. 50: the 5 required plots, each strictly a picture of an already-implemented,
numerically-tested calculation (Sec. 50: "only create plots that correspond to implemented/
validated calculations") -- none of them plots a "validated" holographic result, since none
exists yet; they visualize the object-wave math and the (mesh-derived) distinct-phase-state count."""

import numpy as np


def _require_matplotlib():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "matplotlib is required for the Phase 8B plots (src/phase8b/plots.py); install with "
            "'pip install -r requirements-cad.txt'."
        ) from e
    return plt


def _grid(extent_m=0.05, n=81):
    axis = np.linspace(-extent_m, extent_m, n)
    return np.meshgrid(axis, axis)


def plot_object_wave_single_beam(beam, frequency_hz, path):
    from .object_wave import single_beam_object_field

    plt = _require_matplotlib()
    x, y = _grid()
    field = single_beam_object_field(x, y, beam, frequency_hz)
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    mesh = ax.pcolormesh(
        x * 1e3, y * 1e3, np.degrees(np.angle(field)), shading="nearest", cmap="twilight", vmin=-180, vmax=180
    )
    fig.colorbar(mesh, ax=ax, label="phase(E_obj) [deg]")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title(f"Single-beam object wave (theta={beam.theta_deg:g}, phi={beam.phi_deg:g} deg)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_object_wave_multi_beam(beams, frequency_hz, path):
    from .object_wave import combine_beams

    plt = _require_matplotlib()
    x, y = _grid()
    ex, ey = combine_beams(x, y, beams, frequency_hz)
    field = ex if beams[0].polarization == "x" else ey
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    mesh = ax.pcolormesh(
        x * 1e3, y * 1e3, np.degrees(np.angle(field)), shading="nearest", cmap="twilight", vmin=-180, vmax=180
    )
    fig.colorbar(mesh, ax=ax, label="phase(E_obj) [deg]")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title(f"Multi-beam object wave ({len(beams)} beams, field-level superposition)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_phase_gradient_test(frequency_hz, path):
    """Sec. 37/39: broadside (theta=0) gives zero transverse phase gradient; phi=0/90/180/270 at
    a fixed nonzero theta give gradients along +x/+y/-x/-y respectively -- a direct visual of the
    test `tests/test_phase8b_object_wave.py` checks numerically."""
    from .object_wave import Beam, single_beam_object_field

    plt = _require_matplotlib()
    x, y = _grid(extent_m=0.03, n=61)
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    cases = [("broadside", 0.0, 0.0)] + [(f"phi={p}", 30.0, p) for p in (0, 90, 180, 270)]
    for ax, (label, theta, phi) in zip(axes, cases):
        beam = Beam(theta_deg=theta, phi_deg=phi)
        field = single_beam_object_field(x, y, beam, frequency_hz)
        ax.pcolormesh(
            x * 1e3, y * 1e3, np.degrees(np.angle(field)), shading="nearest", cmap="twilight", vmin=-180, vmax=180
        )
        ax.set_title(label, fontsize=9)
        ax.set_aspect("equal")
    fig.suptitle("Object-wave phase-gradient test (broadside + 4 azimuths, theta=30 deg)")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_phase_control_quantity(evaluations, path):
    plt = _require_matplotlib()
    from . import control_variable as cv

    rows = [f() for f in cv.ALL_CANDIDATES]
    colors = {"VALIDATED": "tab:green", "INVALID": "tab:red", "INSUFFICIENT_EVIDENCE": "tab:orange"}
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    names = [r.name for r in rows]
    ax.barh(names, [1] * len(rows), color=[colors[r.validated] for r in rows])
    for i, r in enumerate(rows):
        ax.text(0.02, i, r.validated, va="center", fontsize=9, color="white")
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("Candidate physical control variables")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_distinct_phase_states(analysis, path):
    plt = _require_matplotlib()
    group_sizes = {}
    for row in analysis.rows:
        group_sizes[row.equivalence_group] = group_sizes.get(row.equivalence_group, 0) + 1
    sizes = sorted(group_sizes.values(), reverse=True)
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.bar(range(1, len(sizes) + 1), sizes)
    ax.set_xlabel("equivalence group (sorted by size)")
    ax.set_ylabel("geometries sharing this phase state")
    ax.set_title(f"{analysis.n_geometries} geometries -> {analysis.n_distinct_phase_states} distinct phase states")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
