"""Phase 26 A: a leakage-aware forward model for ONE monopole-fed hologram (additive; the Phase 21-25
cell-wise pattern path is untouched and every existing result is bit-identical).

The existing forward model (`pattern.cellwise_pattern`, `hologram.butler_forward`) lets every cell radiate
`E_inc x Gamma_state` with a phase-only (unit-modulus) incident field and puts ALL of the monopole's power
through the aperture: its directivity is a pure aperture-illumination number, silent about how much of the
surface wave is actually converted. This module adds the missing accounting, from `dispersion.leakage_model`:

* a local leakage constant `alpha(psi)` along every ray from the feed (the pattern is a two-state hologram,
  so the modulation depth is the same everywhere; `alpha` depends on the ray's azimuth `psi` only through the
  TM/TE split of the radiated harmonic, `zeta = phi_beam - psi`, and through the local modulation period);
* the surface wave decays as `exp(-alpha rho)` in amplitude and spreads as `rho^-1/2`: every cell's incident
  field is multiplied by `rho^-1/2 exp(-alpha rho)`, the aperture illumination that results;
* the fraction of the launched surface-wave power the aperture converts,
  `eta_conv(psi) = 1 - exp(-2 alpha(psi) rho_edge(psi))`, averaged over `psi` (the monopole launches a
  circularly symmetric wave), and the unconverted remainder reaching the aperture edge, reported, not
  hidden. Nothing about what the remainder does at the edge is modelled (no edge diffraction): it is
  assumed to radiate OUTSIDE the main beam, which is an assumption of this model, not a result.

`launch_fraction` is the share of the accepted power that enters the surface wave at all (the rest is the
monopole's direct space-wave radiation). It is NOT computed here -- there is no vertical-pin launching model
in this repo -- and defaults to 1.0, i.e. every prediction is an UPPER BOUND in that respect.

Gain bookkeeping: `directivity_ap` is the directivity of the aperture pattern alone (normalised by its own
radiated power, as every earlier number), and the predicted gain referred to the accepted power is
`directivity_ap * launch_fraction * eta_conv`.
"""

from dataclasses import dataclass

import numpy as np


def polygon_ray_distance(vertices_m, centre_xy_m, psi_rad):
    """Distance from `centre_xy_m` (inside the convex polygon `vertices_m`) to the polygon boundary along
    the ray at azimuth `psi_rad` (array)."""
    v = np.asarray(vertices_m, dtype=float)
    c = np.asarray(centre_xy_m, dtype=float)
    psi = np.atleast_1d(np.asarray(psi_rad, dtype=float))
    d = np.stack([np.cos(psi), np.sin(psi)], axis=-1)  # (P, 2)
    best = np.full(psi.shape, np.inf)
    for a, b in zip(v, np.roll(v, -1, axis=0), strict=True):
        e = b - a
        # solve c + t d = a + s e
        den = d[:, 0] * (-e[1]) - d[:, 1] * (-e[0])
        with np.errstate(divide="ignore", invalid="ignore"):
            t = ((a[0] - c[0]) * (-e[1]) - (a[1] - c[1]) * (-e[0])) / den
            s = (d[:, 0] * (a[1] - c[1]) - d[:, 1] * (a[0] - c[0])) / den
        ok = (t > 0.0) & (s >= -1e-12) & (s <= 1.0 + 1e-12) & np.isfinite(t)
        best = np.where(ok & (t < best), t, best)
    return best


@dataclass(frozen=True)
class LeakageProfile:
    """`alpha_of_psi`: callable azimuth [rad] -> alpha [Np/m]; the per-cell amplitude taper and the converted
    / unconverted power accounting for one feed and aperture."""

    amplitude: np.ndarray  # (N,) rho^-1/2 exp(-alpha rho) per cell, unnormalised
    alpha_cell: np.ndarray  # (N,)
    psi_grid_rad: np.ndarray
    alpha_grid: np.ndarray
    rho_edge_grid_m: np.ndarray
    eta_conv_grid: np.ndarray
    eta_conv: float  # converted fraction of the launched surface-wave power, averaged over psi
    edge_fraction: float  # 1 - eta_conv: reaches the aperture edge unconverted


def leakage_profile(x_m, y_m, feed_xy_m, vertices_m, alpha_of_psi, period_m, n_psi=720):
    """The `LeakageProfile` of the cells `(x_m, y_m)` of one aperture (polygon outline `vertices_m`) fed
    at `feed_xy_m`."""
    x = np.asarray(x_m, dtype=float).ravel()
    y = np.asarray(y_m, dtype=float).ravel()
    dx, dy = x - feed_xy_m[0], y - feed_xy_m[1]
    rho = np.hypot(dx, dy)
    psi = np.arctan2(dy, dx)
    alpha_cell = np.asarray(alpha_of_psi(psi), dtype=float)
    amp = np.maximum(rho, 0.5 * period_m) ** -0.5 * np.exp(-alpha_cell * rho)
    psi_g = np.linspace(-np.pi, np.pi, n_psi, endpoint=False)
    alpha_g = np.asarray(alpha_of_psi(psi_g), dtype=float)
    redge = polygon_ray_distance(vertices_m, feed_xy_m, psi_g)
    eta_g = 1.0 - np.exp(-2.0 * alpha_g * redge)
    eta = float(np.mean(eta_g))
    return LeakageProfile(amp, alpha_cell, psi_g, alpha_g, redge, eta_g, eta, 1.0 - eta)


def scaled_board(board, amplitude):
    """A `butler_forward.ForwardBoard` whose incident field of every feed is multiplied by the per-cell
    `amplitude` (the original untouched)."""
    from .butler_forward import ForwardBoard

    return ForwardBoard(
        board.x_m,
        board.y_m,
        board.gamma,
        board.e_inc * np.asarray(amplitude)[None, :],
        board.feeds_xy_m,
        board.period_m,
        board.frequency_hz,
        board.substrate,
        board.polarization,
        board.label + "+leakage",
        board.amplitude_model,
        board.k_swg_rad_per_m,
        dict(board.meta),
    )


def pattern_budget(theta_deg, phi_deg, e_theta, e_phi, beam_theta_deg, beam_phi_deg, cones_deg=(8.0, 20.0)):
    """Aperture-only pattern metrics on a regular `(theta, phi)` hemisphere grid (`theta` and `phi` 1-D,
    fields `(n_theta, n_phi)`): directivity of the aperture pattern alone, beam peak, RHCP/LHCP split, and the
    fractions of the aperture's radiated power inside cones about the beam / above theta = 70 deg."""
    from ..pattern.vector_aperture_pattern import copol_component

    tg, pg = np.meshgrid(np.asarray(theta_deg), np.asarray(phi_deg), indexing="ij")
    rhcp = np.abs(copol_component(e_theta, e_phi, "rhcp")) ** 2
    lhcp = np.abs(copol_component(e_theta, e_phi, "lhcp")) ** 2
    total = np.abs(e_theta) ** 2 + np.abs(e_phi) ** 2
    dth = np.radians(theta_deg[1] - theta_deg[0])
    dph = np.radians(phi_deg[1] - phi_deg[0])
    w = np.sin(np.radians(tg)) * dth * dph
    p_rad = float(np.sum(total * w))
    t0, p0 = np.radians(beam_theta_deg), np.radians(beam_phi_deg)
    cosd = np.cos(np.radians(tg)) * np.cos(t0) + np.sin(np.radians(tg)) * np.sin(t0) * np.cos(np.radians(pg) - p0)
    dist = np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0)))
    near = dist <= 5.0
    i = np.unravel_index(np.argmax(np.where(near, rhcp, -1.0)), rhcp.shape)
    out = {
        "directivity_ap_rhcp_dbi": float(10 * np.log10(4 * np.pi * rhcp[i] / p_rad)),
        "directivity_ap_total_dbi": float(10 * np.log10(4 * np.pi * total.max() / p_rad)),
        "peak_theta_deg": float(tg[i]),
        "peak_phi_deg": float(pg[i]),
        "xpd_at_peak_db": float(10 * np.log10(rhcp[i] / max(lhcp[i], 1e-300))),
        "rhcp_fraction": float(np.sum(rhcp * w) / p_rad),
        "lhcp_fraction": float(np.sum(lhcp * w) / p_rad),
        "above_70deg_fraction": float(np.sum(total * w * (tg > 70.0)) / p_rad),
    }
    for c in cones_deg:
        out[f"within_{c:g}deg_fraction"] = float(np.sum(total * w * (dist <= c)) / p_rad)
        out[f"within_{c:g}deg_rhcp_fraction"] = float(np.sum(rhcp * w * (dist <= c)) / p_rad)
    return out
