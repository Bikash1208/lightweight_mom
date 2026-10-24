"""Phase 26 C: pure-RHCP rotating-slot tensor hologram on one monopole-fed aperture, designed from the
SOLVER_SUPERCELL table (`scripts/research/phase26_supercell_table.py`). Additive; nothing else is touched.

Model (labels: SOLVER_SUPERCELL, MODEL_ONLY, NOT_HFSS_VALIDATED; the table is collinear, gamma = 0 only):
* The launched wave is cylindrical, the same power per radian `W0` on every ray; along a ray its power per
  radian is `W(rho)`, `dW/drho = -2 alpha W`, with `alpha` the table's `alpha_leak` of the cell chosen at
  that point and the radiated power per unit area `q = 2 alpha W / rho` (energy bookkeeping: radiated +
  remaining = launched, tested).
* Reference phase `phi_ref = int beta drho` along the ray with the table's `beta` (modulated row, complex root
  real part) of the chosen cell. Modulation wavenumber at a point `K = |beta rho_hat - k_t|`, `k_t =
  k0 sin(theta0) (cos phi0, sin phi0)`; the table column is the NEAREST K of the grid (no interpolation
  across periods, as decided), the slot is the table slot whose `alpha_leak` is nearest the required one.
* Rotating-slot cell = rotating-axis tensor sheet; the surface-wave E (along rho_hat) is reflected about the
  slot axis `a`: the induced current direction is `chi = 2 a - phi_ray`, `J = A exp(-j phi_ref) (cos chi, sin chi)`
  with `A = sqrt(2 alpha W / rho)`. `(cos, sin) = exp(j chi)(1,-j)/2 + exp(-j chi)(1,+j)/2`: choosing
  `chi = hand (phi_ref - k_t . r)` puts one circular component in phase with the beam and the other at
  wavevector `2 beta rho_hat - k_t`, |.| > 1.5 k0 -- evanescent, so it does not radiate. Purity is
  COMPUTED from the far field of the resulting currents, not assumed.
* The aperture pattern is the (E_theta, E_phi) of `J` on the grounded slab (`far_field` machinery of the
  Phase 24-26 forward model). Gain referred to the accepted power = `D_ap * launch_fraction * eta_conv`;
  the unconverted remainder and the pin's direct space-wave are NOT given a pattern (assumed outside the
  beam), reported as budget lines.
"""

import glob
import json
from dataclasses import dataclass

import numpy as np

from .leaky_forward import polygon_ray_distance


@dataclass(frozen=True)
class SupercellTable:
    slots: tuple  # ((L_mm, W_mm), ...) ordered by slot area, weakest first
    k_grid: np.ndarray  # K / k0, ascending
    alpha: np.ndarray  # (n_slots, n_k) alpha_leak [Np/m]; nan where missing
    beta: np.ndarray  # (n_slots, n_k) beta / k0
    diameter_mm: float
    pitch_mm: float

    @classmethod
    def from_files(cls, pattern):
        pts = []
        meta = None
        for f in sorted(glob.glob(pattern)):
            d = json.load(open(f, encoding="utf-8"))
            meta = meta or d
            for p in d["points"]:
                if "alpha_leak_np_per_m" in p:
                    pts.append(((float(p["slot_L_W_mm"][0]), float(p["slot_L_W_mm"][1])), p))
        slots = sorted({s for s, _ in pts}, key=lambda s: s[0] * s[1])
        ks = sorted({round(p["K_over_k0"], 4) for _, p in pts})
        a = np.full((len(slots), len(ks)), np.nan)
        b = np.full_like(a, np.nan)
        for s, p in pts:
            i, j = slots.index(s), ks.index(round(p["K_over_k0"], 4))
            a[i, j], b[i, j] = p["alpha_leak_np_per_m"], p["beta_over_k0"]
        return cls(tuple(slots), np.array(ks), a, b, float(meta["diameter_mm"]) if "diameter_mm" in meta else 6.0,
                   float(meta.get("pitch_mm", 6.4)))

    def column(self, k_over_k0):
        """Nearest-grid column index for `K / k0` (array)."""
        k = np.asarray(k_over_k0, dtype=float)
        return np.argmin(np.abs(k[..., None] - self.k_grid), axis=-1)


def modulation_k(beta_over_k0, phi_ray, theta0_rad, phi0_rad):
    """`K / k0` of the radial wave `beta rho_hat` against the beam's transverse wavevector."""
    kt = np.sin(theta0_rad)
    return np.sqrt(
        beta_over_k0**2 + kt**2 - 2.0 * beta_over_k0 * kt * np.cos(phi_ray - phi0_rad)
    )


def march_rays(table, vertices_m, feed_xy_m, k0, theta0_rad, phi0_rad, eta_design, taper_c=0.0,
               rho_min_m=0.007, n_psi=720, drho_m=0.25e-3, beta0_over_k0=1.065, density_fn=None):
    """Ray-by-ray design march (module docstring). Returns a dict of `(n_psi, n_rho)` arrays plus the ray
    axis. `eta_design`: the conversion the target profile asks for (the realised one is lower where the
    table saturates). `taper_c`: amplitude taper `g = 1 - c (rho/rho_edge)^2` of the target density.
    `density_fn(psi, rho, r_edge) -> (n_psi, n_rho)` (optional, additive): an extra amplitude factor multiplying `g`
    (default None: unchanged, bit-identical))."""
    psi = np.linspace(-np.pi, np.pi, n_psi, endpoint=False)
    r_edge = polygon_ray_distance(vertices_m, feed_xy_m, psi)
    n_rho = int(np.ceil(r_edge.max() / drho_m)) + 1
    rho = np.arange(n_rho) * drho_m
    x = rho[None, :] / r_edge[:, None]
    g = np.clip(1.0 - taper_c * x**2, 0.0, None)
    if density_fn is not None:
        g = g * np.asarray(density_fn(psi, rho, r_edge), dtype=float)
    g2 = g**2
    inside = rho[None, :] <= r_edge[:, None]
    # q0 per ray so that the UNSATURATED target converts eta_design: int q rho drho = eta W0
    w0 = 1.0
    i_g = np.sum(g2 * inside * rho[None, :], axis=1) * drho_m
    q0 = eta_design * w0 / np.maximum(i_g, 1e-30)
    alpha = np.zeros((n_psi, n_rho))
    beta = np.zeros_like(alpha)
    w = np.ones_like(alpha)
    phi_ref = np.zeros_like(alpha)
    slot = np.full((n_psi, n_rho), -1, dtype=int)
    kk = np.zeros_like(alpha)
    amax_col = np.nanmax(table.alpha, axis=0)
    cur_w = np.ones(n_psi)
    cur_phi = np.zeros(n_psi)
    cur_beta = np.full(n_psi, beta0_over_k0)
    for k in range(n_rho):
        r = rho[k]
        kvec = modulation_k(cur_beta, psi, theta0_rad, phi0_rad)
        col = table.column(kvec)
        q = q0 * g2[:, k]
        a_t = np.where(cur_w > 1e-12, q * max(r, 1e-9) / (2.0 * np.maximum(cur_w, 1e-12)), amax_col[col])
        if r < rho_min_m:
            a_t = np.zeros(n_psi)
        a_col = table.alpha[:, col]  # (n_slots, n_psi)
        pick = np.nanargmin(np.abs(np.nan_to_num(a_col, nan=1e9) - a_t[None, :]), axis=0)
        a_used = a_col[pick, np.arange(n_psi)]
        b_used = table.beta[pick, col]
        if r < rho_min_m:
            a_used = np.zeros(n_psi)
        a_used = np.where(inside[:, k], a_used, 0.0)
        alpha[:, k], beta[:, k], w[:, k], phi_ref[:, k], slot[:, k], kk[:, k] = a_used, b_used, cur_w, cur_phi, pick, kvec
        cur_w = cur_w * np.exp(-2.0 * a_used * drho_m)
        cur_phi = cur_phi + b_used * k0 * drho_m
        cur_beta = b_used
    w_edge = np.array([w[i, min(int(np.floor(r_edge[i] / drho_m)), n_rho - 1)] for i in range(n_psi)])
    return {"psi": psi, "rho": rho, "r_edge": r_edge, "alpha": alpha, "beta": beta, "w": w, "phi_ref": phi_ref,
            "slot": slot, "K": kk, "w_edge": w_edge, "eta_conv": float(np.mean(1.0 - w_edge)),
            "eta_conv_by_ray": 1.0 - w_edge, "drho": drho_m}


def cell_design(march, table, x_m, y_m, feed_xy_m, k0, theta0_rad, phi0_rad, hand=+1):
    """Per-cell slot (L, W, angle) and aperture current `(Jx, Jy)` from `march`."""
    from scipy.interpolate import RegularGridInterpolator

    dx, dy = x_m - feed_xy_m[0], y_m - feed_xy_m[1]
    rho = np.hypot(dx, dy)
    phi = np.arctan2(dy, dx)
    psi, rg = march["psi"], march["rho"]
    # periodic in psi: pad one ray each side
    pp = np.concatenate([[psi[-1] - 2 * np.pi], psi, [psi[0] + 2 * np.pi]])

    def interp(field):
        f = np.vstack([field[-1:], field, field[:1]])
        return RegularGridInterpolator((pp, rg), f, bounds_error=False, fill_value=None)((phi, np.minimum(rho, rg[-1])))

    w = interp(march["w"])
    phi_ref = interp(march["phi_ref"])
    # nearest ray / step for the discrete choices
    ip = np.argmin(np.abs(((phi[:, None] - psi[None, :] + np.pi) % (2 * np.pi)) - np.pi), axis=1)
    ir = np.minimum(np.round(rho / march["drho"]).astype(int), len(rg) - 1)
    slot = march["slot"][ip, ir]
    alpha = march["alpha"][ip, ir]
    beta = march["beta"][ip, ir]
    amp = np.sqrt(np.maximum(2.0 * alpha * np.maximum(w, 0.0) / np.maximum(rho, 1e-9), 0.0))
    kt = k0 * np.sin(theta0_rad)
    kdotr = kt * (dx * np.cos(phi0_rad) + dy * np.sin(phi0_rad))
    chi = hand * (phi_ref - kdotr)
    axis = 0.5 * (chi + phi)  # slot axis: chi = 2 a - phi_ray
    j = amp * np.exp(-1j * phi_ref)
    jx, jy = j * np.cos(chi), j * np.sin(chi)
    L = np.array([table.slots[s][0] for s in slot])
    W = np.array([table.slots[s][1] for s in slot])
    return {"slot_L_mm": L, "slot_W_mm": W, "slot_angle_rad": np.mod(axis, np.pi), "alpha": alpha, "beta_over_k0": beta,
            "amp": amp, "jx": jx, "jy": jy, "rho": rho, "phi": phi, "w": w, "phi_ref": phi_ref, "slot_index": slot}
