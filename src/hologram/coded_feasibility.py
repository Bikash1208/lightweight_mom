"""Phase 28 item 2: feasibility of a CODED (quantised-state) monopole-fed hologram on the existing Phase 26 forward model.

Nothing here is a new physics model. The cell design (`tensor_leaky_cp.march_rays` / `cell_design`), the aperture far field
(`pattern.vector_aperture_pattern._vector_far_field` on the grounded RT/Duroid slab) and the supercell table are the Phase 26 ones, imported
unchanged. This module only (a) quantises the per-cell state the continuous design asks for, (b) re-evaluates the pattern, (c) derives the
switching contrast a coded cell needs from the same closed-form leakage model the Phase 26 required map used, and (d) counts diodes, bias lines
and driver current. Every result is MODEL_ONLY: the cells are IDEALISED (a state is a tensor / scalar value, no geometry), the supercell table
is the 6.4 mm, 10.54 GHz one and is assumed to hold at every pitch (only the SAMPLING of the hologram changes with pitch).

Two coded families, both read from the continuous design's per-cell quantities (`amp`, `phi_ref`, the cell azimuth `phi`):

* CP, rotating-slot tensor (the Phase 26 family). The state is the slot AXIS angle `a` (period pi: a slot rotated by pi is the same slot). The
  induced current direction is `chi = 2a - phi` (period 2 pi), so N axis states are N equally spaced axes; `quantise_circular` wraps on the circle.
  1-bit = 2 axes (0, 90 deg), 2-bit = 4 axes (0, 45, 90, 135 deg). The continuous design is recovered exactly with no quantisation.
* LP, scalar. Scalar surface reactance m(r) = cos(phi_ref - k_t . r) in [-1, 1] times the radial current direction (the surface wave's E is along rho):
  J = amp exp(-j phi_ref) m (cos phi, sin phi). m is quantised to 2**bits equally spaced levels INCLUDING both end points (1-bit: -1, +1;
  2-bit: -1, -1/3, +1/3, +1), nearest level. The twin term (wavevector 2 beta rho_hat - k_t) is evanescent, as in the Phase 26 CP design.
  A scalar surface does not give a clean polarisation: only total-power metrics are reported for LP, plus the E_theta/E_phi ratio at the peak.

Held fixed from the continuous design (stated assumption, not a result): the surface-wave conversion `eta_conv` and the launch fraction; the gain is
`D_ap + 10 log10(LAUNCH_BARE eta_conv)`, so a quantisation changes D_ap only.
"""

import dataclasses
import math

import numpy as np
from scipy.optimize import brentq

from ..dispersion.leakage_model import (
    leakage_alpha_first_order,
    leakage_alpha_first_order_tensor,
    state_surface_reactances,
)
from ..dispersion.sheet_model import find_beta_from_sheet_admittance
from ..frequency import free_space_wavenumber
from ..green_function import LayeredGreenFunction
from ..pattern.vector_aperture_pattern import _vector_far_field, copol_component
from . import tensor_cp_design as t26
from .tensor_leaky_cp import cell_design, march_rays

F = t26.F
K0 = t26.K0
SUB = t26.SUB
LAUNCH_BARE = t26.LAUNCH_BARE
ETA_DESIGN, TAPER_C = 0.9, 0.5  # the selected Phase 26 single-hexagon design
BOUND_DB = 40.0  # AR / XPD values beyond this are reported as "> 40 dB", never as numbers


# ------------------------------------------------------------------------------------------------ quantisers
def uniform_levels(bits):
    """`2**bits` equally spaced levels on [-1, 1] including both end points (1-bit: -1, 1; 2-bit: -1, -1/3, 1/3, 1)."""
    return np.linspace(-1.0, 1.0, 2**bits)


def quantise_nearest(values, levels):
    """Nearest level; an exact tie goes to the level that comes first in `levels` (floating point: a value that is only nearly half-way is not a tie,
    e.g. 0.0 in the 2-bit set is 1 ulp nearer +1/3). Returns `(index, quantised)`. Non-finite input is an error."""
    v = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(v)):
        raise ValueError("quantise_nearest: non-finite input")
    lv = np.asarray(levels, dtype=float)
    idx = np.argmin(np.abs(v[..., None] - lv), axis=-1)
    return idx, lv[idx]


def quantise_circular(values, n_levels, period):
    """Nearest of `n_levels` equally spaced levels `k period / n_levels` on a circle of the given period (angles wrap: a value just below
    `period` goes to level 0). Exact half-way values go to the higher level. Returns `(index, quantised)` with `quantised` in [0, period)."""
    step = period / n_levels
    v = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(v)):
        raise ValueError("quantise_circular: non-finite input")
    idx = np.floor(v / step + 0.5).astype(int) % n_levels
    return idx, idx * step


def bounded_db(value_db, bound=BOUND_DB):
    """`{"value": x, "text": "x"}` or, beyond `bound`, `{"value": None, "bound": bound, "text": "> bound dB"}` (a clamped number is not a value)."""
    if value_db is None or not np.isfinite(value_db) or value_db > bound:
        return {"value": None, "bound_db": bound, "text": f"> {bound:g} dB"}
    return {"value": float(value_db), "bound_db": None, "text": f"{value_db:.2f}"}


# ------------------------------------------------------------------------------------------------ geometry and cells
def hexagon_m(r_m):
    """Hexagon of circumradius `r_m` with a vertex at 90 deg (the Phase 26 orientation)."""
    return np.array([[r_m * np.cos(np.radians(90 + 60 * i)), r_m * np.sin(np.radians(90 + 60 * i))] for i in range(6)])


def lattice(poly, pitch_m, keepout_m=t26.RHO_CLEAR):
    """Square lattice of pitch `pitch_m` with a node at the feed (origin), inside the convex polygon, outside the feed keep-out."""
    n = int(np.ceil(np.abs(poly).max() / pitch_m)) + 3
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(pitch_m * k, pitch_m * k)
    x, y = xg.ravel(), yg.ravel()
    ok = np.ones(x.shape, dtype=bool)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        ok &= (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) >= -1e-9
    ok &= np.hypot(x, y) >= keepout_m
    return x[ok], y[ok]


def design_cells(table, poly, pitch_m, theta0_deg, phi0_deg, hand=+1, n_psi=720):
    """The continuous Phase 26 design on a lattice of pitch `pitch_m`: `(x, y, d, march)`."""
    x, y = lattice(poly, pitch_m)
    th0, ph0 = np.radians(theta0_deg), np.radians(phi0_deg)
    m = march_rays(table, poly, (0.0, 0.0), K0, th0, ph0, ETA_DESIGN, taper_c=TAPER_C, n_psi=n_psi)
    d = cell_design(m, table, x, y, (0.0, 0.0), K0, th0, ph0, hand=hand)
    return x, y, d, m


def cp_currents(d, axis_rad):
    """CP rotating-slot currents for slot axis angles `axis_rad`: `chi = 2a - phi` (the Phase 26 relation), `J = amp exp(-j phi_ref) (cos chi, sin chi)`."""
    chi = 2.0 * np.asarray(axis_rad) - d["phi"]
    j = d["amp"] * np.exp(-1j * d["phi_ref"])
    return j * np.cos(chi), j * np.sin(chi)


def cp_states(d, bits):
    """Quantised slot axes: `(state_index, axis_rad)`. `bits=None`: the continuous design's own axes."""
    if bits is None:
        return None, d["slot_angle_rad"]
    return quantise_circular(d["slot_angle_rad"], 2**bits, np.pi)


def lp_modulation(d, x, y, theta0_deg, phi0_deg):
    """The continuous scalar modulation `m = cos(phi_ref - k_t . r)` in [-1, 1]."""
    th0, ph0 = np.radians(theta0_deg), np.radians(phi0_deg)
    kdotr = K0 * np.sin(th0) * (x * np.cos(ph0) + y * np.sin(ph0))
    return np.cos(d["phi_ref"] - kdotr)


def lp_currents(d, m):
    """Scalar surface: the induced current follows the surface wave's E, i.e. the radial direction (cos phi, sin phi)."""
    j = d["amp"] * np.exp(-1j * d["phi_ref"]) * m
    return j * np.cos(d["phi"]), j * np.sin(d["phi"])


def lpx_currents(d, m):
    """LP with IDEAL polarisation control (current along x for every cell, whatever the local E direction): the bound a polarisation-controlling
    cell could reach; the scalar radial case above is what a scalar-impedance cell gives."""
    j = d["amp"] * np.exp(-1j * d["phi_ref"]) * m
    return j, np.zeros_like(j)


def lp_states(m, bits):
    """Quantised scalar modulation: `(state_index, m_q)`. `bits=None`: continuous."""
    if bits is None:
        return None, m
    return quantise_nearest(m, uniform_levels(bits))


# ------------------------------------------------------------------------------------------------ pattern and metrics
def far_field(jx, jy, x, y, green, theta_deg, phi_deg, pitch_m, chunk=1500):
    """`(E_theta, E_phi)` on the points `(theta_deg, phi_deg)` (same shape) by the Phase 26 aperture sum, evaluated in chunks."""
    th, ph = np.asarray(theta_deg, dtype=float).ravel(), np.asarray(phi_deg, dtype=float).ravel()
    et = np.zeros(th.size, dtype=complex)
    ep = np.zeros(th.size, dtype=complex)
    for lo in range(0, th.size, chunk):
        sl = slice(lo, lo + chunk)
        et[sl], ep[sl] = _vector_far_field(th[sl], ph[sl], x, y, jx, jy, F, green, pitch_m, pitch_m)
    shape = np.shape(theta_deg)
    return et.reshape(shape), ep.reshape(shape)


def _angle_deg(t1, p1, t2, p2):
    t1, p1, t2, p2 = (np.radians(np.asarray(v, dtype=float)) for v in (t1, p1, t2, p2))
    c = np.cos(t1) * np.cos(t2) + np.sin(t1) * np.sin(t2) * np.cos(p1 - p2)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def refine_grid(theta_deg, phi_deg, half_width_deg=2.5, step_deg=0.25):
    """Direction-cosine grid of step `step_deg` (in sin-angle units) around `(theta, phi)`: valid at broadside too. Returns `(theta, phi, du)` 2-D."""
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    u0, v0 = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)
    s = np.sin(np.radians(step_deg))
    n = int(round(half_width_deg / step_deg))
    ug, vg = np.meshgrid(u0 + s * np.arange(-n, n + 1), v0 + s * np.arange(-n, n + 1), indexing="ij")
    r = np.hypot(ug, vg)
    ok = r < 1.0
    tg = np.degrees(np.arcsin(np.where(ok, r, 0.0)))
    pg = np.degrees(np.arctan2(vg, ug))
    return tg, pg, ok


def _main_power(et, ep, kind):
    return np.abs(copol_component(et, ep, "rhcp")) ** 2 if kind == "cp" else np.abs(et) ** 2 + np.abs(ep) ** 2


def _refined_peak(jx, jy, x, y, green, pitch_m, kind, t_start, p_start):
    tg, pg, ok = refine_grid(t_start, p_start)
    et, ep = far_field(jx, jy, x, y, green, tg, pg, pitch_m)
    pw = np.where(ok, _main_power(et, ep, kind), -1.0)
    i = np.unravel_index(np.argmax(pw), pw.shape)
    return float(pw[i]), float(tg[i]), float(pg[i]), complex(et[i]), complex(ep[i])


def analyse(jx, jy, x, y, green, pitch_m, kind, beam_deg, theta=None, phi=None):
    """Coarse hemisphere grid (1 x 2 deg, the Phase 26 grid) for the radiated power and the sidelobes, then a 0.25 deg refinement of the main beam and of
    the strongest sidelobe. `kind` 'cp' (RHCP co-pol is the main field) or 'lp' (total power). Returns a flat dict of plain numbers."""
    theta = np.arange(0.0, 90.0 + 1e-9, 1.0) if theta is None else theta
    phi = np.arange(0.0, 360.0, 2.0) if phi is None else phi
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    et, ep = far_field(jx, jy, x, y, green, tg, pg, pitch_m)
    total = np.abs(et) ** 2 + np.abs(ep) ** 2
    w = np.sin(np.radians(tg)) * np.radians(theta[1] - theta[0]) * np.radians(phi[1] - phi[0])
    p_rad = float(np.sum(total * w))
    main = _main_power(et, ep, kind)
    near = _angle_deg(tg, pg, beam_deg[0], beam_deg[1]) <= 12.0
    i = np.unravel_index(np.argmax(np.where(near, main, -1.0)), main.shape)
    peak_coarse = float(main[i])
    peak, t_pk, p_pk, e_t, e_p = _refined_peak(jx, jy, x, y, green, pitch_m, kind, tg[i], pg[i])
    d_ap = 10 * math.log10(4 * math.pi * peak / p_rad)
    d_ap_coarse = 10 * math.log10(4 * math.pi * peak_coarse / p_rad)
    dpk = _angle_deg(tg, pg, t_pk, p_pk)
    side = np.where(dpk > 15.0, main, -1.0)
    k = np.unravel_index(np.argmax(side), side.shape)
    s_peak, _, _, _, _ = _refined_peak(jx, jy, x, y, green, pitch_m, kind, tg[k], pg[k])
    s_peak = max(s_peak, float(main[k]))
    out = {
        "d_ap_dbi": d_ap,
        "d_ap_coarse_grid_dbi": d_ap_coarse,
        "peak_theta_deg": t_pk,
        "peak_phi_deg": p_pk % 360.0,
        "pointing_error_deg": float(_angle_deg(t_pk, p_pk, beam_deg[0], beam_deg[1])),
        "sll_db": 10 * math.log10(s_peak / peak),
        "sll_coarse_grid_db": 10 * math.log10(float(main[k]) / peak_coarse),
        "radiated_power_arb": p_rad,
    }
    if kind == "cp":
        rh, lh = np.abs(copol_component(et, ep, "rhcp")) ** 2, np.abs(copol_component(et, ep, "lhcp")) ** 2
        cop, cr = math.sqrt(peak), abs(copol_component(np.array(e_t), np.array(e_p), "lhcp"))
        xpd = 10 * math.log10(peak / cr**2) if cr > 0 else float("inf")
        ar = 20 * math.log10((cop + cr) / abs(cop - cr)) if abs(cop - cr) > 0 else float("inf")
        beam_pts = (rh >= 0.5 * peak_coarse) & (dpk <= 10.0)
        c, x_ = np.sqrt(rh[beam_pts]), np.sqrt(lh[beam_pts])
        ar_worst = float(np.max(20 * np.log10((c + x_) / np.maximum(np.abs(c - x_), 1e-300))))
        out["xpd_at_peak"] = bounded_db(xpd)
        out["ar_at_peak"] = bounded_db(ar)
        out["ar_worst_in_minus3db_beam"] = bounded_db(ar_worst)
        out["rhcp_fraction"] = float(np.sum(rh * w) / p_rad)
        c_c, x_c = (
            math.sqrt(float(rh[i])),
            math.sqrt(float(lh[i])),
        )  # the same point and formulas as Phase 26 `metrics` (coarse-grid peak)
        out["xpd_at_coarse_peak_db"] = float(10 * math.log10(rh[i] / max(lh[i], 1e-30)))
        out["ar_at_coarse_peak_db"] = float(20 * math.log10((c_c + x_c) / max(abs(c_c - x_c), 1e-12)))
        out["wrong_handed_at_peak"] = bool(cr > cop)
    else:
        out["etheta_over_ephi_at_peak_db"] = float(20 * math.log10(max(abs(e_t), 1e-300) / max(abs(e_p), 1e-300)))
    return out


def fixed_source_loss_db(cont, quant):
    """Gain loss of a quantised pattern against the continuous one for the SAME cell currents: the peak field scales as `D_ap * P_rad`, so the loss is
    `(D_c + 10 log10 P_c) - (D_q + 10 log10 P_q)` [dB]. Right for a PHASE-only quantisation (cp: |J| unchanged, the quantisation error falls into the non-radiating
    spectrum and is lost); `D_ap` alone renormalises that lost power away. For an AMPLITUDE quantisation (lp, lpx: |J| changes, which would change alpha) the
    renormalised `D_ap` loss is the energy-consistent one."""
    return (cont["d_ap_dbi"] + 10 * math.log10(cont["radiated_power_arb"])) - (
        quant["d_ap_dbi"] + 10 * math.log10(quant["radiated_power_arb"])
    )


def launch_gain_dbi(d_ap_dbi, eta_conv, launch=LAUNCH_BARE):
    """The Phase 26 bookkeeping: gain referred to the accepted power."""
    return d_ap_dbi + 10 * math.log10(launch * eta_conv)


def table_extrapolation_fraction(march, table):
    """Fraction of the design march's (ray, rho) points inside the aperture whose modulation wavenumber K lies outside the supercell table's K range."""
    k = march["K"]
    inside = march["rho"][None, :] <= march["r_edge"][:, None]
    out = (k < table.k_grid.min()) | (k > table.k_grid.max())
    return float(np.sum(out & inside) / np.sum(inside))


# ------------------------------------------------------------------------------------------------ practicality
def state_bit_stats(state_idx, bits):
    """Mean fraction of bits equal to 1 in the natural binary code of the state indices (bit = 1 is the forward-biased diode)."""
    s = np.asarray(state_idx, dtype=int)
    ones = sum(np.sum((s >> b) & 1) for b in range(bits))
    return float(ones / (s.size * bits))


def bias_budget(
    n_cells,
    bits,
    fraction_on,
    i_forward_a=10e-3,
    v_rail_v=3.3,
    outputs_per_chip=(8, 16),
    chain_chips=8,
    clock_hz=10e6,
    alt_currents_a=(1e-3,),
):
    """Counts and DC power of the control hardware. ONE method everywhere: `P_dc = N_forward * I_F * V_rail` (the supply power, the series resistor included);
    `N_forward = fraction_on * n_diodes`. Diodes per cell = bits (one diode per bit, an ASSUMPTION - the real topology is item 3)."""
    n_diodes = n_cells * bits
    n_fwd = fraction_on * n_diodes
    row = {
        "diodes_per_cell": bits,
        "n_cells": int(n_cells),
        "n_diodes": int(n_diodes),
        "independent_bias_lines": int(n_diodes),
        "forward_diodes_expected": float(n_fwd),
        "forward_diodes_all_on": int(n_diodes),
        "dc_current_expected_a": float(n_fwd * i_forward_a),
        "dc_power_expected_w": float(n_fwd * i_forward_a * v_rail_v),
        "dc_power_all_on_w": float(n_diodes * i_forward_a * v_rail_v),
        "driver_current_per_output_ma": float(i_forward_a * 1e3),
        "dc_at_lower_bias": {
            f"{i * 1e3:g}mA": {
                "dc_current_expected_a": float(n_fwd * i),
                "dc_power_expected_w": float(n_fwd * i * v_rail_v),
                "dc_power_all_on_w": float(n_diodes * i * v_rail_v),
            }
            for i in alt_currents_a
        },
        "bias_line_sharing": {},
        "shift_register": {},
    }
    for share in (1, 2, 4):
        row["bias_line_sharing"][f"diodes_per_line_{share}"] = int(math.ceil(n_diodes / share))
    for outs in outputs_per_chip:
        chips = int(math.ceil(n_diodes / outs))
        chains = int(math.ceil(chips / chain_chips))
        row["shift_register"][f"outputs_per_chip_{outs}"] = {
            "chips": chips,
            "chains_of_8_chips": chains,
            "fpga_io_pins": chains + 3,
            "bits_in_longest_chain": int(min(chips, chain_chips) * outs),
            "update_time_us_at_10mhz": float(min(chips, chain_chips) * outs / clock_hz * 1e6),
            "chip_current_all_on_ma": float(outs * i_forward_a * 1e3),
        }
    return row


# ------------------------------------------------------------------------------------------------ switching contrast
def _y_mean_for_beta(beta_over_k0):
    """The mean patch-sheet admittance `Y = j y` [S] whose TM0 constant is `beta_over_k0 k0` (Phase 26 sheet model, root-solved in y)."""
    k0 = free_space_wavenumber(F)
    f = lambda y: find_beta_from_sheet_admittance(1j * y * 1e-3, SUB, F).beta_rad_per_m / k0 - beta_over_k0  # noqa: E731
    return 1j * brentq(f, 0.5, 20.0, xtol=1e-9) * 1e-3


def first_harmonic_gain(levels, n=4096):
    """`g` such that the quantised cos(theta) has the first-harmonic amplitude `g` (continuous: 1, 2-level: 4/pi)."""
    th = (np.arange(n) + 0.5) * 2 * np.pi / n
    _, q = quantise_nearest(np.cos(th), levels)
    return float(2.0 * np.mean(q * np.cos(th)))


def contrast_from_alpha(alpha_np_per_m, beta_over_k0, theta0_deg, bits):
    """The switching contrast a coded cell needs to deliver the leakage `alpha_np_per_m` (a ceiling: the largest local leakage constant the coded
    surface can realise). DEFINITIONS (the stop-point report states them again):

    * TENSOR (CP, rotating axis): the cell carries a fixed anisotropy `delta` [S], `Y = Ybar I + delta R(a) diag(1,-1) R(a)^T`; the n=-1 harmonic gives
      `alpha = k_T delta^2` (`leakage_alpha_first_order_tensor`, radiated harmonic at the beam polar angle). `dX_aniso` = X(Ybar+delta) - X(Ybar-delta)
      [ohm] is the reactance difference between the two principal axes the surface wave sees; the largest state-to-state change of the sheet tensor
      (two axes 90 deg apart) is `2 delta` [S]. The switch rotates the axis; it does not change dX.
    * SCALAR (LP, levels +-1 full swing): states `Ybar +- j dB/2`; first harmonic `c1 = (dB/2) g / 2`, `alpha = k_S c1^2` (`leakage_alpha_first_order`,
      TM (zeta = 0) and TE (zeta = 90 deg) radiated parts; the larger dB is reported). `dX_switch` = X(Ybar + j dB/2) - X(Ybar - j dB/2) [ohm] with X = Im(1/(Y + Y_slab,TM(beta)))
      (`state_surface_reactances`, the impenetrable-surface reactance of the TM surface wave); `d_beta` = difference of the TM0 constants of the two UNIFORM states.
    """
    k0 = free_space_wavenumber(F)
    th = np.radians(theta0_deg)
    ybar = _y_mean_for_beta(beta_over_k0)
    beta = find_beta_from_sheet_admittance(ybar, SUB, F).beta_rad_per_m
    out = {
        "beta_over_k0_mean": beta / k0,
        "y_mean_ms": float(ybar.imag * 1e3),
        "alpha_target_np_per_m": float(alpha_np_per_m),
    }
    # tensor
    ref = 1e-3
    k_t = leakage_alpha_first_order_tensor(beta, ybar, ref, th, SUB, F) / ref**2
    delta = math.sqrt(alpha_np_per_m / k_t)
    ys = [ybar + 1j * delta, ybar - 1j * delta]  # lossless anisotropy: a susceptance difference
    xs = state_surface_reactances(ys, beta, SUB, F)
    bs = [find_beta_from_sheet_admittance(y, SUB, F).beta_rad_per_m / k0 for y in ys]
    out["tensor"] = {
        "delta_ms": delta * 1e3,
        "state_to_state_tensor_change_ms": 2 * delta * 1e3,
        "dx_aniso_ohm": float(abs(xs[0] - xs[1])),
        "d_beta_over_k0_principal_axes": float(abs(bs[0] - bs[1])),
        "modulation_depth_delta_over_ydown": delta / abs(ybar + _slab_tm(beta)),
    }
    # scalar
    g = first_harmonic_gain(uniform_levels(bits))
    worst = None
    for zeta in (0.0, 90.0):
        k_s = leakage_alpha_first_order(beta, ybar, 1.0, th, np.radians(zeta), SUB, F)  # alpha per c1^2 with c1 = 1 S
        c1 = math.sqrt(alpha_np_per_m / k_s)
        db = 4.0 * c1 / g
        if worst is None or db > worst[0]:
            worst = (db, zeta)
    db, zeta = worst
    ys = [ybar + 0.5j * db, ybar - 0.5j * db]
    xs = state_surface_reactances(ys, beta, SUB, F)
    bs = [find_beta_from_sheet_admittance(y, SUB, F).beta_rad_per_m / k0 for y in ys]
    out["scalar"] = {
        "bits": bits,
        "first_harmonic_gain_g": g,
        "worst_zeta_deg": zeta,
        "d_b_full_swing_ms": db * 1e3,
        "dx_switch_ohm": float(abs(xs[0] - xs[1])),
        "x_states_ohm": [float(xs[0]), float(xs[1])],
        "d_beta_over_k0_uniform_states": float(abs(bs[0] - bs[1])),
        "modulation_depth_db_over_ydown": db / 2 / abs(ybar + _slab_tm(beta)),
    }
    return out


def _slab_tm(beta):
    from ..dispersion.sheet_model import _y_slab
    from ..frequency import angular_frequency

    return _y_slab(beta, angular_frequency(F), SUB.eps_complex, free_space_wavenumber(F), SUB.thickness)


def capped_table(table, alpha_cap_np_per_m):
    """The supercell table with every leakage constant limited to `alpha_cap_np_per_m` (emulates a coded cell whose contrast cannot reach more)."""
    return dataclasses.replace(table, alpha=np.minimum(table.alpha, alpha_cap_np_per_m))


def green_function():
    return LayeredGreenFunction(SUB, F, backing="pec")
