"""Phase 29 design point B2: shared aperture, 7 feeds clustered near the centre, angle-dependent hologram superposition after
Frey et al., IEEE OJAP 4, 392 (2023), 10.1109/OJAP.2023.3261903 (Sections II-III), in the repo's rotating-slot tensor cell model.

Additive; reuses `shared_aperture_multiport` (marches, incident wave, cell current, port fields, scan test) unmodified. Labels of every result:
MODEL_ONLY, NOT_HFSS_VALIDATED, SOLVER_SUPERCELL, INTER_PORT_COUPLING_NOT_MODELLED.

Mapping of the paper's method onto the model (the model's cell law is `J_n = kappa R_a s_n`: the slot reflects the incident vector of feed `n`,
azimuth `psi_n`, about its axis `a`, so the current direction is `2a - psi_n` and the useful circular component carries the phase of that direction):
* ideal sub-hologram `n` (Section II-A, eq. 6-7): objective phase `chi_n = phi_SW,n - k_t,n . r + psi0_n` (hologram offset `psi0_n`, a free constant per
  feed), required slot axis `a_n = (chi_n + psi_n) / 2` and required strength `kappa_n = sqrt(2 alpha_n)` from the single-feed design march;
* angular-dependent superposition (Section II-C, eq. 12-14, 17-18): the shared slot orientation is the circular mean of the required orientations
  (doubled-angle domain, so the modulo-pi ambiguity of a slot is handled: eq. 17/18), the shared slot width the table column with the lowest summed
  strength error over the feeds, each feed seeing the table value at ITS OWN local modulation wavenumber (the model's incident-angle dependence);
* `mode='aia'` is the state-of-the-art comparison (eq. 8-9, pixel-wise arithmetic averaging): the complex unit tensors `kappa_n exp(j (chi_n + psi_n))` are
  averaged, so conflicting orientations shrink the strength instead of being spread as phase errors.
The model has no angle dependence of the strength `kappa` itself (the paper's peanut-shaped Z(phi_inc)); that part of their LUT is NOT reproduced.
"""

import numpy as np

from . import shared_aperture_multiport as sa
from . import tensor_cp_design as t
from .tensor_leaky_cp import modulation_k


def feed_ports(ring_mm, rot_deg=0.0, radial_jitter_mm=None, tangential_jitter_mm=None, tilt_deg=30.0):
    """Centre feed (broadside) + six feeds on a ring, ring azimuth opposite to each beam azimuth (+ `rot_deg`), optional per-feed jitters (GA genes).
    Returns `[(x_mm, y_mm, theta_deg, phi_deg)]` (the `shared_aperture_multiport` port format)."""
    rj = np.zeros(6) if radial_jitter_mm is None else np.asarray(radial_jitter_mm, dtype=float)
    tj = np.zeros(6) if tangential_jitter_mm is None else np.asarray(tangential_jitter_mm, dtype=float)
    out = [(0.0, 0.0, 0.0, 0.0)]
    for k in range(6):
        ph = 60.0 * k
        az = np.radians(ph + 180.0 + rot_deg)
        r = ring_mm + rj[k]
        out.append((r * np.cos(az) - tj[k] * np.sin(az), r * np.sin(az) + tj[k] * np.cos(az), tilt_deg, ph))
    return out


def design_frey(
    table, r_mm, ports, offsets_deg, w0_mm, q0, mode="frey", n_iter=2, n_psi=240, weight_exp=0.0, defocus_rad=0.0, march_fn=None
):
    """Shared-aperture design by angle-dependent superposition (module docstring). `ports`: `[(x_mm, y_mm, theta_deg, phi_deg)]` (any number of feeds),
    `offsets_deg`: hologram offset per feed. `weight_exp`: orientation weights `(|s_n| kappa_n)^weight_exp` (0 = the paper's uniform mean).
    `defocus_rad` (scalar or one per feed): extra objective phase `defocus_rad (|r| / r_mm)^2` (r from the aperture centre) added to every feed's objective: a quadratic (defocus) phase that
    broadens the beam; 0 = the paper's pencil beams.
    `march_fn(j, feed_m, beam_rad) -> march dict` (optional): replaces `march_design` as the per-feed ideal design (e.g. the Phase 26 `march_rays` with its
    per-ray conversion and taper); default None = unchanged.
    Returns the dict `shared_aperture_multiport.port_fields` / `scan_table` consume (`x_m, y_m, J, eta_conv, ...`)."""
    k0 = t.K0
    xm, ym = (a * 1e-3 for a in sa.aperture_cells(r_mm, ports))
    poly = sa.hexagon_m(r_mm)
    n_p = len(ports)
    feeds = [(p[0] * 1e-3, p[1] * 1e-3) for p in ports]
    beams = [(np.radians(p[2]), np.radians(p[3])) for p in ports]
    off = np.radians(np.asarray(offsets_deg, dtype=float))
    defoc = np.broadcast_to(np.asarray(defocus_rad, dtype=float), (n_p,))

    def window(x, y):
        return np.exp(-((x**2 + y**2) / (w0_mm * 1e-3) ** 2))

    s_list, kappa_des, psi_c, dxy, ph_des = [], [], [], [], []
    for j in range(n_p):
        m = (
            march_fn(j, feeds[j], beams[j])
            if march_fn is not None
            else sa.march_design(table, poly, feeds[j], k0, beams[j][0], beams[j][1], q0, window, n_psi)
        )
        s, alpha, _, ps, ph = sa.incident(m, feeds[j], xm, ym)
        s_list.append(s)
        ph_des.append(ph)
        kappa_des.append(np.sqrt(2.0 * alpha))
        psi_c.append(ps)
        dxy.append((xm - feeds[j][0], ym - feeds[j][1]))
    s_list, kappa_des, psi_c = np.array(s_list), np.array(kappa_des), np.array(psi_c)
    ph_cur = np.array(ph_des)
    col = np.array([table.column(modulation_k(sa.BETA0, psi_c[j], *beams[j])) for j in range(n_p)])
    slot_kappa = np.sqrt(2.0 * np.nan_to_num(table.alpha))  # (n_slots, n_k)
    ks = np.stack([slot_kappa[:, col[j]] for j in range(n_p)])  # (n_ports, n_slots, n_cells)
    ana, out_j, hist = None, None, []
    for _ in range(n_iter + 1):
        chi = np.array(
            [
                ph_cur[j]
                - k0 * np.sin(beams[j][0]) * (dxy[j][0] * np.cos(beams[j][1]) + dxy[j][1] * np.sin(beams[j][1]))
                + off[j]
                + defoc[j] * (xm**2 + ym**2) / (r_mm * 1e-3) ** 2
                for j in range(n_p)
            ]
        )
        tgt = np.exp(1j * (chi + psi_c))  # required unit tensor phasor per feed (angle 2a_n)
        w = (
            (np.sqrt(np.sum(np.abs(s_list) ** 2, axis=-1)) * kappa_des) ** weight_exp
            if weight_exp
            else np.ones_like(kappa_des)
        )
        if mode == "frey":
            z = np.sum(w * tgt, axis=0)
            axis = np.mod(0.5 * np.angle(z), np.pi)
            err = np.sum(w[:, None, :] * (ks - kappa_des[:, None, :]) ** 2, axis=0)  # (n_slots, n_cells)
            slot = np.argmin(err, axis=0)
        else:  # AiA: complex average of the unit tensors with their strengths
            z = np.sum(kappa_des * tgt, axis=0) / n_p
            axis = np.mod(0.5 * np.angle(z), np.pi)
            kmean = np.mean(ks, axis=0)  # (n_slots, n_cells)
            slot = np.argmin(np.abs(kmean - np.abs(z)[None, :]), axis=0)
        ana = [
            sa.march_analysis(table, poly, feeds[j], k0, beams[j][0], beams[j][1], np.c_[xm, ym], slot, n_psi)
            for j in range(n_p)
        ]
        out_j = [sa.incident(ana[j], feeds[j], xm, ym) for j in range(n_p)]
        s_list = np.array([o[0] for o in out_j])
        ph_cur = np.array([o[4] for o in out_j])
        hist.append(
            float(np.mean([np.mean(np.abs(np.angle(np.exp(1j * (2 * axis - psi_c[j] - chi[j]))))) for j in range(n_p)]))
        )
    cur = [sa.tensor_current(out_j[j][0], np.sqrt(2.0 * out_j[j][1]), axis) for j in range(n_p)]
    return {
        "r_mm": r_mm,
        "ports": ports,
        "x_m": xm,
        "y_m": ym,
        "slot": slot,
        "slot_L_mm": np.array([table.slots[s][0] for s in slot]),
        "slot_W_mm": np.array([table.slots[s][1] for s in slot]),
        "axis_rad": axis,
        "J": cur,
        "eta_conv": [a["eta_conv"] for a in ana],
        "march": ana,
        "kappa_des": kappa_des,
        "feeds_m": feeds,
        "beams": beams,
        "mean_abs_orientation_error_rad_history": hist,
        "mode": mode,
    }
