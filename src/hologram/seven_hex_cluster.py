"""Phase 29: seven-hexagon tensor-CP cluster (centre + 6 around), per-tile design and pair-excitation forward model.

Additive; reuses `tensor_cp_design` / `tensor_leaky_cp` unmodified. Labels of every result: MODEL_ONLY, NOT_HFSS_VALIDATED,
SOLVER_SUPERCELL, INTER_TILE_COUPLING_NOT_MODELLED. Each tile is the Phase 26 flow (rotating-slot tensor hologram, monopole at the
tile centre) on ONE shared 6.4 mm lattice; its own field is the main aperture field only (the neighbour-pickup term of
`tensor_cp_design.far_fields` is left out on purpose: no coupling between tiles). Gain is referred to the accepted power
(1 W per excited port) and anchored to the single-hexagon HFSS run exactly as `tensor_cp_design` does (`LAUNCH`, `OFFSET_DB`).
"""

import numpy as np
from scipy.ndimage import maximum_filter

from ..green_function import LayeredGreenFunction
from ..pattern.vector_aperture_pattern import copol_component
from . import tensor_cp_design as t
from .tensor_leaky_cp import march_rays

SQRT3 = np.sqrt(3.0)


def tile_centres_mm(r_mm, gap_mm=0.0):
    """Centre (0, 0) + six tiles at distance sqrt(3) R + gap along azimuths 0, 60, ... 300 deg (vertex-at-+/-y hexagons share flat
    edges facing those directions). Returns `[(x_mm, y_mm, azimuth_deg or None)]`, index 0 = centre."""
    d = SQRT3 * r_mm + gap_mm
    out = [(0.0, 0.0, None)]
    for k in range(6):
        a = 60.0 * k
        out.append((d * np.cos(np.radians(a)), d * np.sin(np.radians(a)), a))
    return out


def _hex_poly_mm(c, r_mm):
    return np.array(
        [
            [c[0] + r_mm * np.cos(np.radians(90 + 60 * i)), c[1] + r_mm * np.sin(np.radians(90 + 60 * i))]
            for i in range(6)
        ]
    )


def tile_cells_mm(c, r_mm, keepout_mm=6.5, pitch_mm=t.P_MM):
    """Shared-lattice cells (global mm, lattice through the cluster centre) inside the hexagon of circumradius `r_mm` centred at `c`."""
    n = int(np.ceil((np.hypot(*c) + r_mm) / pitch_mm)) + 2
    k = np.arange(-n, n + 1)
    xg, yg = np.meshgrid(pitch_mm * k, pitch_mm * k)
    x, y = xg.ravel(), yg.ravel()
    m = t.inside(_hex_poly_mm(c, r_mm), x, y) & (np.hypot(x - c[0], y - c[1]) >= keepout_mm)
    return x[m], y[m]


def owned_cells_mm(centres, r_mm):
    """`tile_cells_mm` for every tile with each lattice point on a shared edge given to ONE tile (the lowest index): a point exactly
    on an edge (e.g. the x = 0 lattice column between outer tiles 2 and 3) would otherwise be built twice."""
    taken, out = set(), []
    for q in centres:
        x, y = tile_cells_mm(q[:2], r_mm)
        keys = [(round(a, 6), round(b, 6)) for a, b in zip(x, y, strict=True)]
        keep = np.array([k not in taken for k in keys])
        taken.update(k for k, kp in zip(keys, keep, strict=True) if kp)
        out.append((x[keep], y[keep]))
    return out


def design_tile(table, c, r_mm, theta_deg, phi_deg, eta_design=0.9, taper_c=0.5, n_psi=720, cells=None):
    """One tile: cells + design (Phase 26 march/`cell_design`, relative to the tile's own feed). Beam azimuth is the global one.
    `cells`: optional `(x_mm, y_mm)` from `owned_cells_mm` (default: every lattice point inside the hexagon)."""
    xs, ys = tile_cells_mm(c, r_mm) if cells is None else cells
    poly = np.array(
        [
            [r_mm * 1e-3 * np.cos(np.radians(90 + 60 * i)), r_mm * 1e-3 * np.sin(np.radians(90 + 60 * i))]
            for i in range(6)
        ]
    )
    th0, ph0 = np.radians(theta_deg), np.radians(phi_deg)
    m = march_rays(table, poly, (0.0, 0.0), t.K0, th0, ph0, eta_design, taper_c=taper_c, n_psi=n_psi)
    d = t.cell_design(m, table, (xs - c[0]) * 1e-3, (ys - c[1]) * 1e-3, (0.0, 0.0), t.K0, th0, ph0, hand=+1)
    return {
        "x_mm": xs,
        "y_mm": ys,
        "design": d,
        "eta_conv": m["eta_conv"],
        "beam": (theta_deg, phi_deg),
        "centre_mm": c,
    }


def tile_field(tile, green, theta, phi):
    """Main aperture field `(E_theta, E_phi)` at the tile's real position (global frame, origin = cluster centre), scaled so that
    its radiated power is `LAUNCH * eta_conv` of 1 W accepted (no pickup term). Phase reference: the tile's own feed."""
    em = t.pattern(tile["design"], tile["x_mm"] * 1e-3, tile["y_mm"] * 1e-3, green, theta, phi)
    s = np.sqrt(t.LAUNCH * tile["eta_conv"] / t.prad(*em, theta, phi))
    return em[0] * s, em[1] * s


def rhcp_gain_dbi(et, ep, p_acc=1.0):
    g = np.abs(copol_component(et, ep, "rhcp")) ** 2 * 4 * np.pi * 10 ** (t.OFFSET_DB / 10) / p_acc
    return 10 * np.log10(np.maximum(g, 1e-30))


def _unit(theta_deg, phi_deg):
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=-1)


def angle_between(a, b):
    """Angle in degrees between two (theta, phi) directions."""
    return float(np.degrees(np.arccos(np.clip(np.dot(_unit(*a), _unit(*b)), -1, 1))))


def lobes(g_db, theta, phi, within_db=3.0, min_sep_deg=3.0):
    """RHCP lobes within `within_db` of the global peak: `[(gain_dBi, theta, phi)]`, strongest first (lobes closer than
    `min_sep_deg` merged). Sub-grid peak by a parabola in theta and phi."""
    mx = g_db == maximum_filter(g_db, size=(5, 5), mode=("nearest", "wrap"))
    pk = g_db.max()
    cand = sorted(zip(*np.nonzero(mx & (g_db >= pk - within_db)), strict=True), key=lambda k: -g_db[k])
    kept = []
    for it, ip in cand:
        th, ph = float(theta[it]), float(phi[ip])
        if any(angle_between((th, ph), (q[1], q[2])) < min_sep_deg for q in kept):
            continue
        a = g_db[max(it - 1, 0), ip], g_db[it, ip], g_db[min(it + 1, g_db.shape[0] - 1), ip]
        b = g_db[it, (ip - 1) % g_db.shape[1]], g_db[it, ip], g_db[it, (ip + 1) % g_db.shape[1]]

        def para(v):
            den = v[0] - 2 * v[1] + v[2]
            return 0.0 if abs(den) < 1e-12 else 0.5 * (v[0] - v[2]) / den

        th, ph = th + para(a) * (theta[1] - theta[0]), ph + para(b) * (phi[1] - phi[0])
        if th < 0:  # parabola through the pole: mirror to the other side of broadside
            th, ph = -th, ph + 180.0
        kept.append((float(g_db[it, ip]), th, ph % 360.0))
    return kept


def classify_pair(lb, beam_a, beam_b, frac_edge=0.15):
    """Named state of one pair pattern. `lb`: `lobes` of the combined pattern; `beam_a/b`: each tile's own beam direction.
    SINGLE_BEAM_IN_BETWEEN: one lobe within 3 dB of the peak and it lies between the two beams (off both by > `frac_edge` of
    their separation, and no more than 1.1 x separation along the way); SINGLE_BEAM_AT_ONE_TILE: one lobe but on a tile's beam;
    MULTI_LOBE: two or more lobes within 3 dB (two beams, or one beam with fringe/grating lobes)."""
    if len(lb) != 1:
        return "MULTI_LOBE"
    p = (lb[0][1], lb[0][2])
    sep = angle_between(beam_a, beam_b)
    da, db = angle_between(p, beam_a), angle_between(p, beam_b)
    if min(da, db) > frac_edge * sep and da + db <= 1.1 * sep:
        return "SINGLE_BEAM_IN_BETWEEN"
    return "SINGLE_BEAM_AT_ONE_TILE" if min(da, db) <= frac_edge * sep else "SINGLE_BEAM_OFF_AXIS"


def pair_scan(fa, fb, theta, phi, beam_a, beam_b, phases_deg=range(0, 360, 30), amp_ratio=1.0):
    """Excite tiles a, b together (amplitudes 1 : `amp_ratio`, power 1 W + `amp_ratio`^2 W accepted), relative phase of b swept.
    Returns one row per phase: lobes, class, combined peak gain referred to the total accepted power."""
    p_acc = 1.0 + amp_ratio**2
    rows = []
    for ph_deg in phases_deg:
        w = amp_ratio * np.exp(1j * np.radians(ph_deg))
        g = rhcp_gain_dbi(fa[0] + w * fb[0], fa[1] + w * fb[1], p_acc)
        lb = lobes(g, theta, phi)
        rows.append(
            {
                "phase_deg": float(ph_deg),
                "peak_gain_dbi": lb[0][0],
                "peak_theta_phi_deg": [lb[0][1], lb[0][2]],
                "n_lobes_3db": len(lb),
                "lobes": [[round(a, 2), round(b, 2), round(c, 2)] for a, b, c in lb[:4]],
                "state": classify_pair(lb, beam_a, beam_b),
            }
        )
    return rows


def green_function():
    return LayeredGreenFunction(t.SUB, t.F, backing="pec")
