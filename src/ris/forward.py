"""Phase 20 R5-2: plane-wave forward model of a coded reflecting RIS (MODEL_ONLY).

A rectangular array of `rows x cols` cells at pitch P, each cell in one of the four CP-2 states, is illuminated by a plane wave of propagation-relative RHCP
at (theta_i, phi_i). The reflected field is the physical-optics (aperture) field of the cell responses

    a_c(m, n) = C[c, R](state(m, n)) * exp(-j k0 (ux_i x_m + uy_i y_n)),   c in {R (co-handed), L (opposite-handed)},

with `C` the structure-frame circular response of `CalibratedLibrary.response` at the incident angle (every cell sees the same local incidence), `u = sin(theta)
(cos phi, sin phi)` direction cosines, and `(theta_i, phi_i)` the direction of the SPECULAR beam (the incident tangential wavevector). Each cell is a
piecewise-constant P x P patch, so the element factor sinc(k0 ux P/2) sinc(k0 uy P/2) and the grating / quantisation lobes are inherent:

    E_c(u) = obliquity(theta) * sum_mn a_c(m, n) e^{+j k0 (ux x_m + uy y_n)} * P^2 sinc(.) sinc(.)

Per-direction circular frames and the structure-frame turn exp(+-j phi) multiply all cells alike, so they change no magnitude and no relative level
(this is why the scalar channel treatment is exact for |E|, gain, SLL and axial ratio). Gain is relative to the power intercepted by the aperture,
G = 4 pi U / P_inc = 4 pi |E|^2 / (lambda^2 A cos(theta_i)) with E in units of the incident field times area; at normal incidence a uniform aperture gives
4 pi A / lambda^2. `obliquity` is the documented MODEL choice `cos(theta)` (the convention of `src/pattern/aperture_pattern.py`) or `none`.

No coupling between cells, no edge diffraction beyond the finite aperture sum, no feed: labels MODEL_ONLY; labels of the cell responses are carried through.
"""

from dataclasses import dataclass, field

import numpy as np

from .polarization import L, R

C0 = 299_792_458.0
OBLIQUITY_MODES = ("cos", "none")


def _sinc(x):
    return np.sinc(x / np.pi)  # numpy sinc is sin(pi x)/(pi x)


@dataclass(frozen=True)
class ApertureField:
    """Per-cell reflected amplitudes for unit incident RHCP."""

    a_co: np.ndarray  # (rows, cols) complex, reflected R (handedness preserved)
    a_cross: np.ndarray  # (rows, cols) complex, reflected L
    pitch_m: float
    f_ghz: float
    theta_i_deg: float
    phi_i_deg: float
    labels: tuple = ()
    cell_power: float = 0.0  # mean over cells of |co|^2 + |cross|^2 (before the incidence phase)

    @property
    def shape(self):
        return self.a_co.shape

    @property
    def area_m2(self):
        return self.a_co.size * self.pitch_m**2


def cell_positions(rows, cols, pitch_m):
    """Cell-centre coordinates (x along columns, y along rows), array centred on the origin."""
    x = (np.arange(cols) - (cols - 1) / 2.0) * pitch_m
    y = (np.arange(rows) - (rows - 1) / 2.0) * pitch_m
    return x, y


def direction_cosines(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)


def aperture_from_states(cal, states, f_ghz, theta_i_deg=0.0, phi_i_deg=0.0, pitch_m=15e-3, with_tolerance_coff=None):
    """Aperture field of a coding. `states` is a (rows, cols) integer array with values 1..4 (CP-2 one-hot states)."""
    st = np.asarray(states, dtype=int)
    if st.ndim != 2:
        raise ValueError("states must be a 2-D array")
    table, labels = {}, set()
    for s in np.unique(st):
        r = cal.response(f_ghz, theta_i_deg, phi_i_deg, int(s), with_tolerance_coff=with_tolerance_coff)
        table[int(s)] = r.c[:, R]  # reflected (R, L) for unit incident R
        labels.update(r.labels)
    co = np.vectorize(lambda s: table[int(s)][R])(st)
    cross = np.vectorize(lambda s: table[int(s)][L])(st)
    return aperture_from_amplitudes(co, cross, f_ghz, theta_i_deg, phi_i_deg, pitch_m, tuple(sorted(labels)))


def aperture_from_amplitudes(co, cross, f_ghz, theta_i_deg=0.0, phi_i_deg=0.0, pitch_m=15e-3, labels=()):
    """Aperture field from per-cell reflection amplitudes (co-handed, opposite-handed); applies the incident-wave phase."""
    co = np.asarray(co, dtype=complex)
    cross = np.broadcast_to(np.asarray(cross, dtype=complex), co.shape).copy()
    k0 = 2.0 * np.pi * f_ghz * 1e9 / C0
    x, y = cell_positions(co.shape[0], co.shape[1], pitch_m)
    uxi, uyi = direction_cosines(theta_i_deg, phi_i_deg)
    inc = np.exp(-1j * k0 * (uyi * y[:, None] + uxi * x[None, :]))
    power = float(np.mean(np.abs(co) ** 2 + np.abs(cross) ** 2))
    return ApertureField(co * inc, cross * inc, pitch_m, f_ghz, theta_i_deg, phi_i_deg, tuple(labels), power)


def _obliquity(theta, mode):
    if mode not in OBLIQUITY_MODES:
        raise ValueError(f"obliquity must be one of {OBLIQUITY_MODES}")
    return np.cos(theta) if mode == "cos" else np.ones_like(theta)


def _fields(ap, ux, uy, obliquity):
    """Complex far-field amplitudes (co, cross) at direction cosines `ux`, `uy` (same-shape arrays). Directions outside the unit disc are NaN."""
    ux, uy = np.asarray(ux, float), np.asarray(uy, float)
    shape = ux.shape
    ux, uy = ux.ravel(), uy.ravel()
    k0 = 2.0 * np.pi * ap.f_ghz * 1e9 / C0
    x, y = cell_positions(ap.shape[0], ap.shape[1], ap.pitch_m)
    inside = ux**2 + uy**2 <= 1.0 + 1e-12
    uxs, uys = np.where(inside, ux, 0.0), np.where(inside, uy, 0.0)
    ex = np.exp(1j * k0 * uxs[:, None] * x[None, :]) * _sinc(0.5 * k0 * uxs[:, None] * ap.pitch_m)  # (P, cols)
    ey = np.exp(1j * k0 * uys[:, None] * y[None, :]) * _sinc(0.5 * k0 * uys[:, None] * ap.pitch_m)  # (P, rows)
    theta = np.arcsin(np.clip(np.sqrt(uxs**2 + uys**2), 0.0, 1.0))
    ob = _obliquity(theta, obliquity) * ap.pitch_m**2
    out = []
    for a in (ap.a_co, ap.a_cross):
        e = np.einsum("pn,nm,pm->p", ey, a, ex) * ob
        out.append(np.where(inside, e, np.nan).reshape(shape))
    return out[0], out[1]


def gain_from_field(ap, e):
    """Gain (linear) of a complex far-field amplitude `e` (output of `_fields`), relative to the power intercepted by the aperture."""
    lam = C0 / (ap.f_ghz * 1e9)
    return 4.0 * np.pi * np.abs(e) ** 2 / (lam**2 * ap.area_m2 * np.cos(np.radians(ap.theta_i_deg)))


@dataclass(frozen=True)
class Pattern:
    ux: np.ndarray  # (nv, nu) grid
    uy: np.ndarray
    g_co: np.ndarray  # gain, linear, NaN outside the visible disc
    g_cross: np.ndarray
    e_co: np.ndarray
    e_cross: np.ndarray


def pattern_grid(ap, step=0.004, obliquity="cos"):
    """Gain pattern on a regular direction-cosine grid over the visible disc (separable evaluation, cheap)."""
    u = np.arange(-1.0, 1.0 + 0.5 * step, step)
    ux, uy = np.meshgrid(u, u)
    ec, ex = _fields(ap, ux, uy, obliquity)
    return Pattern(ux, uy, gain_from_field(ap, ec), gain_from_field(ap, ex), ec, ex)


def pattern_cut(ap, phi_deg, theta_deg, obliquity="cos"):
    """Co / cross gain (linear) along the cut `phi_deg` for signed-free polar angles `theta_deg` (array, 0..90)."""
    ux, uy = direction_cosines(np.asarray(theta_deg, float), phi_deg)
    ec, ex = _fields(ap, ux, uy, obliquity)
    return gain_from_field(ap, ec), gain_from_field(ap, ex)


def _db(g):
    return 10.0 * np.log10(np.maximum(g, 1e-300))


def _angles(ux, uy):
    return np.degrees(np.arcsin(np.clip(np.hypot(ux, uy), 0, 1))), np.degrees(np.arctan2(uy, ux)) % 360.0


def _main_lobe_region(g, i0, j0):
    """Pixels reachable from the peak (i0, j0) along non-increasing paths: the main lobe down to its first null / saddle."""
    mask = np.zeros(g.shape, bool)
    stack = [(i0, j0)]
    mask[i0, j0] = True
    ny, nx = g.shape
    while stack:
        i, j = stack.pop()
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                a, b = i + di, j + dj
                if (
                    (di or dj)
                    and 0 <= a < ny
                    and 0 <= b < nx
                    and not mask[a, b]
                    and np.isfinite(g[a, b])
                    and g[a, b] <= g[i, j]
                ):
                    mask[a, b] = True
                    stack.append((a, b))
    return mask


def _local_maxima(g, floor):
    pad = np.pad(g, 1, constant_values=-np.inf)
    ok = np.isfinite(g) & (g > floor)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di or dj:
                nb = pad[1 + di : 1 + di + g.shape[0], 1 + dj : 1 + dj + g.shape[1]]
                ok &= g >= np.where(np.isfinite(nb), nb, -np.inf)
    return np.argwhere(ok)


@dataclass(frozen=True)
class Metrics:
    peak_gain_dbi: float
    peak_theta_deg: float
    peak_phi_deg: float
    directivity_dbi: float
    aperture_efficiency: float  # D / (4 pi A / lambda^2)
    radiated_over_incident: float  # P_rad / P_inc from the pattern integral
    sll_db: float  # strongest secondary lobe (any lobe outside the main lobe region) re the peak; -inf if none
    sll_theta_deg: float
    sll_phi_deg: float
    cross_peak_db: float  # opposite-handed peak re the co-pol peak
    cross_peak_theta_deg: float
    cross_peak_phi_deg: float
    cross_at_specular_db: float  # opposite-handed level at the specular direction re the co-pol peak
    axial_ratio_db: float  # on the co-pol peak
    pointing_error_deg: float  # vs the requested direction (NaN when none was given)
    lobes: tuple = field(default_factory=tuple)  # strongest secondary co-pol lobes: (level_db_re_peak, theta, phi)


def metrics(ap, pat=None, target=None, obliquity="cos", step=0.004, n_lobes=6):
    """Pattern metrics. `target` = (theta, phi) of the requested beam for the pointing error."""
    pat = pat if pat is not None else pattern_grid(ap, step, obliquity)
    g = pat.g_co
    flat = np.nan_to_num(g, nan=-1.0)
    i0, j0 = np.unravel_index(np.argmax(flat), g.shape)
    ux0, uy0 = _refine_peak(ap, pat.ux[i0, j0], pat.uy[i0, j0], pat.ux[0, 1] - pat.ux[0, 0], obliquity)
    ec0, ex0 = _fields(ap, np.array([ux0]), np.array([uy0]), obliquity)
    gpk = float(gain_from_field(ap, ec0)[0])
    th0, ph0 = (float(v[0]) for v in _angles(np.array([ux0]), np.array([uy0])))
    # radiated power: integral of the (co + cross) gain over the hemisphere, d(Omega) = du dv / cos(theta)
    du = pat.ux[0, 1] - pat.ux[0, 0]
    cos_t = np.sqrt(np.clip(1.0 - pat.ux**2 - pat.uy**2, 1e-12, None))
    gt = np.nan_to_num(pat.g_co + pat.g_cross, nan=0.0)
    p_rel = float(np.sum(gt / cos_t) * du * du / (4.0 * np.pi))  # = P_rad / P_inc  (G = 4 pi U / P_inc)
    d = gpk / p_rel if p_rel > 0 else float("nan")
    lam = C0 / (ap.f_ghz * 1e9)
    eff = d / (4.0 * np.pi * ap.area_m2 / lam**2)

    region = _main_lobe_region(flat, i0, j0)
    cand = _local_maxima(np.where(region, -np.inf, g), 1e-6 * gpk)
    lobes = sorted(((float(g[i, j]), i, j) for i, j in cand), reverse=True)[:n_lobes]
    lobe_t = tuple((_db(v / gpk), *(float(a) for a in _angles(pat.ux[i, j], pat.uy[i, j]))) for v, i, j in lobes)
    sll = lobe_t[0] if lobe_t else (-np.inf, float("nan"), float("nan"))

    gx = np.nan_to_num(pat.g_cross, nan=-1.0)
    ic, jc = np.unravel_index(np.argmax(gx), gx.shape)
    tc, pc = (float(a) for a in _angles(pat.ux[ic, jc], pat.uy[ic, jc]))
    uxs, uys = direction_cosines(ap.theta_i_deg, ap.phi_i_deg)
    ecs, exs = _fields(ap, np.array([uxs]), np.array([uys]), obliquity)
    cross_spec = float(gain_from_field(ap, exs)[0])

    a, b = abs(ec0[0]), abs(ex0[0])
    ar = float("inf") if abs(a - b) < 1e-15 else 20.0 * np.log10((a + b) / abs(a - b))
    pe = float("nan")
    if target is not None:
        pe = _great_circle_deg(th0, ph0, *target)
    return Metrics(
        _db(gpk), th0, ph0, _db(d), float(eff), p_rel, sll[0], sll[1], sll[2],
        _db(gx[ic, jc] / gpk), tc, pc, _db(cross_spec / gpk), ar, pe, lobe_t,
    )  # fmt: skip


def _great_circle_deg(t1, p1, t2, p2):
    a = np.array(direction_cosines(t1, p1) + (np.cos(np.radians(t1)),))
    b = np.array(direction_cosines(t2, p2) + (np.cos(np.radians(t2)),))
    return float(np.degrees(np.arccos(np.clip(a @ b, -1.0, 1.0))))


def _refine_peak(ap, ux, uy, h, obliquity, iters=40):
    """Local refinement of the co-pol peak (pattern search around the grid maximum)."""
    best = float(np.abs(_fields(ap, np.array([ux]), np.array([uy]), obliquity)[0][0]))
    for _ in range(iters):
        moved = False
        for dx, dy in ((h, 0), (-h, 0), (0, h), (0, -h)):
            cx, cy = ux + dx, uy + dy
            if cx * cx + cy * cy > 1.0:
                continue
            v = float(np.abs(_fields(ap, np.array([cx]), np.array([cy]), obliquity)[0][0]))
            if v > best:
                best, ux, uy, moved = v, cx, cy, True
        if not moved:
            h *= 0.5
        if h < 1e-6:
            break
    return ux, uy
