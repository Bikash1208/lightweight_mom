"""Phase 20 R5-3: coding synthesis for the plane-wave CP-2 RIS (MODEL_ONLY, local-periodicity assumption).

Given a calibrated library (`CalibratedLibrary`), a frequency, the incidence and the wanted beam (theta_s, phi_s), choose one of the four one-hot states for
every cell. The design uses the ACTUAL complex co-handed responses of the states (their measured amplitude and step distortion), not the ideal 0/90/180/270
ladder. Variants (all deterministic; seeds recorded):

* `plain`   - required gradient phase psi(m, n) = -k0 (u_s - u_i) . r, state = argmax_s Re(C_RR(s) e^{-j(psi + psi0)}), global phase psi0 scanned for the best
              coherent sum at the target.
* `dither`  - the same with a seeded random offset of +-half a quantisation step added to psi (breaks the periodic staircase that creates the strongest
              quantisation lobe, at some gain cost).
* `tilt`    - the gradient is re-aimed at nearby directions (period / tilt search, +-`tilt_u` in direction cosine) and the coding with the highest gain at the
              REQUESTED direction is kept (a period optimisation: the staircase period changes, the lobe positions move).
* `refined` - coordinate ascent on the true pattern objective from the plain coding: gain at the target minus hinge penalties on the strongest co-pol sidelobe
              and the opposite-handed level (the far field is linear in the cells, so a state change updates all directions in O(directions)).

Assumption flagged on every result: the amplitude each cell sends toward the beam is its response at the INCIDENT angle (no bistatic cell pattern);
the array-level supercell runs of R5-4 test exactly this. Sensitivity band: the same coding evaluated with C_off = 0.0225 / 0.025 / 0.030 pF
(`ToleranceBand`); no optimisation against the corners.
"""

from dataclasses import dataclass, field

import numpy as np

from . import forward as fw
from .polarization import R

STATES = (1, 2, 3, 4)
ASSUMPTIONS = ("MODEL_ONLY", "LOCAL_RESPONSE_AT_INCIDENT_ANGLE_ASSUMED", "NOT_HFSS_VALIDATED")


@dataclass(frozen=True)
class Coding:
    states: np.ndarray
    variant: str
    f_ghz: float
    theta_s: float
    phi_s: float
    theta_i: float
    phi_i: float
    seed: int | None = None
    notes: tuple = ()
    log: dict = field(default_factory=dict)


def required_phase(rows, cols, pitch_m, f_ghz, theta_s, phi_s, theta_i=0.0, phi_i=0.0):
    """Phase each cell must add (arg of its reflection coefficient) so the reflected field of the incident plane wave adds up toward (theta_s, phi_s)."""
    k0 = 2.0 * np.pi * f_ghz * 1e9 / fw.C0
    x, y = fw.cell_positions(rows, cols, pitch_m)
    uxs, uys = fw.direction_cosines(theta_s, phi_s)
    uxi, uyi = fw.direction_cosines(theta_i, phi_i)
    return -k0 * ((uxs - uxi) * x[None, :] + (uys - uyi) * y[:, None])


def _state_table(cal, f_ghz, theta_i, phi_i, coff=None):
    c = np.array(
        [cal.response(f_ghz, theta_i, phi_i, s, with_tolerance_coff=coff).c[:, R] for s in STATES]
    )  # (4, 2): reflected (R, L)
    return c[:, 0], c[:, 1]


def nearest_states(co_table, psi, offsets=72):
    """argmax_s Re(C_s e^{-j(psi + psi0)}) per cell with the global offset psi0 (scan of `offsets` values) that maximises the total coherent sum."""
    best, best_val = None, -np.inf
    for psi0 in np.linspace(0.0, 2 * np.pi, offsets, endpoint=False):
        proj = np.real(co_table[:, None, None] * np.exp(-1j * (psi + psi0))[None])
        pick = np.argmax(proj, axis=0)
        val = float(np.sum(np.take_along_axis(proj, pick[None], 0)))
        if val > best_val:
            best, best_val = pick, val
    return best + 1


def _gain_at(cal, states, f_ghz, theta_i, phi_i, pitch_m, target, coff=None, obliquity="cos"):
    ap = fw.aperture_from_states(cal, states, f_ghz, theta_i, phi_i, pitch_m, with_tolerance_coff=coff)
    ux, uy = fw.direction_cosines(*target)
    ec, _ = fw._fields(ap, np.array([ux]), np.array([uy]), obliquity)
    return float(fw.gain_from_field(ap, ec)[0])


def design(
    cal,
    variant,
    f_ghz,
    theta_s,
    phi_s,
    theta_i=0.0,
    phi_i=0.0,
    rows=20,
    cols=20,
    pitch_m=15e-3,
    seed=0,
    tilt_u=0.03,
    **refine_kw,
):
    table, _ = _state_table(cal, f_ghz, theta_i, phi_i)
    psi = required_phase(rows, cols, pitch_m, f_ghz, theta_s, phi_s, theta_i, phi_i)
    kw = dict(f_ghz=f_ghz, theta_s=theta_s, phi_s=phi_s, theta_i=theta_i, phi_i=phi_i)
    if variant == "plain":
        return Coding(nearest_states(table, psi), "plain", **kw, notes=ASSUMPTIONS)
    if variant == "dither":
        rng = np.random.default_rng(seed)
        d = rng.uniform(-np.pi / 4, np.pi / 4, psi.shape)  # half a quantisation step (2-bit)
        return Coding(nearest_states(table, psi + d), "dither", **kw, seed=seed, notes=ASSUMPTIONS)
    if variant == "tilt":
        uxs, uys = fw.direction_cosines(theta_s, phi_s)
        best, best_g, best_off = None, -1.0, (0.0, 0.0)
        for dx in np.linspace(-tilt_u, tilt_u, 9):
            for dy in np.linspace(-tilt_u, tilt_u, 9):
                ux, uy = uxs + dx, uys + dy
                if ux * ux + uy * uy >= 1.0:
                    continue
                t2 = np.degrees(np.arcsin(np.hypot(ux, uy)))
                p2 = np.degrees(np.arctan2(uy, ux))
                st = nearest_states(
                    table, required_phase(rows, cols, pitch_m, f_ghz, t2, p2, theta_i, phi_i), offsets=72
                )
                g = _gain_at(cal, st, f_ghz, theta_i, phi_i, pitch_m, (theta_s, phi_s))
                if g > best_g:
                    best, best_g, best_off = st, g, (float(dx), float(dy))
        return Coding(best, "tilt", **kw, notes=ASSUMPTIONS, log={"tilt_u": best_off})
    if variant == "refined":
        base = design(cal, "plain", f_ghz, theta_s, phi_s, theta_i, phi_i, rows, cols, pitch_m)
        st, log = refine(cal, base.states, f_ghz, (theta_s, phi_s), theta_i, phi_i, pitch_m, seed=seed, **refine_kw)
        return Coding(st, "refined", **kw, seed=seed, notes=ASSUMPTIONS, log=log)
    raise ValueError(f"unknown variant {variant!r}")


def refine(
    cal,
    states,
    f_ghz,
    target,
    theta_i=0.0,
    phi_i=0.0,
    pitch_m=15e-3,
    seed=0,
    sweeps=4,
    sll_limit_db=-12.0,
    cross_limit_db=-15.0,
    du=0.02,
    weight=1.0,
):
    """Coordinate ascent: J = G_target[dB] - weight * (max(0, SLL - sll_limit) + max(0, cross - cross_limit)), SLL = strongest co-pol level outside the main
    lobe (radius 2 lambda / D in direction cosine) re the target, cross = strongest opposite-handed level re the target. Visits the cells in a seeded random
    order each sweep, takes the best of the other three states, stops when a sweep changes nothing. Returns (states, log)."""
    st = np.array(states, dtype=int)
    rows, cols = st.shape
    co_t, cx_t = _state_table(cal, f_ghz, theta_i, phi_i)
    k0 = 2.0 * np.pi * f_ghz * 1e9 / fw.C0
    lam = fw.C0 / (f_ghz * 1e9)
    x, y = fw.cell_positions(rows, cols, pitch_m)
    uxi, uyi = fw.direction_cosines(theta_i, phi_i)
    inc = np.exp(-1j * k0 * (uyi * y[:, None] + uxi * x[None, :])).ravel()
    u = np.arange(-1.0, 1.0 + 0.5 * du, du)
    ux, uy = np.meshgrid(u, u)
    vis = ux**2 + uy**2 <= 1.0
    ux, uy = np.append(ux[vis], fw.direction_cosines(*target)[0]), np.append(uy[vis], fw.direction_cosines(*target)[1])
    nt = len(ux) - 1  # index of the exact target direction
    xs, ys = np.meshgrid(x, y)
    xs, ys = xs.ravel(), ys.ravel()
    th = np.arcsin(np.clip(np.hypot(ux, uy), 0, 1))
    basis = (
        np.exp(1j * k0 * (ux[:, None] * xs[None] + uy[:, None] * ys[None]))
        * np.sinc(0.5 * k0 * ux[:, None] * pitch_m / np.pi)
        * np.sinc(0.5 * k0 * uy[:, None] * pitch_m / np.pi)
        * (np.cos(th) * pitch_m**2)[:, None]
        * inc[None]
    )  # (directions, cells)
    far = np.hypot(ux - ux[nt], uy - uy[nt]) > 2.0 * lam / (max(rows, cols) * pitch_m)
    flat = st.ravel().copy()
    e_co = basis @ co_t[flat - 1]
    e_cx = basis @ cx_t[flat - 1]
    scale = 4 * np.pi / (lam**2 * rows * cols * pitch_m**2 * np.cos(np.radians(theta_i)))

    def score(ec, ex):
        pt = np.abs(ec[nt]) ** 2
        if pt <= 0:
            return -1e9
        sll = 10 * np.log10(max(np.max(np.abs(ec[far]) ** 2), 1e-300) / pt)
        cr = 10 * np.log10(max(np.max(np.abs(ex) ** 2), 1e-300) / pt)
        return 10 * np.log10(pt * scale) - weight * (max(0.0, sll - sll_limit_db) + max(0.0, cr - cross_limit_db))

    j0 = score(e_co, e_cx)
    rng = np.random.default_rng(seed)
    cur, flips = j0, 0
    for _ in range(sweeps):
        changed = False
        for c in rng.permutation(flat.size):
            s0 = flat[c]
            best_s, best_j, best_e = s0, cur, None
            for s in (1, 2, 3, 4):
                if s == s0:
                    continue
                ec = e_co + basis[:, c] * (co_t[s - 1] - co_t[s0 - 1])
                ex = e_cx + basis[:, c] * (cx_t[s - 1] - cx_t[s0 - 1])
                j = score(ec, ex)
                if j > best_j + 1e-12:
                    best_s, best_j, best_e = s, j, (ec, ex)
            if best_s != s0:
                flat[c], cur, (e_co, e_cx) = best_s, best_j, best_e
                changed, flips = True, flips + 1
        if not changed:
            break
    return flat.reshape(st.shape), {
        "objective_before": float(j0),
        "objective_after": float(cur),
        "flips": flips,
        "sll_limit_db": sll_limit_db,
        "cross_limit_db": cross_limit_db,
    }


def evaluate(cal, coding, f_ghz=None, pitch_m=15e-3, coff=None, step=0.004, obliquity="cos"):
    """Forward-model metrics of a coding at `f_ghz` (default the design frequency), optionally with the diode C_off (pF) of the tolerance band.
    Returns (Metrics, gain at the requested direction in dBi)."""
    f = coding.f_ghz if f_ghz is None else f_ghz
    ap = fw.aperture_from_states(cal, coding.states, f, coding.theta_i, coding.phi_i, pitch_m, with_tolerance_coff=coff)
    m = fw.metrics(ap, target=(coding.theta_s, coding.phi_s), obliquity=obliquity, step=step)
    ux, uy = fw.direction_cosines(coding.theta_s, coding.phi_s)
    ec, _ = fw._fields(ap, np.array([ux]), np.array([uy]), obliquity)
    g_t = float(10 * np.log10(max(fw.gain_from_field(ap, ec)[0], 1e-300)))
    return m, g_t


def ideal_references(theta_s, phi_s, f_ghz, theta_i=0.0, phi_i=0.0, rows=20, cols=20, pitch_m=15e-3, step=0.004):
    """Metrics of the ideal lossless unit-amplitude apertures with the same required phases: continuous phase and exact 2-bit ladder (nearest of 0/90/180/270
    after the best global offset). The cos(theta) obliquity convention is used as for the real codings."""
    psi = required_phase(rows, cols, pitch_m, f_ghz, theta_s, phi_s, theta_i, phi_i)
    out = {}
    cont = fw.aperture_from_amplitudes(np.exp(1j * psi), np.zeros_like(psi), f_ghz, theta_i, phi_i, pitch_m)
    out["ideal_continuous"] = fw.metrics(cont, target=(theta_s, phi_s), step=step)
    best = None
    for psi0 in np.linspace(0.0, np.pi / 2, 18, endpoint=False):
        q = np.round((psi + psi0) / (np.pi / 2)) * (np.pi / 2)
        ap = fw.aperture_from_amplitudes(np.exp(1j * q), np.zeros_like(psi), f_ghz, theta_i, phi_i, pitch_m)
        ux, uy = fw.direction_cosines(theta_s, phi_s)
        ec, _ = fw._fields(ap, np.array([ux]), np.array([uy]), "cos")
        g = float(fw.gain_from_field(ap, ec)[0])
        if best is None or g > best[0]:
            best = (g, ap)
    out["ideal_2bit"] = fw.metrics(best[1], target=(theta_s, phi_s), step=step)
    return out


def grating_onset_deg(f_ghz, pitch_m=15e-3):
    """Smallest beam angle at which a first grating lobe of the cell lattice enters the visible region (1-D: sin theta_s = lambda/P - 1)."""
    lam = fw.C0 / (f_ghz * 1e9)
    s = lam / pitch_m - 1.0
    return 90.0 if s >= 1.0 else float(np.degrees(np.arcsin(max(s, 0.0))))


def gain_bandwidth(freqs, gain_db, drops=(1.0, 3.0)):
    """{drop_dB: (f_lo, f_hi, width)}: contiguous band around the maximum of the sweep where the gain is within `drop` dB of that maximum (edges by linear
    interpolation in dB; None edge = the band reaches the end of the sweep, reported with `open_low` / `open_high`)."""
    f, g = np.asarray(freqs, float), np.asarray(gain_db, float)
    k = int(np.argmax(g))
    out = {}
    for d in drops:
        thr = g[k] - d
        lo = k
        while lo > 0 and g[lo - 1] >= thr:
            lo -= 1
        hi = k
        while hi < len(g) - 1 and g[hi + 1] >= thr:
            hi += 1
        flo = f[lo] if lo == 0 else f[lo - 1] + (thr - g[lo - 1]) * (f[lo] - f[lo - 1]) / (g[lo] - g[lo - 1])
        fhi = f[hi] if hi == len(g) - 1 else f[hi] + (thr - g[hi]) * (f[hi + 1] - f[hi]) / (g[hi + 1] - g[hi])
        out[d] = {
            "f_lo": float(flo),
            "f_hi": float(fhi),
            "width": float(fhi - flo),
            "open_low": lo == 0,
            "open_high": hi == len(g) - 1,
        }
    return out
