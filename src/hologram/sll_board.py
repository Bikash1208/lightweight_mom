"""Phase 23 C/D glue: a `/synthesize_ports` job result <-> the sidelobe optimiser, and the held-out
evaluation that every reported SLL / directivity number comes from.

* `build_problem(result, ...)` turns a coded 1-bit multi-port result (rectangles or hexagons) into the
  optimiser's `OptPort`s: every cell's two possible aperture values `s_n(state) = E_ref,n * exp(j psi_state)`
  (`E_ref` recovered from the stored `aperture_field`, `psi` the exact achieved phases in
  `result["two_state"]["achieved_phase_deg"]`), the cell-wise co-/cross-pol contributions, the requested
  beams and the scenarios (each port alone; all ports together).
* `evaluate(problem, levels, alpha)` evaluates every scenario on a grid the optimiser never saw
  (theta 0.5 deg x phi 1 deg over the hemisphere by default) with the project's own vector far field,
  and reports directivity, pointing, co-pol SLL (global, and on the board x / y cuts through each beam),
  XPD and level balance. `at_frequency` re-evaluates the SAME states at another frequency (nothing is
  re-optimised there).
* `apply_levels(result, problem, levels, alpha)` returns a new job-result dict carrying the optimised
  states (patches, aperture fields, combined grid, per-port report) so the existing exports work on it.

SLL definitions (documented because the repo's "Phase 18" metric text could not be found): the global SLL
of a beam is `beam_polarization_metrics`' (highest co-pol level outside the beam's first-null cone and
nearer to it than to any other requested beam), relative to the beam's own fine-resolution peak. The
x / y cut SLL is the same quantity along the line of directions through the beam parallel to the board x
axis (v = v0) / y axis (u = u0) in the u-v plane. This is the Phase 18 cut definition (constant v / constant
u through the beam peak: `sll_cut_along_x_db` / `sll_cut_along_y_db` in `scripts/adhoc/coded4.py` on branch
`phase18-four-port-300mm`), with the main lobe excluded by the first-null cone.

Frequency behaviour: `build_problem(extra_frequencies_hz=...)` RE-SOLVES the two cells at each frequency
(Phase 22 `auto_states`, same square side) and moves the reference wave with the TM0 phase constant.
"""

from dataclasses import dataclass, field

import numpy as np

from ..frequency import free_space_wavelength
from ..green_function import LayeredGreenFunction
from ..materials import Substrate
from ..pattern.cellwise_pattern import (
    cell_contributions,
    uv_direction_set,
)
from ..pattern.vector_aperture_pattern import (
    ORTHOGONAL_POLARIZATION,
    _vector_far_field,
    aperture_current_from_scalar,
    copol_component,
    first_null_half_width_deg,
)
from .coded_synthesis import FEED_KEEPOUT_FRACTION
from .multi_port_synthesis import region_cell_mask
from .ports import PortRegion
from .sll_optimizer import OptBeam, OptimiserConfig, OptPort, Scenario, optimise_sidelobes

_C0 = 299792458.0


@dataclass
class BoardPort:
    port_id: int
    region: PortRegion
    x_m: np.ndarray  # (ny, nx) the port's design grid
    y_m: np.ndarray
    feed_xy_m: tuple
    owned: np.ndarray  # (ny, nx) bool: cells that belong to this port
    e_ref: np.ndarray  # (ny, nx) complex reference-wave value per cell
    levels0: np.ndarray  # (ny, nx) int starting state (0 outside `owned`)
    cell_ij: np.ndarray  # (N, 2) grid indices of the owned cells (row-major order)
    free: np.ndarray  # (N,) bool


@dataclass
class BoardProblem:
    ports: list
    beams: list  # (theta, phi, polarization) in beam-id order
    beam_port: list  # port id of each beam
    psi_rad: tuple  # achieved phase of state 0 / state 1
    frequency_hz: float
    substrate: Substrate
    period_m: float
    polarization: str
    state_geometry: tuple  # geometry dict of state 0 / state 1
    opt_ports: list = field(default_factory=list)
    opt_beams: list = field(default_factory=list)
    scenarios: list = field(default_factory=list)
    directions: object = None


def regions_from_result(result):
    """The `PortRegion`s of a job result (hexagon outlines when present)."""
    regions = []
    for p in result["ports"]:
        x0, x1, y0, y1 = (p[k] * 1e-3 for k in ("x_min_mm", "x_max_mm", "y_min_mm", "y_max_mm"))
        verts = None
        if p.get("outline_mm"):
            verts = tuple((x * 1e-3, y * 1e-3) for x, y in p["outline_mm"])
        regions.append(
            PortRegion(p["port_id"], x0, x1, y0, y1, p["feed_x_mm"] * 1e-3, p["feed_y_mm"] * 1e-3, vertices_m=verts)
        )
    return regions


def _extent_m(ports, period_m):
    xs = np.concatenate([p.x_m[p.owned] for p in ports])
    ys = np.concatenate([p.y_m[p.owned] for p in ports])
    return float(max(xs.max() - xs.min(), ys.max() - ys.min()) + period_m)


def _dbeta(problem, frequency_hz):
    """Change of the TM0 surface-wave phase constant from the design frequency to `frequency_hz`."""
    if abs(frequency_hz - problem.frequency_hz) <= 1e-3:
        return 0.0
    from ..dispersion.analytic_tm0 import tm0_beta

    h = problem.substrate.thickness
    return tm0_beta(frequency_hz, problem.substrate.eps_r, h) - tm0_beta(
        problem.frequency_hz, problem.substrate.eps_r, h
    )


def build_problem(
    result,
    frequency_hz,
    eps_r,
    tan_delta,
    thickness_mm,
    period_mm,
    port_beams,
    polarization="rhcp",
    amplitudes=None,
    du_far=0.025,
    du_near=0.012,
    dtype=np.complex64,
    extra_frequencies_hz=(),
    states_cache_dir=None,
    n_phase_draws=0,
    phase_error_deg=10.0,
    draw_seed=0,
    progress=None,
):
    """`port_beams`: `[[(theta, phi), ...], ...]` the requested beams of each port, in port order.

    `extra_frequencies_hz`: more frequencies to model (frequency checks, robust optimisation). At each,
    the two cells' reflection coefficients are RE-SOLVED with the Phase 22 automatic-states machinery
    (the same square side as the design's state 1, the empty slab's exact value) and the reference wave
    is moved with the TM0 phase constant; the design frequency is index 0. `n_phase_draws` seeded draws
    of independent port-phase errors (uniform within +-`phase_error_deg`; the first port is the phase
    reference) are kept for the all-ports scenario."""
    two_state = result.get("two_state")
    if not two_state or "achieved_phase_deg" not in two_state:
        raise ValueError("the sidelobe optimiser needs a coded 1-bit multi-port design result")
    psi = tuple(float(np.radians(v)) for v in two_state["achieved_phase_deg"])
    period_m = period_mm * 1e-3
    regions = regions_from_result(result)
    ports = []
    for p, region in zip(result["ports"], regions, strict=True):
        x_m, y_m = np.asarray(p["x_m"]), np.asarray(p["y_m"])
        earlier = [r for r in regions if r.port_id < region.port_id]
        owned = region_cell_mask(region, x_m, y_m, earlier) if region.is_polygon else np.ones(x_m.shape, dtype=bool)
        field_ap = np.asarray(p["aperture_field_real"]) + 1j * np.asarray(p["aperture_field_imag"])
        levels = np.zeros(x_m.shape, dtype=int)
        for c in p["patches_mm"]["cells"]:
            levels[c["i"], c["j"]] = 1
        levels = np.where(owned, levels, 0)
        psi_cell = np.where(levels == 1, psi[1], psi[0])
        e_ref = np.where(owned, field_ap / np.exp(1j * psi_cell), 0.0)
        ij = np.argwhere(owned)
        fx, fy = p["feed_x_mm"] * 1e-3, p["feed_y_mm"] * 1e-3
        limit = FEED_KEEPOUT_FRACTION * period_m
        under_feed = (np.abs(x_m[owned] - fx) < limit) & (np.abs(y_m[owned] - fy) < limit)
        ports.append(BoardPort(p["port_id"], region, x_m, y_m, (fx, fy), owned, e_ref, levels, ij, ~under_feed))
    geometry = tuple(s["geometry"] for s in two_state["states"])
    beams, beam_port = [], []
    for pid, plist in enumerate(port_beams):
        for th, ph in plist:
            beams.append((float(th), float(ph), polarization))
            beam_port.append(pid)
    prob = BoardProblem(
        ports,
        beams,
        beam_port,
        psi,
        frequency_hz,
        Substrate("user", eps_r, tan_delta, thickness_mm * 1e-3),
        period_m,
        polarization,
        geometry,
    )
    prob.frequencies_hz = (frequency_hz, *[float(f) for f in extra_frequencies_hz])
    prob.freq_psi = {frequency_hz: psi}
    if len(prob.frequencies_hz) > 1:
        from .auto_states import auto_states

        g1 = geometry[1]
        if g1.get("kind") != "rectangle":
            raise ValueError("per-frequency re-solving needs a square state 1 (rectangle geometry)")
        for f in prob.frequencies_hz[1:]:
            if progress:
                progress(f"re-solving the two cell states at {f / 1e9:.3f} GHz")
            a = auto_states(1, prob.substrate, f, period_m, fixed_L_mm=float(g1["L_mm"]), cache_dir=states_cache_dir)
            prob.freq_psi[f] = tuple(float(np.angle(st.gamma)) for st in a.states)
    lam = free_space_wavelength(frequency_hz)
    prob.directions = uv_direction_set(
        du_far=du_far,
        beams_deg=[(t, p) for t, p, _ in beams],
        du_near=du_near,
        near_radius_deg=20.0,
        patch_radius_deg=3.5,
        patch_step=0.004,
    )
    patch_idx = np.nonzero(prob.directions.tag >= 0)[0]
    cross_pol = ORTHOGONAL_POLARIZATION[polarization]
    from ..pattern.cellwise_pattern import DirectionSet

    patch_dirs = DirectionSet(
        *(getattr(prob.directions, k)[patch_idx] for k in ("theta_deg", "phi_deg", "u", "v", "tag"))
    )
    greens = {f: LayeredGreenFunction(prob.substrate, f, backing="pec") for f in prob.frequencies_hz}
    for bp in ports:
        x = bp.x_m[bp.owned]
        y = bp.y_m[bp.owned]
        rho = np.hypot(x - bp.feed_xy_m[0], y - bp.feed_xy_m[1])
        c_cos, c_xs, s_list = [], [], []
        for f in prob.frequencies_hz:
            c_cos.append(
                cell_contributions(
                    x, y, bp.feed_xy_m, prob.directions, f, greens[f], period_m, period_m, polarization, dtype=dtype
                )
            )
            c_xs.append(
                cell_contributions(
                    x, y, bp.feed_xy_m, patch_dirs, f, greens[f], period_m, period_m, cross_pol, dtype=dtype
                )
            )
            e = bp.e_ref[bp.owned] * np.exp(-1j * _dbeta(prob, f) * rho)
            psi_f = prob.freq_psi[f]
            s_list.append(np.stack([e * np.exp(1j * psi_f[0]), e * np.exp(1j * psi_f[1])], axis=1))
        amp = 1.0 if amplitudes is None else float(amplitudes[bp.port_id])
        prob.opt_ports.append(
            OptPort(bp.port_id, c_cos, c_xs, s_list, bp.levels0[bp.owned].astype(int), bp.free.copy(), amp)
        )
    # beams / scenarios (a beam's main-lobe cone follows the extent of the sub-aperture that radiates it)
    own_extent = {bp.port_id: _extent_m([bp], period_m) for bp in ports}
    for k, (th, ph, _pol) in enumerate(beams):
        prob.opt_beams.append(
            OptBeam(k, beam_port[k], th, ph, first_null_half_width_deg(own_extent[beam_port[k]], lam, th))
        )
    rng = np.random.default_rng(draw_seed)
    prob.draws = [
        {
            bp.port_id: (
                0.0 if bp.port_id == ports[0].port_id else float(rng.uniform(-phase_error_deg, phase_error_deg))
            )
            for bp in ports
        }
        for _ in range(n_phase_draws)
    ]
    nf = len(prob.frequencies_hz)

    def suffix(fi):
        return "" if fi == 0 else f"@{prob.frequencies_hz[fi] / 1e9:.2f}GHz"

    for bp in ports:
        ids = tuple(k for k, pid in enumerate(beam_port) if pid == bp.port_id)
        cones = tuple(first_null_half_width_deg(own_extent[bp.port_id], lam, beams[k][0]) for k in ids)
        for fi in range(nf):
            prob.scenarios.append(Scenario(f"port{bp.port_id + 1}{suffix(fi)}", (bp.port_id,), ids, 1.0, cones, fi))
    if len(ports) > 1:
        ids = tuple(range(len(beams)))
        cones = tuple(first_null_half_width_deg(own_extent[beam_port[k]], lam, beams[k][0]) for k in ids)
        all_ports = tuple(bp.port_id for bp in ports)
        for fi in range(nf):
            prob.scenarios.append(Scenario(f"all{suffix(fi)}", all_ports, ids, 1.0, cones, fi))
        for j, errs in enumerate(prob.draws):
            prob.scenarios.append(
                Scenario(f"all@draw{j}", all_ports, ids, 1.0, cones, j % nf, tuple(errs.items()), False)
            )
    return prob


def optimise(problem, mode, config=None, alpha_deg=None, robust=False):
    """`mode`: 'single' (each port alone, independent per-port problems), 'all' (all ports together,
    cells + port phases), 'both' (weighted worst case over the single-port scenarios and the all-ports
    one). `robust=True` minimises the worst case over every modelled frequency (cells re-solved there)
    and, for 'all'/'both', over the seeded port-phase-error draws. Returns `(levels_by_port, alpha_deg,
    OptimiserResult list)`; `levels_by_port` maps port id -> (N,) int in `OptPort` order."""
    cfg = config or OptimiserConfig()
    levels = {p.port_id: np.array(p.state, copy=True) for p in problem.opt_ports}
    alpha = {p.port_id: 0.0 for p in problem.opt_ports}
    if alpha_deg:
        alpha.update(alpha_deg)

    def pick(pred):
        out = [s for s in problem.scenarios if pred(s)]
        if not robust:
            out = [s for s in out if s.freq_index == 0 and s.constrain]
        else:
            # nominal instances first (they rank the candidates), then the phase-error draws
            out = sorted(out, key=lambda s: (not s.constrain, s.freq_index))
        return out

    results = []
    if mode == "single":
        for p in problem.opt_ports:
            scs = pick(lambda s, pid=p.port_id: s.ports == (pid,))
            sub_cfg = OptimiserConfig(**{**cfg.__dict__, "optimise_alpha": False})
            res = optimise_sidelobes([p], problem.opt_beams, scs, problem.directions, sub_cfg)
            levels[p.port_id] = res.state[p.port_id]
            results.append(res)
    elif mode in ("all", "both"):
        scs = pick(lambda s: len(s.ports) > 1) if mode == "all" else pick(lambda s: True)
        res = optimise_sidelobes(problem.opt_ports, problem.opt_beams, scs, problem.directions, cfg, alpha)
        levels = res.state
        alpha = res.alpha_deg
        results.append(res)
    else:
        raise ValueError(f"mode must be 'single', 'all' or 'both', got {mode!r}")
    return levels, alpha, results


def optimise_alpha_only(problem, config=None):
    """Port phases only (cells untouched) -- the fair baseline for the all-ports mode."""
    cfg = OptimiserConfig(
        **{**(config or OptimiserConfig()).__dict__, "max_flips": 0, "optimise_alpha": True, "anneal_steps": 0}
    )
    sc = [s for s in problem.scenarios if s.name == "all"]
    res = optimise_sidelobes(problem.opt_ports, problem.opt_beams, sc, problem.directions, cfg)
    return res.alpha_deg


# ------------------------------------------------------------------------------- held-out evaluation


@dataclass
class HeldOutGrid:
    theta_step_deg: float = 0.5
    phi_step_deg: float = 1.0
    cut_step_u: float = 0.0015
    patch_radius_deg: float = 3.0
    patch_step_deg: float = 0.1


def _aperture(problem, bp: BoardPort, levels_flat, dbeta_rad_per_m=0.0, psi=None):
    """Scalar aperture field of port `bp` (zero outside its cells) for per-cell `levels_flat`.
    `dbeta_rad_per_m` moves the reference wave to another frequency (`exp(-j dbeta rho)`); `psi` is
    the (state 0, state 1) phase pair at that frequency (default: the design frequency's)."""
    psi = problem.psi_rad if psi is None else psi
    lev = np.zeros(bp.x_m.shape, dtype=int)
    lev[bp.owned] = levels_flat
    psi_cell = np.where(lev == 1, psi[1], psi[0])
    e_ref = bp.e_ref
    if dbeta_rad_per_m:
        rho = np.hypot(bp.x_m - bp.feed_xy_m[0], bp.y_m - bp.feed_xy_m[1])
        e_ref = e_ref * np.exp(-1j * dbeta_rad_per_m * rho)
    return np.where(bp.owned, e_ref * np.exp(1j * psi_cell), 0.0)


def _psi_at(problem, frequency_hz):
    """The cells' phase pair at `frequency_hz`: re-solved values when the problem carries them for that
    frequency, else the design-frequency pair (the caller should know which)."""
    for f, psi in problem.freq_psi.items():
        if abs(f - frequency_hz) <= 1.0:
            return psi
    return problem.psi_rad


class _FarField:
    """E_theta, E_phi of the board for a per-port complex excitation (memoised per point set)."""

    def __init__(self, problem, levels, alpha_deg, frequency_hz, active, amplitudes=None, weights=None):
        self.problem = problem
        self.f = frequency_hz
        self.green = LayeredGreenFunction(problem.substrate, frequency_hz, backing="pec")
        self.parts = []
        self._memo = {}
        dbeta = _dbeta(problem, frequency_hz)
        psi = _psi_at(problem, frequency_hz)
        for bp in problem.ports:
            if bp.port_id not in active:
                continue
            if weights is not None:
                w = weights[bp.port_id]
            elif len(active) > 1:
                a = 1.0 if amplitudes is None else amplitudes[bp.port_id]
                w = a * np.exp(1j * np.radians(alpha_deg[bp.port_id]))
            else:
                w = 1.0
            s = _aperture(problem, bp, levels[bp.port_id], dbeta, psi) * w
            jx, jy = aperture_current_from_scalar(bp.x_m, bp.y_m, s, feed_xy_m=bp.feed_xy_m)
            sel = bp.owned
            self.parts.append((bp.x_m[sel], bp.y_m[sel], jx[sel], jy[sel]))

    def __call__(self, theta_deg, phi_deg, chunk=3000):
        theta = np.asarray(theta_deg, dtype=float).ravel()
        phi = np.asarray(phi_deg, dtype=float).ravel()
        key = (theta.size, float(theta[0]), float(theta[-1]), float(theta.sum()), float(phi.sum()))
        if key in self._memo:
            return self._memo[key]
        et = np.zeros(theta.size, dtype=complex)
        ep = np.zeros(theta.size, dtype=complex)
        dx = self.problem.period_m
        for x, y, jx, jy in self.parts:
            for lo in range(0, theta.size, chunk):
                t, p = theta[lo : lo + chunk], phi[lo : lo + chunk]
                a, b = _vector_far_field(t, p, x, y, jx, jy, self.f, self.green, dx, dx)
                et[lo : lo + chunk] += a
                ep[lo : lo + chunk] += b
        self._memo[key] = (et, ep)
        return et, ep


class _ComboFF:
    """Weighted sum of single-port far fields (each memoised), for many port-phase draws at no extra cost."""

    def __init__(self, port_ffs, weights):
        self.port_ffs, self.weights = port_ffs, weights

    def __call__(self, theta_deg, phi_deg):
        et = ep = None
        for pid, ff in self.port_ffs.items():
            a, b = ff(theta_deg, phi_deg)
            w = self.weights[pid]
            et = w * a if et is None else et + w * a
            ep = w * b if ep is None else ep + w * b
        return et, ep


def _uv(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.sin(t) * np.cos(p), np.sin(t) * np.sin(p)


def _ang(t1, p1, t0, p0):
    a, b = np.radians(t1), np.radians(p1)
    c, d = np.radians(t0), np.radians(p0)
    return np.degrees(np.arccos(np.clip(np.cos(a) * np.cos(c) + np.sin(a) * np.sin(c) * np.cos(b - d), -1, 1)))


def _scenario_eval(problem, sc, levels, alpha_deg, frequency_hz, grid, amplitudes, ff=None):
    lam = free_space_wavelength(frequency_hz)
    active = set(sc.ports)
    if ff is None:
        ff = _FarField(problem, levels, alpha_deg, frequency_hz, active, amplitudes)
    active_ports = [bp for bp in problem.ports if bp.port_id in active]
    extent = _extent_m(active_ports, problem.period_m)  # reported only
    pol = problem.polarization
    cross_pol = ORTHOGONAL_POLARIZATION[pol]
    theta = np.arange(0.0, 90.0 + 1e-9, grid.theta_step_deg)
    phi = np.arange(0.0, 360.0, grid.phi_step_deg)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    et, ep = ff(tg, pg)
    et, ep = et.reshape(tg.shape), ep.reshape(tg.shape)
    co = np.abs(copol_component(et, ep, pol)) ** 2
    total = np.abs(et) ** 2 + np.abs(ep) ** 2
    w = np.sin(np.radians(tg)) * np.radians(grid.theta_step_deg) * np.radians(grid.phi_step_deg)
    integral = float(np.sum(total * w))
    dir_peak_total = float(total.max())
    d_ref_db = 10.0 * np.log10(4.0 * np.pi * dir_peak_total / integral)

    beams = [(k, problem.beams[k]) for k in sc.beams]
    # fine local peaks (co-pol and total) around every requested beam
    out_beams = []
    peak_co = {}
    for k, (th, ph, _pol) in beams:
        n = int(round(grid.patch_radius_deg / grid.patch_step_deg))
        off = np.arange(-n, n + 1) * grid.patch_step_deg
        dth, dph = np.meshgrid(off, off, indexing="ij")
        s0 = max(np.sin(np.radians(th)), 0.05)
        lt, lp = th + dth, (ph + dph / s0) % 360.0
        pet, pep = ff(lt, lp)
        pco = np.abs(copol_component(pet, pep, pol)) ** 2
        ptot = np.abs(pet) ** 2 + np.abs(pep) ** 2
        pcr = np.abs(copol_component(pet, pep, cross_pol)) ** 2
        i = int(np.argmax(pco))
        peak_co[k] = float(pco[i])
        got_t, got_p = float(lt.ravel()[i]), float(lp.ravel()[i])
        out_beams.append(
            {
                "beam_id": k,
                "port_id": problem.beam_port[k],
                "requested_theta_deg": th,
                "requested_phi_deg": ph,
                "theta_deg": got_t,
                "phi_deg": got_p,
                "pointing_error_deg": float(_ang(got_t, got_p, th, ph)),
                "directivity_dbi": float(d_ref_db + 10.0 * np.log10(max(ptot.ravel()[i], 1e-300) / dir_peak_total)),
                "xpd_db": float(10.0 * np.log10(max(pco.ravel()[i], 1e-300) / max(pcr.ravel()[i], 1e-300))),
                "peak_co": peak_co[k],
            }
        )
    # global SLL per beam: the beam_polarization_metrics partition, peak = fine local co-pol peak
    ids = [k for k, _ in beams]
    # a beam's main-lobe cone follows the extent of the sub-aperture that radiates it
    own = {bp.port_id: _extent_m([bp], problem.period_m) for bp in problem.ports}
    cone = {k: first_null_half_width_deg(own[problem.beam_port[k]], lam, problem.beams[k][0]) for k in ids}
    dist = {k: _ang(tg, pg, problem.beams[k][0], problem.beams[k][1]) for k in ids}
    for ob in out_beams:
        k = ob["beam_id"]
        outside = dist[k] > cone[k]
        for k2 in ids:
            if k2 != k:
                outside &= dist[k] <= dist[k2]
        top = float(co[outside].max()) if outside.any() else 0.0
        peak = max(peak_co[k], float(co[dist[k] <= cone[k]].max()))
        ob["sll_global_db"] = float(10.0 * np.log10(top / peak)) if top > 0 else None
        ob["cone_deg"] = cone[k]
    # x / y cuts through each beam: v = v0 (x-cut) and u = u0 (y-cut), same partition and cone
    u_line = np.arange(-1.0, 1.0 + 1e-9, grid.cut_step_u)
    for ob in out_beams:
        k = ob["beam_id"]
        u0, v0 = _uv(problem.beams[k][0], problem.beams[k][1])
        for name, (uu, vv) in (("x", (u_line, np.full_like(u_line, v0))), ("y", (np.full_like(u_line, u0), u_line))):
            ok = np.hypot(uu, vv) < 1.0 - 1e-6
            uu, vv = uu[ok], vv[ok]
            t_line = np.degrees(np.arcsin(np.hypot(uu, vv)))
            p_line = np.degrees(np.arctan2(vv, uu)) % 360.0
            lt_, lp_ = ff(t_line, p_line)
            pc = np.abs(copol_component(lt_, lp_, pol)) ** 2
            d_own = _ang(t_line, p_line, problem.beams[k][0], problem.beams[k][1])
            outside = d_own > cone[k]
            for k2 in ids:
                if k2 != k:
                    outside &= d_own <= _ang(t_line, p_line, problem.beams[k2][0], problem.beams[k2][1])
            ob[f"sll_{name}cut_db"] = (
                float(
                    10.0
                    * np.log10(pc[outside].max() / max(ob["peak_co"], pc[~outside].max() if (~outside).any() else 0))
                )
                if outside.any() and pc[outside].max() > 0
                else None
            )
    levels_db = [10.0 * np.log10(ob["peak_co"]) for ob in out_beams]
    return {
        "scenario": sc.name,
        "ports": list(sc.ports),
        "peak_directivity_dbi": float(d_ref_db),
        "beams": out_beams,
        "level_spread_db": float(max(levels_db) - min(levels_db)) if len(levels_db) > 1 else 0.0,
        "extent_m": extent,
    }


def _nominal_scenarios(problem):
    return [s for s in problem.scenarios if s.freq_index == 0 and s.constrain]


def evaluate(problem, levels=None, alpha_deg=None, frequency_hz=None, grid=None, scenarios=None, amplitudes=None):
    """Held-out evaluation of the nominal scenarios (each port alone; all ports together with
    `alpha_deg`) -- or of the named ones. `levels`: `{port_id: (N,) int}` (default: the starting
    states). Everything is computed with the vector far field on `grid` (`HeldOutGrid`), never on the
    optimisation directions. At `frequency_hz` != the design frequency the reference wave is moved and,
    when the problem was built with that frequency in `extra_frequencies_hz`, the two cells use their
    RE-SOLVED phases there (otherwise they are held at the design-frequency values)."""
    levels = levels or {p.port_id: p.state for p in problem.opt_ports}
    alpha_deg = alpha_deg or {p.port_id: 0.0 for p in problem.opt_ports}
    grid = grid or HeldOutGrid()
    f = frequency_hz or problem.frequency_hz
    names = scenarios or [s.name for s in _nominal_scenarios(problem)]
    return {
        sc.name: _scenario_eval(problem, sc, levels, alpha_deg, f, grid, amplitudes)
        for sc in problem.scenarios
        if sc.name in names
    }


def evaluate_phase_error_draws(problem, levels, alpha_deg, frequency_hz=None, grid=None, errors=None):
    """The all-ports scenario re-evaluated for each draw of independent port-phase errors (default:
    `problem.draws`) on the held-out grid. The per-port far fields are computed once and recombined.
    Returns a list of `_scenario_eval` dicts, one per draw."""
    levels = levels or {p.port_id: p.state for p in problem.opt_ports}
    alpha_deg = alpha_deg or {p.port_id: 0.0 for p in problem.opt_ports}
    grid = grid or HeldOutGrid()
    f = frequency_hz or problem.frequency_hz
    errors = problem.draws if errors is None else errors
    sc = next(s for s in problem.scenarios if s.name == "all")
    port_ffs = {bp.port_id: _FarField(problem, levels, alpha_deg, f, {bp.port_id}) for bp in problem.ports}
    out = []
    for err in errors:
        w = {pid: np.exp(1j * np.radians(alpha_deg[pid] + err.get(pid, 0.0))) for pid in port_ffs}
        out.append(_scenario_eval(problem, sc, levels, alpha_deg, f, grid, None, ff=_ComboFF(port_ffs, w)))
    return out


def cells_changed(problem, levels):
    return {p.port_id: int(np.sum(np.asarray(levels[p.port_id]) != p.state)) for p in problem.opt_ports}


def patch_count(levels):
    return {pid: int(np.sum(v == 1)) for pid, v in levels.items()}


def apply_levels(result, problem, levels, alpha_deg=None, evaluation=None):
    """A copy of the job `result` with the optimised states: per-port `patches_mm` / `patch_count` /
    `aperture_field_*` / `n_cells` unchanged, the combined `lw_grid_m`, and (when `evaluation` is given)
    each port's `report` from the held-out single-port evaluation. The chosen port phases are stored
    in `result['sidelobe_optimisation']`."""
    import copy

    out = dict(result)
    out["ports"] = []
    state_geom = problem.state_geometry
    for p_res, bp in zip(result["ports"], problem.ports, strict=True):
        p = dict(p_res)
        lev = np.asarray(levels[bp.port_id])
        field_ap = _aperture(problem, bp, lev)
        cells = [
            {"i": int(i), "j": int(j), "geometry": _geom_mm(state_geom[1])}
            for (i, j), lv in zip(bp.cell_ij, lev, strict=True)
            if lv == 1
        ]
        layout = dict(p_res["patches_mm"])
        layout["cells"] = cells
        p["patches_mm"] = layout
        p["patch_count"] = len(cells)
        p["aperture_field_real"] = field_ap.real
        p["aperture_field_imag"] = field_ap.imag
        out["ports"].append(p)
    out["lw_grid_m"] = _combined_grid(result, problem, levels)
    if evaluation is not None:
        for p in out["ports"]:
            ev = evaluation.get(f"port{p['port_id'] + 1}")
            if ev:
                p["report"] = {
                    "peak_directivity_dbi": ev["peak_directivity_dbi"],
                    "aperture_efficiency_upper_bound": (
                        10 ** (ev["peak_directivity_dbi"] / 10)
                        / (
                            4
                            * np.pi
                            * problem.ports[p["port_id"]].region.area_m2
                            / free_space_wavelength(problem.frequency_hz) ** 2
                        )
                    ),
                    "beams": [
                        {
                            "requested_theta_deg": b["requested_theta_deg"],
                            "requested_phi_deg": b["requested_phi_deg"],
                            "theta_deg": b["theta_deg"],
                            "phi_deg": b["phi_deg"],
                            "pointing_error_deg": b["pointing_error_deg"],
                            "level_db": 0.0,
                            "directivity_dbi": b["directivity_dbi"],
                            "polarization": problem.polarization,
                            "cross_pol_level_db": None,
                            "xpd_db": b["xpd_db"],
                            "sll_global_db": b["sll_global_db"],
                        }
                        for b in ev["beams"]
                    ],
                }
    out["sidelobe_optimisation"] = {
        "alpha_deg": {int(k): float(v) for k, v in (alpha_deg or {}).items()},
        "cells_changed": cells_changed(problem, levels),
    }
    del copy
    return out


def _geom_mm(geometry):
    if geometry.get("kind") == "empty":
        return None
    return [round(float(geometry["L_mm"]), 9), round(float(geometry.get("W_mm", geometry["L_mm"])), 9)]


def _combined_grid(result, problem, levels):
    """The combined board grid (`lw_grid_m`) rebuilt from the optimised states, with the same
    owner / nearest-cell rule as `multi_port_synthesis.combined_board_lw_grid_m`."""
    grid = result["lw_grid_m"]
    ny, nx = len(grid), len(grid[0])
    xs = [bp.region.x_min_m for bp in problem.ports]
    ys = [bp.region.y_min_m for bp in problem.ports]
    x_min, y_min = min(xs), min(ys)
    x_max = max(bp.region.x_max_m for bp in problem.ports)
    y_max = max(bp.region.y_max_m for bp in problem.ports)
    px, py = (x_max - x_min) / nx, (y_max - y_min) / ny
    geom1 = problem.state_geometry[1]
    one = None if geom1.get("kind") == "empty" else (geom1["L_mm"] * 1e-3, geom1.get("W_mm", geom1["L_mm"]) * 1e-3)
    new = [[(0.0, 0.0)] * nx for _ in range(ny)]
    for i in range(ny):
        cy = y_min + (i + 0.5) * py
        for j in range(nx):
            cx = x_min + (j + 0.5) * px
            owner = next((bp for bp in problem.ports if bp.region.contains(cx, cy)), None)
            if owner is None:
                continue
            xc, yc = owner.x_m[owner.owned], owner.y_m[owner.owned]
            if xc.size == 0:
                continue
            n = int(np.argmin((xc - cx) ** 2 + (yc - cy) ** 2))
            if int(np.asarray(levels[owner.port_id])[n]) == 1 and one is not None:
                new[i][j] = one
    return new
