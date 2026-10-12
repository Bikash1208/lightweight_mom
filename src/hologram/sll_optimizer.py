"""Phase 23 C: a sidelobe-cancelling optimiser for coded (1-bit) multi-port boards.

The 1-bit state of every cell is optimised DIRECTLY to lower the peak co-polar sidelobe, instead of
only quantising a continuous hologram (Johnson et al., "Sidelobe canceling for optimization of
reconfigurable holographic metamaterial antenna", arXiv 1408.3032 / IEEE TAP 2015). The radiated
field is linear in the per-cell aperture contributions (`pattern.cellwise_pattern`), so flipping one
cell updates the field of every direction in O(directions):

    F(d) <- F(d) + w_p * (s_n(new state) - s_n(old state)) * c[n, d]

**Objective** (per scenario): the peak co-pol sidelobe level relative to the beam's own CURRENT peak,
over every direction outside the main-lobe cone (first-null radius, `first_null_half_width_deg`) of
its nearest requested beam -- the nearest-beam partition of `beam_polarization_metrics`, so one beam's
lobe is never counted as another's sidelobe. The non-smooth maximum is minimised through the
`q`-norm `(sum_d (P_d / peak_b(d))^q)^(1/q)` with `q` raised along a schedule.

**Constraints** (per scenario, checked on every accepted change): each beam's peak drops by at most
`max_gain_loss_db` against the starting design; its peak direction moves by at most
`max_pointing_shift_deg`; the beam levels stay within `max_level_spread_db` of each other (or the
starting spread, if that was wider); the point XPD at the peak is not worse by more than
`max_xpd_loss_db`.

**Scenarios**: a scenario is a set of ports excited together and the beams they radiate -- one port
alone (HFSS Edit Sources one at a time), all ports together with optimisable port phases
`alpha_p` (`alpha` of the first active port is the 0 reference), or both at once (the objective is then
the `q`-norm over the union of all scenarios' sidelobe points, a weighted worst case).

Deterministic and seeded; no learned model (repo rule).
"""

from dataclasses import dataclass, field

import numpy as np

_PATCH_FLOOR = 1e-300


@dataclass(frozen=True)
class OptBeam:
    beam_id: int  # index into the direction set's `tag` / the beam list given to `uv_direction_set`
    port_id: int
    theta_deg: float
    phi_deg: float
    cone_deg: float  # main-lobe exclusion radius (first-null half-width)


@dataclass
class OptPort:
    port_id: int
    c_co: np.ndarray  # (N, D) co-pol contribution of each cell per unit s (cellwise_pattern.cell_contributions)
    c_cross: np.ndarray  # (N, Dp) cross-pol contribution on the patch directions only
    s_states: np.ndarray  # (N, 2) complex: the aperture value of each cell in state 0 / state 1
    state: np.ndarray  # (N,) int, the starting state of every cell
    free: np.ndarray  # (N,) bool, False = fixed (e.g. under the feed)
    amplitude: float = 1.0  # excitation magnitude when several ports are excited together


@dataclass(frozen=True)
class Scenario:
    name: str
    ports: tuple  # port ids excited in this scenario
    beams: tuple  # beam ids radiated (the requested beams of those ports)
    weight: float = 1.0  # weight in the combined objective
    cone_deg: "tuple | None" = None  # main-lobe exclusion radius per beam in this scenario (default: OptBeam.cone_deg)


@dataclass
class OptimiserConfig:
    max_gain_loss_db: float = 0.5
    max_pointing_shift_deg: float = 0.3
    max_level_spread_db: float = 0.5
    max_xpd_loss_db: float = 1.0
    max_flips: int = 600  # iteration budget: accepted single-cell changes
    q_schedule: tuple = (6.0, 12.0, 24.0, 48.0)
    accepts_per_round: int = 12
    optimise_alpha: bool = False  # shared port phases of multi-port scenarios
    alpha_step_deg: float = 15.0
    anneal_steps: int = 0  # simulated-annealing proposals after the greedy phase (0 = off)
    anneal_t0: float = 0.002  # initial temperature, in units of the relative q-norm change
    seed: int = 0
    chunk: int = 192  # candidates scored at once


@dataclass
class ScenarioReport:
    name: str
    sll_db: float  # peak sidelobe / own beam peak on the optimisation directions
    beam_peak_db: dict  # beam_id -> peak level change vs the starting design [dB]
    beam_shift_deg: dict  # beam_id -> peak-direction change vs the starting design
    level_spread_db: float
    xpd_db: dict  # beam_id -> point XPD at the peak


@dataclass
class OptimiserResult:
    state: dict  # port_id -> (N,) int, the optimised states
    alpha_deg: dict  # port_id -> phase [deg] (0 for ports not in a multi-port scenario)
    flips: dict  # port_id -> list of cell indices whose state changed
    before: dict  # scenario name -> ScenarioReport
    after: dict
    objective_before: float
    objective_after: float
    history: list = field(default_factory=list)  # (accepted flips, q, objective) after each round
    n_evaluations: int = 0


# --------------------------------------------------------------------------------------- geometry


def _angular_distance_deg(theta_deg, phi_deg, theta0_deg, phi0_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    t0, p0 = np.radians(theta0_deg), np.radians(phi0_deg)
    c = np.cos(t) * np.cos(t0) + np.sin(t) * np.sin(t0) * np.cos(p - p0)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


class _ScenarioContext:
    """Per-scenario index sets on the direction set: sidelobe points and each beam's fine patch."""

    def __init__(self, scenario, beams_by_id, directions):
        self.scenario = scenario
        self.beam_ids = list(scenario.beams)
        beams = [beams_by_id[b] for b in self.beam_ids]
        theta, phi = directions.theta_deg, directions.phi_deg
        dist = np.stack([_angular_distance_deg(theta, phi, b.theta_deg, b.phi_deg) for b in beams])  # (B, D)
        nearest = np.argmin(dist, axis=0)
        d_near = dist[nearest, np.arange(theta.size)]
        cones = scenario.cone_deg if scenario.cone_deg is not None else [b.cone_deg for b in beams]
        cone = np.array(cones, dtype=float)[nearest]
        outside = d_near > cone
        self.s_idx = np.nonzero(outside)[0]
        self.s_beam = nearest[outside]  # local beam index of each sidelobe point
        self.patch_idx = [np.nonzero(directions.tag == bid)[0] for bid in self.beam_ids]
        if any(p.size == 0 for p in self.patch_idx):
            raise ValueError(f"scenario {scenario.name!r}: a beam has no fine patch points (tag the direction set)")
        # patch lattice neighbours (for the parabolic peak refinement), in the patch's own indexing
        self.patch_u = [directions.u[p] for p in self.patch_idx]
        self.patch_v = [directions.v[p] for p in self.patch_idx]
        self.neigh = [self._neighbours(u, v) for u, v in zip(self.patch_u, self.patch_v, strict=True)]

    @staticmethod
    def _neighbours(u, v):
        """For each patch point the index of its +u, -u, +v, -v lattice neighbour (-1 if absent)."""
        step = np.min(np.diff(np.unique(np.round(u, 9)))) if np.unique(np.round(u, 9)).size > 1 else 1.0
        key = {(int(round(a / step)), int(round(b / step))): k for k, (a, b) in enumerate(zip(u, v, strict=True))}
        out = np.full((u.size, 4), -1, dtype=int)
        for k, (a, b) in enumerate(zip(u, v, strict=True)):
            iu, iv = int(round(a / step)), int(round(b / step))
            for j, (du, dv) in enumerate(((1, 0), (-1, 0), (0, 1), (0, -1))):
                out[k, j] = key.get((iu + du, iv + dv), -1)
        return out, step


# ------------------------------------------------------------------------------------ evaluation


def _refined_peak(power_patch, neigh_pack, u, v):
    """`(peak_power, u, v)` of one beam's fine patch: the lattice maximum refined by a parabola through
    its neighbours (in ln power), separately in u and v."""
    neigh, step = neigh_pack
    k = int(np.argmax(power_patch))
    p00 = np.log(max(power_patch[k], _PATCH_FLOOR))
    uu, vv = float(u[k]), float(v[k])
    ln_peak = p00
    for (jm, jp), axis in (((1, 0), 0), ((3, 2), 1)):  # (minus, plus) neighbour columns of `neigh`
        a, b = neigh[k, jm], neigh[k, jp]
        if a < 0 or b < 0:
            continue
        pm = np.log(max(power_patch[a], _PATCH_FLOOR))
        pp = np.log(max(power_patch[b], _PATCH_FLOOR))
        denom = pm - 2.0 * p00 + pp
        if denom < -1e-12:
            off = float(np.clip(0.5 * (pm - pp) / denom, -1.0, 1.0))
            ln_peak += -((pp - pm) ** 2) / (8.0 * denom)
            if axis == 0:
                uu += off * step
            else:
                vv += off * step
    return float(np.exp(ln_peak)), uu, vv


def _uv_to_theta_phi_deg(u, v):
    s = float(np.hypot(u, v))
    return float(np.degrees(np.arcsin(min(s, 1.0)))), float(np.degrees(np.arctan2(v, u)) % 360.0)


class _Evaluator:
    """Holds the per-port fields and evaluates scenarios from them."""

    def __init__(self, ports, scenarios, contexts, directions, alpha_deg):
        self.ports = {p.port_id: p for p in ports}
        self.scenarios = scenarios
        self.ctx = contexts
        self.directions = directions
        self.alpha = dict(alpha_deg)
        self.patch_all = np.nonzero(directions.tag >= 0)[0]
        self.patch_pos = {int(d): k for k, d in enumerate(self.patch_all)}
        self.F = {}
        self.Fx = {}

    # -- fields
    def set_states(self, state_by_port):
        for pid, p in self.ports.items():
            s = p.s_states[np.arange(p.state.size), state_by_port[pid]]
            self.F[pid] = s @ p.c_co
            self.Fx[pid] = s @ p.c_cross

    def weight(self, scenario, pid):
        if len(scenario.ports) == 1:
            return 1.0
        return self.ports[pid].amplitude * np.exp(1j * np.radians(self.alpha[pid]))

    def scenario_field(self, scenario, F=None):
        F = self.F if F is None else F
        return sum(self.weight(scenario, pid) * F[pid] for pid in scenario.ports)

    def scenario_cross_patch(self, scenario, Fx=None):
        Fx = self.Fx if Fx is None else Fx
        return sum(self.weight(scenario, pid) * Fx[pid] for pid in scenario.ports)

    # -- metrics
    def peaks(self, scenario, field, cross):
        """Per beam of the scenario: (refined peak power, theta, phi, point XPD in dB)."""
        ctx = self.ctx[scenario.name]
        power = np.abs(field) ** 2
        out = []
        for b, pidx in enumerate(ctx.patch_idx):
            pp = power[pidx]
            peak, uu, vv = _refined_peak(pp, ctx.neigh[b], ctx.patch_u[b], ctx.patch_v[b])
            th, ph = _uv_to_theta_phi_deg(uu, vv)
            k = int(np.argmax(pp))
            pc = abs(cross[self.patch_pos[int(pidx[k])]]) ** 2
            xpd = 10.0 * np.log10(max(pp[k], _PATCH_FLOOR) / max(pc, _PATCH_FLOOR))
            out.append((peak, th, ph, xpd))
        return out

    def q_term(self, scenario, field, q):
        """`sum_d (P_d / peak_b(d))^q` over the scenario's sidelobe points."""
        ctx = self.ctx[scenario.name]
        power = np.abs(field) ** 2
        peaks = np.array([power[pidx].max() for pidx in ctx.patch_idx])
        ratio = power[ctx.s_idx] / np.maximum(peaks[ctx.s_beam], _PATCH_FLOOR)
        return float(np.sum(ratio**q))

    def sll_db(self, scenario, field):
        ctx = self.ctx[scenario.name]
        power = np.abs(field) ** 2
        peaks = np.array([power[pidx].max() for pidx in ctx.patch_idx])
        ratio = power[ctx.s_idx] / np.maximum(peaks[ctx.s_beam], _PATCH_FLOOR)
        return float(10.0 * np.log10(max(ratio.max(), _PATCH_FLOOR)))

    def total_term(self, q, F=None):
        return sum(sc.weight * self.q_term(sc, self.scenario_field(sc, F), q) for sc in self.scenarios)

    def report(self, baseline=None):
        out = {}
        for sc in self.scenarios:
            field = self.scenario_field(sc)
            cross = self.scenario_cross_patch(sc)
            pk = self.peaks(sc, field, cross)
            ctx = self.ctx[sc.name]
            levels = [10.0 * np.log10(max(p[0], _PATCH_FLOOR)) for p in pk]
            if baseline is None:
                dpk = {bid: 0.0 for bid in ctx.beam_ids}
                shift = {bid: 0.0 for bid in ctx.beam_ids}
            else:
                b0 = baseline[sc.name]
                dpk = {bid: levels[i] - b0["levels"][i] for i, bid in enumerate(ctx.beam_ids)}
                shift = {
                    bid: float(_angular_distance_deg(pk[i][1], pk[i][2], b0["peaks"][i][1], b0["peaks"][i][2]))
                    for i, bid in enumerate(ctx.beam_ids)
                }
            out[sc.name] = {
                "report": ScenarioReport(
                    sc.name,
                    self.sll_db(sc, field),
                    dpk,
                    shift,
                    float(max(levels) - min(levels)) if len(levels) > 1 else 0.0,
                    {bid: pk[i][3] for i, bid in enumerate(ctx.beam_ids)},
                ),
                "levels": levels,
                "peaks": pk,
            }
        return out


# --------------------------------------------------------------------------------------- optimiser


def optimise_sidelobes(ports, beams, scenarios, directions, config=None, alpha_deg=None):
    """Optimises the states of every free cell of `ports` (and, with `config.optimise_alpha`, the port
    phases of multi-port scenarios) to minimise the combined `q`-norm of the scenarios' sidelobes under
    the module's constraints. Returns an `OptimiserResult`. Reproducible: same inputs and `seed` give
    the same states.

    `beams`: list of `OptBeam`; `directions`: a `DirectionSet` whose `tag` marks beam `k`'s fine patch
    (`uv_direction_set(beams_deg=...)` in the same order as `beams`)."""
    cfg = config or OptimiserConfig()
    rng = np.random.default_rng(cfg.seed)
    beams_by_id = {b.beam_id: b for b in beams}
    contexts = {sc.name: _ScenarioContext(sc, beams_by_id, directions) for sc in scenarios}
    port_ids = [p.port_id for p in ports]
    alpha = {pid: 0.0 for pid in port_ids}
    if alpha_deg:
        alpha.update({pid: float(a) for pid, a in alpha_deg.items()})
    ev = _Evaluator(ports, scenarios, contexts, directions, alpha)
    state = {p.port_id: np.array(p.state, dtype=int, copy=True) for p in ports}
    ev.set_states(state)

    base = ev.report(baseline=None)
    baseline = {name: {"levels": v["levels"], "peaks": v["peaks"]} for name, v in base.items()}
    spread0 = {name: (max(v["levels"]) - min(v["levels"])) if len(v["levels"]) > 1 else 0.0 for name, v in base.items()}
    xpd0 = {name: [p[3] for p in v["peaks"]] for name, v in base.items()}
    obj_before = ev.total_term(cfg.q_schedule[-1]) ** (1.0 / cfg.q_schedule[-1])
    n_eval = [0]

    def feasible(F=None, Fx=None, alpha_override=None):
        """All constraints of every scenario for the given per-port fields."""
        saved = ev.alpha
        if alpha_override is not None:
            ev.alpha = alpha_override
        try:
            for sc in scenarios:
                field = ev.scenario_field(sc, F)
                cross = ev.scenario_cross_patch(sc, Fx)
                pk = ev.peaks(sc, field, cross)
                b0 = baseline[sc.name]
                levels = [10.0 * np.log10(max(p[0], _PATCH_FLOOR)) for p in pk]
                for i in range(len(pk)):
                    if levels[i] - b0["levels"][i] < -cfg.max_gain_loss_db:
                        return False
                    shift = _angular_distance_deg(pk[i][1], pk[i][2], b0["peaks"][i][1], b0["peaks"][i][2])
                    if shift > cfg.max_pointing_shift_deg:
                        return False
                    if pk[i][3] < xpd0[sc.name][i] - cfg.max_xpd_loss_db:
                        return False
                if len(levels) > 1 and max(levels) - min(levels) > max(cfg.max_level_spread_db, spread0[sc.name]):
                    return False
            return True
        finally:
            ev.alpha = saved

    def candidate_terms(pid, idx, q):
        """Total q-term if each cell in `idx` of port `pid` were flipped (others unchanged)."""
        p = ev.ports[pid]
        delta = p.s_states[idx, 1 - state[pid][idx]] - p.s_states[idx, state[pid][idx]]
        involved = [sc for sc in scenarios if pid in sc.ports]
        others = sum(sc.weight * ev.q_term(sc, ev.scenario_field(sc), q) for sc in scenarios if pid not in sc.ports)
        totals = np.full(idx.size, others, dtype=float)
        for sc in involved:
            ctx = contexts[sc.name]
            base_field = ev.scenario_field(sc)
            w = ev.weight(sc, pid)
            cols = ctx.s_idx
            patch_cat = np.concatenate(ctx.patch_idx)
            patch_split = np.cumsum([pi.size for pi in ctx.patch_idx])[:-1]
            for lo in range(0, idx.size, cfg.chunk):
                sl = idx[lo : lo + cfg.chunk]
                rows = p.c_co[sl]  # (K, D)
                dlt = (w * delta[lo : lo + cfg.chunk])[:, None]
                p_patch = np.abs(base_field[patch_cat][None, :] + dlt * rows[:, patch_cat]) ** 2
                peaks = np.stack([blk.max(axis=1) for blk in np.split(p_patch, patch_split, axis=1)], axis=1)  # (K, B)
                p_s = np.abs(base_field[cols][None, :] + dlt * rows[:, cols]) ** 2
                ratio = p_s / np.maximum(peaks[:, ctx.s_beam], _PATCH_FLOOR)
                totals[lo : lo + sl.size] += sc.weight * np.sum(ratio**q, axis=1)
            n_eval[0] += idx.size
        return totals

    def apply_flip(pid, n):
        p = ev.ports[pid]
        old, new = state[pid][n], 1 - state[pid][n]
        d = p.s_states[n, new] - p.s_states[n, old]
        ev.F[pid] = ev.F[pid] + d * p.c_co[n]
        ev.Fx[pid] = ev.Fx[pid] + d * p.c_cross[n]
        state[pid][n] = new

    def undo_flip(pid, n):
        apply_flip(pid, n)  # a flip is its own inverse

    history = []
    accepted = 0
    changed = {pid: set() for pid in port_ids}

    def optimise_alpha_step(q):
        """Coordinate search of the shared port phases of multi-port scenarios."""
        multi = [sc for sc in scenarios if len(sc.ports) > 1]
        if not multi:
            return False
        movers = sorted({pid for sc in multi for pid in sc.ports})[1:]  # the first port is the 0 reference
        improved = False
        for pid in movers:
            cur = ev.alpha[pid]
            best_t = ev.total_term(q)
            best_a = cur
            for a in np.arange(0.0, 360.0, cfg.alpha_step_deg):
                if abs(((a - cur + 180.0) % 360.0) - 180.0) < 1e-9:
                    continue
                trial = dict(ev.alpha)
                trial[pid] = float(a)
                saved = ev.alpha
                ev.alpha = trial
                t = ev.total_term(q)
                ok = t < best_t * (1.0 - 1e-9) and feasible()
                ev.alpha = saved
                if ok:
                    best_t, best_a = t, float(a)
            if best_a != cur:
                ev.alpha[pid] = best_a
                improved = True
        return improved

    for q in cfg.q_schedule:
        t_cur = ev.total_term(q)
        for _round in range(10_000):
            if accepted >= cfg.max_flips:
                break
            ranked = []
            for pid in port_ids:
                p = ev.ports[pid]
                idx = np.nonzero(p.free)[0]
                if idx.size == 0:
                    continue
                terms = candidate_terms(pid, idx, q)
                good = np.nonzero(terms < t_cur * (1.0 - 1e-12))[0]
                ranked += [(float(terms[g]), pid, int(idx[g])) for g in good]
            ranked.sort(key=lambda r: (r[0], r[1], r[2]))
            got = 0
            for tried, (_t, pid, n) in enumerate(ranked):
                if got >= cfg.accepts_per_round or accepted >= cfg.max_flips or tried > 6 * cfg.accepts_per_round:
                    break
                apply_flip(pid, n)
                t_new = ev.total_term(q)
                if t_new < t_cur * (1.0 - 1e-12) and feasible():
                    t_cur = t_new
                    accepted += 1
                    got += 1
                    changed[pid].symmetric_difference_update({n})
                else:
                    undo_flip(pid, n)
            if cfg.optimise_alpha and got == 0 and optimise_alpha_step(q):
                t_cur = ev.total_term(q)
                got = 1
            history.append((accepted, float(q), float(t_cur ** (1.0 / q))))
            if got == 0:
                break
        # refresh the exact fields once per q (drops accumulated rounding of the incremental updates)
        ev.set_states(state)

    if cfg.optimise_alpha and any(len(sc.ports) > 1 for sc in scenarios):
        optimise_alpha_step(cfg.q_schedule[-1])

    if cfg.anneal_steps > 0:
        q = cfg.q_schedule[-1]
        t_cur = ev.total_term(q)
        free = [(pid, n) for pid in port_ids for n in np.nonzero(ev.ports[pid].free)[0]]
        best_state = {pid: state[pid].copy() for pid in port_ids}
        best_t = t_cur
        for k in range(cfg.anneal_steps):
            temp = cfg.anneal_t0 * (1.0 - k / cfg.anneal_steps) + 1e-9
            pid, n = free[int(rng.integers(len(free)))]
            apply_flip(pid, int(n))
            t_new = ev.total_term(q)
            n_new, n_cur = t_new ** (1.0 / q), t_cur ** (1.0 / q)
            rel = (n_new - n_cur) / max(n_cur, _PATCH_FLOOR)  # relative change of the q-norm
            if (t_new < t_cur or rng.random() < np.exp(-rel / temp)) and feasible():
                t_cur = t_new
                changed[pid].symmetric_difference_update({int(n)})
                if t_cur < best_t:
                    best_t = t_cur
                    best_state = {p_: state[p_].copy() for p_ in port_ids}
            else:
                undo_flip(pid, int(n))
        if best_t < ev.total_term(q):
            for pid in port_ids:
                state[pid][:] = best_state[pid]
            ev.set_states(state)
            changed = {pid: set(np.nonzero(state[pid] != ev.ports[pid].state)[0].tolist()) for pid in port_ids}

    ev.set_states(state)
    after = ev.report(baseline=baseline)
    q_last = cfg.q_schedule[-1]
    obj_after = ev.total_term(q_last) ** (1.0 / q_last)
    return OptimiserResult(
        state={pid: state[pid].copy() for pid in port_ids},
        alpha_deg=dict(ev.alpha),
        flips={pid: sorted(int(i) for i in np.nonzero(state[pid] != ev.ports[pid].state)[0]) for pid in port_ids},
        before={k: v["report"] for k, v in base.items()},
        after={k: v["report"] for k, v in after.items()},
        objective_before=float(obj_before),
        objective_after=float(obj_after),
        history=history,
        n_evaluations=n_eval[0],
    )
