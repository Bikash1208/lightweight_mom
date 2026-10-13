"""Phase 25 A: the far field of a FIXED cell map under an ARBITRARY monopole excitation vector.

Every earlier pattern path (`sll_board._FarField`, `compute_vector_pattern`) evaluates the design-time
`aperture_field`, i.e. the field of exactly the excitation the board was synthesised for. A Butler-fed
sector scanner needs the pattern of the SAME cells under a different complex excitation
`a = (a_1 ... a_M)` of the M monopoles (Butler input k, a blend of two inputs, one monopole alone).

Model (the one the synthesis already assumes, made linear in `a`):

    s_m,n  = E_m(x_n, y_n) * Gamma_n                scalar aperture value of cell n for feed m
    J_m,n  = rho_hat_m(x_n, y_n) * s_m,n            radial current about feed m
    E(d)   = sum_m a_m * G(d) . FT[J_m](d) cos(theta)    existing vector far field, `_vector_far_field`

* `E_m` is the surface wave of monopole `m` alone (`phase8b.reference_wave.monopole_reference_field`
  about the feed, the same function `hologram.multi_feed` wraps), `Gamma_n` the per-cell reflection
  phasor of the state the cell carries (`exp(j psi_state)`, the exact achieved phases the synthesis used).
* `amplitude_model='phase_only'` (default) uses `E_m = exp(-j k_swg rho_m)`, the convention of the
  single-port synthesis (`ports.port_reference_field` normalises the monopole to unit modulus, so every
  existing `aperture_field` is `phase-only reference x exp(j psi)`). `'physical'` keeps the cylindrical
  spreading, `E_m = rho_m^(-1/2) exp(-j k_swg rho_m)` (normalised to 1 at a distance `reference_rho_m`,
  `rho` floored at half a period at the feed cell). They differ only in how the feeds' relative levels fall
  off across the board; both are reported where it matters.
* The current of feed m is radial about feed m, so a multi-feed excitation is the vector sum of the
  feeds' radial currents (for one feed this IS the single-port model, bit for bit).
* A cell may be illuminated by only some feeds (`footprint` zeros in `e_inc`): the segmented layouts
  (one region per port, each lit only by its own monopole) are the special case where feed m's footprint
  is region m.

Approximation (documented, not hidden): first order in the surface-wave/aperture interaction -- the cell
radiates `Gamma_state x E_inc` (local periodicity: an isolated periodic cell's reflection phasor, the
state's own `achieved_phase`), the board does not perturb the surface wave, there is NO feed-feed coupling
(S_mn of the monopoles is not modelled, a sum of independent feeds), no surface-wave loss across the board,
and only the radiated field is modelled (no scattering back into the feeds). What a plane-wave full-array
MoM can and cannot check is spelled out in `scripts/research/phase25_forward_validation.py`.

Linearity: `E(a) = sum_m a_m E_m` exactly (`FeedFarField.pattern`), so the far field of any excitation
-- every Butler input, every blend, a sweep of the relative phase psi -- costs M far-field evaluations in
total, computed once per direction set.
"""

from dataclasses import dataclass, field

import numpy as np

from ..green_function import LayeredGreenFunction
from ..pattern.vector_aperture_pattern import (
    _vector_far_field,
    aperture_current_from_scalar,
)
from ..phase8b.reference_wave import monopole_reference_field

AMPLITUDE_MODELS = ("phase_only", "physical")


def monopole_incident_field(
    x_m, y_m, feed_xy_m, k_swg, amplitude_model="phase_only", period_m=None, reference_rho_m=1.0
):
    """Complex incident surface wave of ONE monopole at `feed_xy_m` on the points `(x_m, y_m)` (module
    docstring). `period_m` is only used by the `'physical'` model, to floor `rho` at half a cell."""
    if amplitude_model not in AMPLITUDE_MODELS:
        raise ValueError(f"amplitude_model must be one of {AMPLITUDE_MODELS}, got {amplitude_model!r}")
    x, y = np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
    dx, dy = x - feed_xy_m[0], y - feed_xy_m[1]
    if amplitude_model == "phase_only":
        return np.exp(-1j * k_swg * np.hypot(dx, dy))
    rho = np.hypot(dx, dy)
    if period_m is not None:
        # floor rho at half a period for the amplitude only (the singular cell under the feed)
        floor = 0.5 * period_m
        amp = np.sqrt(reference_rho_m / np.maximum(rho, floor))
        return amp * np.exp(-1j * k_swg * rho)
    return monopole_reference_field(dx, dy, k_swg, np.sqrt(reference_rho_m))


@dataclass(frozen=True)
class ForwardBoard:
    """A fixed cell map as flat arrays over its N cells.

    `gamma`: per-cell reflection phasor of the state it carries (complex, N). `e_inc`: `(M, N)`
    complex, the incident field of each feed alone at each cell (0 where the feed does not light the
    cell). `feeds_xy_m`: `(M, 2)`. `period_m`: the square cell pitch (the aperture-sum weight)."""

    x_m: np.ndarray
    y_m: np.ndarray
    gamma: np.ndarray
    e_inc: np.ndarray
    feeds_xy_m: np.ndarray
    period_m: float
    frequency_hz: float
    substrate: object
    polarization: str = "rhcp"
    label: str = ""
    amplitude_model: str = "phase_only"
    k_swg_rad_per_m: float = float("nan")
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        n = np.asarray(self.x_m).size
        if np.asarray(self.y_m).size != n or np.asarray(self.gamma).size != n:
            raise ValueError("x_m, y_m and gamma must have the same number of cells")
        if np.asarray(self.e_inc).shape != (np.asarray(self.feeds_xy_m).shape[0], n):
            raise ValueError("e_inc must be (n_feeds, n_cells)")

    @property
    def n_cells(self):
        return int(np.asarray(self.x_m).size)

    @property
    def n_feeds(self):
        return int(np.asarray(self.feeds_xy_m).shape[0])

    def aperture_values(self, excitation):
        """`s_n = sum_m a_m E_m,n Gamma_n` -- the scalar aperture value of every cell (the quantity the
        existing `aperture_field` stores when `excitation` is the designed one)."""
        a = np.asarray(excitation, dtype=complex)
        if a.shape != (self.n_feeds,):
            raise ValueError(f"excitation must have {self.n_feeds} entries, got shape {a.shape}")
        return (a @ self.e_inc) * np.asarray(self.gamma)

    def with_states(self, gamma):
        """The same board with other per-cell states (a re-synthesised or perturbed map)."""
        return ForwardBoard(
            self.x_m,
            self.y_m,
            np.asarray(gamma, dtype=complex),
            self.e_inc,
            self.feeds_xy_m,
            self.period_m,
            self.frequency_hz,
            self.substrate,
            self.polarization,
            self.label,
            self.amplitude_model,
            self.k_swg_rad_per_m,
            self.meta,
        )


def board_from_cells(
    x_m,
    y_m,
    gamma,
    feeds_xy_m,
    period_m,
    frequency_hz,
    substrate,
    k_swg_rad_per_m,
    polarization="rhcp",
    amplitude_model="phase_only",
    footprints=None,
    label="",
    meta=None,
):
    """A `ForwardBoard` with the analytic monopole field of every feed on the cells. `footprints`:
    optional `(M, N)` bool, which cells feed m lights (default: every cell, every feed)."""
    x = np.asarray(x_m, dtype=float).ravel()
    y = np.asarray(y_m, dtype=float).ravel()
    feeds = np.asarray(feeds_xy_m, dtype=float).reshape(-1, 2)
    e = np.stack([monopole_incident_field(x, y, tuple(f), k_swg_rad_per_m, amplitude_model, period_m) for f in feeds])
    if footprints is not None:
        fp = np.asarray(footprints, dtype=bool)
        if fp.shape != e.shape:
            raise ValueError(f"footprints must be {e.shape}, got {fp.shape}")
        e = np.where(fp, e, 0.0)
    return ForwardBoard(
        x,
        y,
        np.asarray(gamma, dtype=complex).ravel(),
        e,
        feeds,
        float(period_m),
        float(frequency_hz),
        substrate,
        polarization,
        label,
        amplitude_model,
        float(k_swg_rad_per_m),
        dict(meta or {}),
    )


def board_from_ports_result(
    result, period_mm, frequency_hz, eps_r, tan_delta, thickness_mm, polarization="rhcp", incident="analytic"
):
    """A `ForwardBoard` from a coded 1-bit `/synthesize_ports` job result: every port's owned cells,
    each lit only by its own monopole (the segmented layout), `Gamma = exp(j psi_state)` with the exact
    achieved phases `result['two_state']['achieved_phase_deg']`. The incident field is rebuilt
    ANALYTICALLY at the cell centres (phase-only TM0, bare-slab beta -- the convention
    `aperture_field` was built with) when `incident='analytic'`, so comparing `aperture_values(e_p)` with
    the stored field measures how closely the forward model's own incident field reproduces the synthesis
    (the synthesis interpolates a finer native grid, which differs from the analytic phasor mostly in the
    few cells next to the feed). `incident='stored'` recovers each port's own `E_ref = aperture_field /
    Gamma` instead -- the exact incident field of the design, used to prove the pattern path reproduces
    the synthesis prediction to rounding."""
    from ..dispersion.analytic_tm0 import tm0_beta
    from ..materials import Substrate
    from .multi_port_synthesis import region_cell_mask
    from .sll_board import regions_from_result

    two_state = result.get("two_state")
    if not two_state or "achieved_phase_deg" not in two_state:
        raise ValueError("needs a coded 1-bit multi-port design result (two_state.achieved_phase_deg)")
    psi = [float(np.radians(v)) for v in two_state["achieved_phase_deg"]]
    period_m = period_mm * 1e-3
    substrate = Substrate("user", eps_r, tan_delta, thickness_mm * 1e-3)
    k_swg = tm0_beta(frequency_hz, eps_r, thickness_mm * 1e-3)
    regions = regions_from_result(result)
    xs, ys, gs, owner = [], [], [], []
    feeds = []
    stored = []
    for p, region in zip(result["ports"], regions, strict=True):
        x_m, y_m = np.asarray(p["x_m"]), np.asarray(p["y_m"])
        earlier = [r for r in regions if r.port_id < region.port_id]
        owned = region_cell_mask(region, x_m, y_m, earlier) if region.is_polygon else np.ones(x_m.shape, dtype=bool)
        levels = np.zeros(x_m.shape, dtype=int)
        for c in p["patches_mm"]["cells"]:
            levels[c["i"], c["j"]] = 1
        gamma = np.where(levels == 1, np.exp(1j * psi[1]), np.exp(1j * psi[0]))
        xs.append(x_m[owned])
        ys.append(y_m[owned])
        gs.append(gamma[owned])
        owner.append(np.full(int(owned.sum()), region.port_id))
        feeds.append((p["feed_x_mm"] * 1e-3, p["feed_y_mm"] * 1e-3))
        stored.append((np.asarray(p["aperture_field_real"]) + 1j * np.asarray(p["aperture_field_imag"]))[owned])
    x, y, g, own = np.concatenate(xs), np.concatenate(ys), np.concatenate(gs), np.concatenate(owner)
    footprints = np.stack([own == k for k in range(len(feeds))])
    if incident not in ("analytic", "stored"):
        raise ValueError(f"incident must be 'analytic' or 'stored', got {incident!r}")
    board = board_from_cells(
        x,
        y,
        g,
        feeds,
        period_m,
        frequency_hz,
        substrate,
        k_swg,
        polarization,
        "phase_only",
        footprints,
        label="ports_result",
        meta={"stored_aperture_field": np.concatenate(stored), "owner_port": own},
    )
    if incident == "stored":
        e = np.zeros_like(board.e_inc)
        s_all = np.concatenate(stored)
        for k in range(len(feeds)):
            e[k] = np.where(own == k, s_all / g, 0.0)
        board = ForwardBoard(
            board.x_m,
            board.y_m,
            board.gamma,
            e,
            board.feeds_xy_m,
            board.period_m,
            board.frequency_hz,
            board.substrate,
            board.polarization,
            "ports_result_stored_incident",
            board.amplitude_model,
            board.k_swg_rad_per_m,
            board.meta,
        )
    return board


class FeedFarField:
    """`(E_theta, E_phi)` of every feed ALONE (unit excitation) at arbitrary points, memoised per point
    set -- every excitation's pattern is then the weighted sum `sum_m a_m E_m` (module docstring,
    linearity), so a sweep over excitations (every Butler input, a psi sweep) costs M transforms in
    total per point set."""

    def __init__(self, board: ForwardBoard, chunk=2000):
        self.board = board
        self.chunk = chunk
        self._green = LayeredGreenFunction(board.substrate, board.frequency_hz, backing="pec")
        self._memo = {}
        self._current = []
        for m in range(board.n_feeds):
            lit = np.abs(board.e_inc[m]) > 0.0
            s = board.e_inc[m] * board.gamma
            jx, jy = aperture_current_from_scalar(board.x_m, board.y_m, s, feed_xy_m=tuple(board.feeds_xy_m[m]))
            self._current.append((board.x_m[lit], board.y_m[lit], jx[lit], jy[lit]))

    def per_feed(self, theta_deg, phi_deg):
        """`(E_theta, E_phi)`, each `(M, P)`: feed m alone, unit excitation, at the P points."""
        theta = np.asarray(theta_deg, dtype=float).ravel()
        phi = np.asarray(phi_deg, dtype=float).ravel()
        key = (theta.size, float(theta.sum()), float(phi.sum()), float(theta[0]), float(theta[-1]), float(phi[-1]))
        hit = self._memo.get(key)
        if hit is not None:
            return hit
        b = self.board
        et = np.zeros((b.n_feeds, theta.size), dtype=complex)
        ep = np.zeros((b.n_feeds, theta.size), dtype=complex)
        for m, (x, y, jx, jy) in enumerate(self._current):
            for lo in range(0, theta.size, self.chunk):
                sl = slice(lo, lo + self.chunk)
                a, c = _vector_far_field(
                    theta[sl], phi[sl], x, y, jx, jy, b.frequency_hz, self._green, b.period_m, b.period_m
                )
                et[m, sl], ep[m, sl] = a, c
        self._memo[key] = (et, ep)
        return et, ep

    def fields(self, excitation, theta_deg, phi_deg):
        """`(E_theta, E_phi)` for the complex excitation vector `a`: `sum_m a_m E_m`."""
        a = np.asarray(excitation, dtype=complex)
        if a.shape != (self.board.n_feeds,):
            raise ValueError(f"excitation must have {self.board.n_feeds} entries, got shape {a.shape}")
        et, ep = self.per_feed(theta_deg, phi_deg)
        return a @ et, a @ ep

    def as_callable(self, excitation):
        """`(theta, phi) -> (E_theta, E_phi)`, the `ff=` argument of `sll_board._scenario_eval`."""
        a = np.asarray(excitation, dtype=complex)
        return lambda theta_deg, phi_deg: self.fields(a, theta_deg, phi_deg)


def far_field(board: ForwardBoard, excitation, theta_deg, phi_deg, chunk=2000):
    """One-off convenience: `(E_theta, E_phi)` of `board` under `excitation`, summing the excitation's
    CURRENT first (one transform instead of one per feed). Equal to `FeedFarField.fields` by linearity
    (tested)."""
    green = LayeredGreenFunction(board.substrate, board.frequency_hz, backing="pec")
    a = np.asarray(excitation, dtype=complex)
    jx = np.zeros(board.n_cells, dtype=complex)
    jy = np.zeros(board.n_cells, dtype=complex)
    for m in range(board.n_feeds):
        if a[m] == 0:
            continue
        s = a[m] * board.e_inc[m] * board.gamma
        ax, ay = aperture_current_from_scalar(board.x_m, board.y_m, s, feed_xy_m=tuple(board.feeds_xy_m[m]))
        jx += ax
        jy += ay
    theta = np.asarray(theta_deg, dtype=float).ravel()
    phi = np.asarray(phi_deg, dtype=float).ravel()
    et = np.zeros(theta.size, dtype=complex)
    ep = np.zeros(theta.size, dtype=complex)
    lit = (jx != 0) | (jy != 0)
    for lo in range(0, theta.size, chunk):
        sl = slice(lo, lo + chunk)
        a_, b_ = _vector_far_field(
            theta[sl],
            phi[sl],
            board.x_m[lit],
            board.y_m[lit],
            jx[lit],
            jy[lit],
            board.frequency_hz,
            green,
            board.period_m,
            board.period_m,
        )
        et[sl], ep[sl] = a_, b_
    return et, ep


def evaluate_excitation(
    board: ForwardBoard, ffcache: FeedFarField, excitation, beams, cell_port=None, beam_port=None, grid=None, name="exc"
):
    """Held-out metrics of one excitation with the repo's own evaluation code (`sll_board._scenario_eval`:
    directivity, beam peak and pointing on a fine local grid, global and x/y-cut co-pol SLL with the
    first-null cone and nearest-beam partition, XPD, level spread), on `grid` (`sll_board.HeldOutGrid`,
    default 0.5 deg x 1 deg over the hemisphere).

    `beams`: `[(theta, phi, polarization), ...]`; `beam_port[k]`: which aperture radiates beam k
    (default all 0); `cell_port`: `(N,)` aperture id of every cell (default all 0) -- a beam's main-lobe
    cone follows the extent of ITS aperture, so a sub-aperture board passes its port map. The reported
    `extent_m` is that of the apertures in `beam_port`."""
    from types import SimpleNamespace

    from .sll_board import HeldOutGrid, _scenario_eval

    n = board.n_cells
    cell_port = np.zeros(n, dtype=int) if cell_port is None else np.asarray(cell_port, dtype=int)
    beam_port = [0] * len(beams) if beam_port is None else list(beam_port)
    ports = [
        SimpleNamespace(
            port_id=int(p),
            x_m=board.x_m[cell_port == p],
            y_m=board.y_m[cell_port == p],
            owned=np.ones(int(np.sum(cell_port == p)), dtype=bool),
        )
        for p in sorted(set(cell_port.tolist()))
    ]
    problem = SimpleNamespace(
        ports=ports,
        period_m=board.period_m,
        polarization=board.polarization,
        beams=[tuple(b) for b in beams],
        beam_port=beam_port,
    )
    sc = SimpleNamespace(name=name, ports=tuple(sorted(set(beam_port))), beams=tuple(range(len(beams))))
    grid = grid or HeldOutGrid()
    return _scenario_eval(problem, sc, None, None, board.frequency_hz, grid, None, ff=ffcache.as_callable(excitation))
