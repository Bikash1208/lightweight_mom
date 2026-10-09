"""Phase 22: automatic coded-cell states. The user gives the substrate, cell size, frequency and the
number of bits; this module finds the cell geometries (empty + square patches) and their reflection
coefficient at the patch plane, with the project's own periodic MoM.

How a square's phase is found (and why it is not a plain mesh solve):
* The rooftop mesh is cell-centre based, so a symmetric patch is always an EVEN number of mesh cells
  wide: at dx = P/64 (0.1 mm) the effective side steps by 0.2 mm and a requested 5.85 mm silently becomes
  5.8 mm. A side on the 0.05 mm fabrication grid cannot be put on the mesh directly (dx = 0.025 mm
  would need ~10 GB of dense matrix).
* So the phase-vs-side curve is built from mesh-aligned solves only -- a cheap 0.2 mm-mesh scan, then
  three 0.1 mm-mesh solves around the target -- each corrected by a first-order Richardson extrapolation
  in the mesh size (v_inf = v(dx) + (v(dx) - v(2dx)), the Phase 12 behaviour near resonance), and
  interpolated by a quadratic. The side is snapped to `L_grid_mm` and the reported Gamma is that curve
  evaluated at the snapped side (labelled `GAMMA_INTERPOLATED`), never the value at the unsnapped
  optimum. The Richardson shift is reported as the mesh uncertainty.
* Solves are cached on disk by (substrate, period, frequency, mesh, side in mesh cells)."""

import hashlib
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..phase import principal_phase_deg, wrapped_difference_deg
from .two_state_cells import TwoStateCell, TwoStateChoice

COARSE_CELLS = 32  # dx = P/32 (0.2 mm for the 6.4 mm reference cell)
FINE_CELLS = 64  # dx = P/64 (0.1 mm)
SCAN_MIN_FRACTION = 0.3
SCAN_MAX_FRACTION = 0.97
NEAR_RESONANCE_SENSITIVITY_DEG_PER_0P01MM = 1.0
CACHE_VERSION = "auto_states_v1"

LABEL_NEAR_RESONANCE = "NEAR_RESONANCE"
LABEL_UNREACHABLE_LEVEL = "UNREACHABLE_LEVEL"
LABEL_GAMMA_INTERPOLATED = "GAMMA_INTERPOLATED"
LABEL_EXACT_EMPTY = "EXACT_GROUNDED_SLAB"


@dataclass(frozen=True)
class AutoState:
    level: int
    label: str
    geometry: dict  # {"kind": "empty"} or {"kind": "rectangle", "L_mm", "W_mm"}
    gamma: complex  # at the patch plane
    target_phase_deg: float  # absolute (empty's phase + level * 360 / 2**n_bits)
    phase_error_deg: float  # |achieved - target|
    sensitivity_deg_per_0p01mm: "float | None"
    mesh_uncertainty_deg: "float | None"  # |Richardson shift| of the bracketing solves
    labels: tuple = ()

    @property
    def phase_deg(self):
        return float(np.degrees(np.angle(self.gamma)))

    @property
    def magnitude_db(self):
        return float(20.0 * np.log10(abs(self.gamma)))


@dataclass(frozen=True)
class AutoStates:
    n_bits: int
    states: "tuple[AutoState, ...]"
    frequency_hz: float
    period_mm: float
    solve_seconds: float
    n_solves: int  # MoM solves actually run (cache misses)
    n_cached: int  # solves served from the cache

    def two_state_choice(self):
        """The `TwoStateChoice` the 1-bit coded synthesis consumes."""
        if self.n_bits != 1:
            raise ValueError("two_state_choice is for 1-bit states")
        a, b = self.cells()
        sep = abs(wrapped_difference_deg(principal_phase_deg(a.value), principal_phase_deg(b.value)))
        return TwoStateChoice(state_0=a, state_1=b, achieved_separation_deg=float(sep))

    def cells(self):
        """The states as `TwoStateCell`s (patch-plane values, geometry for export)."""
        return [_as_cell(s) for s in self.states]

    def library_candidates(self, library_id="auto_states"):
        """`(cell_id, L_mm, W_mm, library_id, value)` tuples for `synthesize_coded_n_bit`: each state's
        value is its Gamma with the empty cell's phase removed, so level k sits at k * 360 / 2**n_bits
        (a constant phase offset does not change the far field)."""
        ref = np.exp(-1j * np.angle(self.states[0].gamma))
        out = []
        for s in self.states:
            geom = s.geometry
            length, width = (0.0, 0.0) if geom["kind"] == "empty" else (geom["L_mm"], geom["W_mm"])
            out.append((f"auto_level_{s.level}", length, width, library_id, complex(s.gamma * ref)))
        return out


def _as_cell(state):
    return TwoStateCell(state.label, f"auto_level_{state.level}", complex(state.gamma), state.geometry)


# ----------------------------------------------------------------------------- solves

def _substrate_key(substrate):
    return (float(substrate.eps_r), float(substrate.tan_delta), float(substrate.thickness))


def effective_cells(n_cells_requested_side_mm, n, period_m, substrate=None):
    """How many mesh cells wide a centred square of `n_cells_requested_side_mm` really is on an n x n
    mesh (cell-centre mask: always even). Computed from the solver's own mask, not assumed."""
    from ..basis_functions import RooftopBasis
    from ..geometry import Geometry, UnitCell
    from ..materials import Metal, Substrate
    from ..periodic_bc import RectangularLattice

    sub = substrate or Substrate("probe", 2.2, 0.0, 1.0e-3)
    uc = UnitCell(
        RectangularLattice(period_m, period_m),
        Geometry("rectangle", n_cells_requested_side_mm * 1e-3, n_cells_requested_side_mm * 1e-3),
        sub,
        Metal(),
        backing="pec",
    )
    mask = RooftopBasis(uc, n, n).metal_mask
    return int(mask[mask.shape[0] // 2].sum())


def mom_square_gamma(side_cells, n, substrate, frequency_hz, period_m):
    """One real MoM solve: Gamma at the patch plane of a centred square `side_cells` mesh cells wide on
    an n x n mesh (PEC ground, normal incidence)."""
    from ..basis_functions import MeshSettings
    from ..excitation import ReferencePlanes
    from ..geometry import Geometry, UnitCell
    from ..materials import Metal
    from ..periodic_bc import RectangularLattice
    from ..sparameters import compute_s_parameters

    side_m = side_cells * period_m / n
    uc = UnitCell(RectangularLattice(period_m, period_m), Geometry("rectangle", side_m, side_m), substrate, Metal(), backing="pec")
    planes = ReferencePlanes(top_air_thickness=0.0, bottom_air_thickness=None, substrate_thickness=substrate.thickness)
    sp = compute_s_parameters(uc, frequency_hz, MeshSettings(n, n, n), planes, warn_grating=False)
    return complex(sp.s11[0, 0])


def empty_state_gamma(substrate, frequency_hz):
    """State 0: the exact grounded slab at the patch plane (`scattering_from_fields` with zero current)."""
    from ..excitation import ReferencePlanes
    from ..green_function import LayeredGreenFunction
    from ..sparameters import scattering_from_fields

    planes = ReferencePlanes(top_air_thickness=0.0, bottom_air_thickness=None, substrate_thickness=substrate.thickness)
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    return complex(scattering_from_fields(green, planes, np.zeros((1, 2, 2), complex)).s11[0, 0])


class _Solves:
    """A cached MoM-solve front end that counts misses (`n_solves`) and hits (`n_cached`)."""

    def __init__(self, substrate, frequency_hz, period_m, cache_dir=None, solver=None, progress=None):
        self.substrate, self.frequency_hz, self.period_m = substrate, frequency_hz, period_m
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.solver = solver or mom_square_gamma
        self.progress = progress
        self.n_solves = 0
        self.n_cached = 0
        self._memo = {}

    def _path(self, side_cells, n):
        key = json.dumps(
            [CACHE_VERSION, *_substrate_key(self.substrate), self.period_m, self.frequency_hz, n, side_cells]
        )
        return self.cache_dir / (hashlib.sha1(key.encode()).hexdigest() + ".json")

    def gamma(self, side_cells, n):
        memo_key = (side_cells, n)
        if memo_key in self._memo:
            return self._memo[memo_key]
        value = None
        if self.cache_dir:
            path = self._path(side_cells, n)
            if path.is_file():
                try:
                    d = json.loads(path.read_text())
                    value = complex(d["re"], d["im"])
                    self.n_cached += 1
                except (ValueError, KeyError, OSError):
                    value = None
        if value is None:
            if self.progress:
                self.progress(f"MoM solve: square {side_cells * self.period_m / n * 1e3:.2f} mm on the {n}-cell mesh")
            value = complex(self.solver(side_cells, n, self.substrate, self.frequency_hz, self.period_m))
            self.n_solves += 1
            if self.cache_dir:
                self._path(side_cells, n).write_text(json.dumps({"re": value.real, "im": value.imag}))
        self._memo[memo_key] = value
        return value


# ----------------------------------------------------------------------------- the phase-vs-side curve

def _unwrapped_deg(gammas):
    return np.degrees(np.unwrap(np.angle(np.asarray(gammas))))


@dataclass(eq=False)
class _Curve:
    """Richardson-corrected Gamma at three mesh-aligned sides, with a quadratic through them."""

    sides_mm: np.ndarray
    phase_deg: np.ndarray  # unwrapped, corrected
    log_mag: np.ndarray
    shift_deg: float  # mean |Richardson shift| at the end points

    def __post_init__(self):
        x = self.sides_mm - self.sides_mm[1]
        self._phase_poly = np.polyfit(x, self.phase_deg, 2)
        self._mag_poly = np.polyfit(x, self.log_mag, 2)

    def align_branch(self, reference_phase_deg):
        """Move the whole curve by a multiple of 360 degrees so its first point is nearest
        `reference_phase_deg` (the scan's unwrapped phase at the same side)."""
        shift = 360.0 * round((reference_phase_deg - self.phase_deg[0]) / 360.0)
        if shift:
            self.phase_deg = self.phase_deg + shift
            self.__post_init__()

    def phase_at(self, side_mm):
        return float(np.polyval(self._phase_poly, side_mm - self.sides_mm[1]))

    def mag_at(self, side_mm):
        return float(np.exp(np.polyval(self._mag_poly, side_mm - self.sides_mm[1])))

    def slope_deg_per_mm(self, side_mm):
        return float(np.polyval(np.polyder(self._phase_poly), side_mm - self.sides_mm[1]))

    def gamma_at(self, side_mm):
        return self.mag_at(side_mm) * np.exp(1j * np.radians(self.phase_at(side_mm)))


def _build_curve(solves, la_mm, period_m):
    """`la_mm`: a multiple of two coarse cells; the curve spans [la, la + 2 coarse cells]."""
    dx_c = period_m / COARSE_CELLS * 1e3
    dx_f = period_m / FINE_CELLS * 1e3
    m_a = int(round(la_mm / dx_f))  # fine cells at La
    sides = np.array([la_mm, la_mm + dx_c, la_mm + 2 * dx_c])
    fine = [solves.gamma(m_a + k * 2, FINE_CELLS) for k in range(3)]
    coarse_a = solves.gamma(int(round(la_mm / dx_c)), COARSE_CELLS)
    coarse_b = solves.gamma(int(round(la_mm / dx_c)) + 2, COARSE_CELLS)
    f_ph = _unwrapped_deg(fine)
    c_ph = _unwrapped_deg([coarse_a, coarse_b])
    # keep each coarse phase on the fine solve's 360-degree branch
    shifts = []
    corrected = f_ph.copy()
    for idx, (c_val, f_val) in zip((0, 2), ((c_ph[0], f_ph[0]), (c_ph[1], f_ph[2])), strict=True):
        c_val = c_val + 360.0 * round((f_val - c_val) / 360.0)
        shift = f_val - c_val
        shifts.append(shift)
        corrected[idx] = f_val + shift
    corrected[1] = f_ph[1] + 0.5 * (shifts[0] + shifts[1])
    mags = np.array([abs(g) for g in fine])
    mag_c = [abs(coarse_a), abs(coarse_b)]
    mag_corr = mags.copy()
    mag_corr[0] = mags[0] + (mags[0] - mag_c[0])
    mag_corr[2] = mags[2] + (mags[2] - mag_c[1])
    mag_corr[1] = mags[1] + 0.5 * ((mags[0] - mag_c[0]) + (mags[2] - mag_c[1]))
    mag_corr = np.clip(mag_corr, 1e-6, 1.0)
    return _Curve(sides, corrected, np.log(mag_corr), float(np.mean(np.abs(shifts))))


def _scan(solves, period_m):
    """Coarse scan: centred squares 0.4 mm apart (even coarse-cell widths) from 0.3 P to 0.97 P.
    Returns `(sides_mm, unwrapped_phase_deg)` ordered by increasing side."""
    dx_c = period_m / COARSE_CELLS * 1e3
    lo = int(np.ceil(SCAN_MIN_FRACTION * period_m * 1e3 / dx_c / 2.0)) * 2
    hi = int(np.floor(SCAN_MAX_FRACTION * period_m * 1e3 / dx_c / 2.0)) * 2
    cells = list(range(lo, hi + 1, 2))
    sides = np.array([c * dx_c for c in cells])
    gammas = [solves.gamma(c, COARSE_CELLS) for c in cells]
    return sides, _unwrapped_deg(gammas), gammas


def _round_to_grid(value_mm, grid_mm):
    return round(value_mm / grid_mm) * grid_mm


def find_square_for_phase(
    target_phase_deg,
    substrate,
    frequency_hz,
    period_m,
    L_grid_mm=0.05,
    cache_dir=None,
    progress=None,
    solver=None,
    solves=None,
):
    """The square whose patch-plane phase is closest to `target_phase_deg` (absolute): returns
    `(side_mm, gamma, info)`, `info` = {sensitivity, mesh_uncertainty, labels, phase_error}."""
    solves = solves or _Solves(substrate, frequency_hz, period_m, cache_dir, solver, progress)
    dx_c = period_m / COARSE_CELLS * 1e3
    sides, phases, gammas = _scan(solves, period_m)
    # find a bracket of the (monotonic) curve around some branch of the target
    best = None
    for k in range(-3, 4):
        t = target_phase_deg + 360.0 * k
        for i in range(len(sides) - 1):
            a, b = phases[i], phases[i + 1]
            if (a - t) * (b - t) <= 0.0 and a != b:
                best = (i, t)
                break
        if best:
            break
    if best is None:  # unreachable with squares: the closest scanned end point, Richardson-corrected
        errs = [abs(wrapped_difference_deg(float(np.degrees(np.angle(g))), target_phase_deg)) for g in gammas]
        i = int(np.argmin(errs))
        la = float(sides[min(i, len(sides) - 2)])
        curve = _build_curve(solves, la, period_m)
        grid_sides = np.arange(la, la + 2 * dx_c + 1e-9, L_grid_mm)
        errs_fine = [
            abs(wrapped_difference_deg(float(np.degrees(np.angle(curve.gamma_at(s)))), target_phase_deg))
            for s in grid_sides
        ]
        side = float(grid_sides[int(np.argmin(errs_fine))])
        labels = (LABEL_UNREACHABLE_LEVEL, LABEL_GAMMA_INTERPOLATED)
    else:
        i, t = best
        la = float(sides[i])
        curve = _build_curve(solves, la, period_m)
        curve.align_branch(phases[i])
        # solve curve.phase_at(L) = t on [la, la + 2 dx_c] (the quadratic is monotonic on the bracket)
        grid = np.linspace(curve.sides_mm[0], curve.sides_mm[2], 4001)
        vals = np.array([curve.phase_at(s) for s in grid]) - t
        crossings = np.where(vals[:-1] * vals[1:] <= 0.0)[0]
        idx = int(crossings[0]) if len(crossings) else int(np.argmin(np.abs(vals)))
        side_star = float(grid[idx])
        side = float(np.clip(_round_to_grid(side_star, L_grid_mm), curve.sides_mm[0], curve.sides_mm[2]))
        labels = (LABEL_GAMMA_INTERPOLATED,)
    gamma = curve.gamma_at(side)
    sens = abs(curve.slope_deg_per_mm(side)) * 0.01
    if sens > NEAR_RESONANCE_SENSITIVITY_DEG_PER_0P01MM:
        labels = (LABEL_NEAR_RESONANCE, *labels)
    err = abs(wrapped_difference_deg(float(np.degrees(np.angle(gamma))), target_phase_deg))
    return round(side, 6), complex(gamma), {
        "sensitivity_deg_per_0p01mm": float(sens),
        "mesh_uncertainty_deg": curve.shift_deg,
        "labels": labels,
        "phase_error_deg": float(err),
    }


def gamma_for_square(side_mm, substrate, frequency_hz, period_m, cache_dir=None, progress=None, solver=None, solves=None):
    """The user's own square (the "fix L" option): its Gamma from the same corrected, interpolated curve.
    Returns `(gamma, info)` like `find_square_for_phase`."""
    solves = solves or _Solves(substrate, frequency_hz, period_m, cache_dir, solver, progress)
    dx_c = period_m / COARSE_CELLS * 1e3
    if not (dx_c * 2 <= side_mm <= period_m * 1e3 * SCAN_MAX_FRACTION + dx_c):
        raise ValueError(f"square side {side_mm} mm is outside the supported range for a {period_m * 1e3:g} mm cell")
    la = np.floor(side_mm / (2 * dx_c)) * 2 * dx_c
    la = min(la, (np.floor(SCAN_MAX_FRACTION * period_m * 1e3 / (2 * dx_c)) - 1) * 2 * dx_c)
    curve = _build_curve(solves, float(la), period_m)
    gamma = curve.gamma_at(side_mm)
    sens = abs(curve.slope_deg_per_mm(side_mm)) * 0.01
    labels = (LABEL_GAMMA_INTERPOLATED,)
    if sens > NEAR_RESONANCE_SENSITIVITY_DEG_PER_0P01MM:
        labels = (LABEL_NEAR_RESONANCE, *labels)
    return complex(gamma), {
        "sensitivity_deg_per_0p01mm": float(sens),
        "mesh_uncertainty_deg": curve.shift_deg,
        "labels": labels,
    }


def auto_states(
    n_bits,
    substrate,
    frequency_hz,
    period_m,
    L_grid_mm=0.05,
    fixed_L_mm=None,
    cache_dir=None,
    progress=None,
    solver=None,
):
    """Automatic coded states. 1-bit: empty + the square closest to empty + 180 deg. 2-/3-bit: 4/8 states
    at empty + k * 360 / 2**n_bits degrees (level 0 is the empty cell); a level no square can reach gets
    the closest square and an `UNREACHABLE_LEVEL` label with its phase error. `fixed_L_mm` (1-bit only)
    uses that square for state 1 instead of searching."""
    if n_bits < 1:
        raise ValueError(f"n_bits must be >= 1, got {n_bits!r}")
    if fixed_L_mm is not None and n_bits != 1:
        raise ValueError("fixed_L_mm is for 1-bit designs")
    t0 = time.time()
    solves = _Solves(substrate, frequency_hz, period_m, cache_dir, solver, progress)
    g_empty = empty_state_gamma(substrate, frequency_hz)
    phase_empty = float(np.degrees(np.angle(g_empty)))
    states = [
        AutoState(0, "empty (bare substrate)", {"kind": "empty"}, g_empty, phase_empty, 0.0, None, None, (LABEL_EXACT_EMPTY,))
    ]
    step = 360.0 / 2**n_bits
    for level in range(1, 2**n_bits):
        target = phase_empty + level * step
        if fixed_L_mm is not None:
            gamma, info = gamma_for_square(fixed_L_mm, substrate, frequency_hz, period_m, solves=solves)
            side = float(fixed_L_mm)
            info["phase_error_deg"] = abs(wrapped_difference_deg(float(np.degrees(np.angle(gamma))), target))
        else:
            side, gamma, info = find_square_for_phase(target, substrate, frequency_hz, period_m, L_grid_mm, solves=solves)
        states.append(
            AutoState(
                level,
                f"square patch {side:g} mm",
                {"kind": "rectangle", "L_mm": side, "W_mm": side},
                gamma,
                target,
                info["phase_error_deg"],
                info["sensitivity_deg_per_0p01mm"],
                info["mesh_uncertainty_deg"],
                tuple(info["labels"]),
            )
        )
    return AutoStates(n_bits, tuple(states), frequency_hz, period_m * 1e3, time.time() - t0, solves.n_solves, solves.n_cached)
