"""Phase 12A "Physics (verify, don't trust)" and "Validation > Analytic grounded slab" checks,
as *measured, reported* results (not just hidden unit tests): the bare (unpatterned) grounded
slab's reflection coefficient against the textbook closed form, and the TM0 surface-wave-pole
proximity for the project's two periods actually used anywhere (4.8 mm, the Phase 7 library;
8 mm, the original Phase 1-6 examples), swept across the 10-15 GHz band at normal incidence.
"""

from dataclasses import dataclass

import numpy as np

from ..constants import ETA0
from ..excitation import ReferencePlanes
from ..geometry import Geometry, UnitCell
from ..green_function import LayeredGreenFunction
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..sparameters import stack_s_parameters, surface_wave_proximity

FREQUENCIES_HZ = (10.0e9, 11.25e9, 12.5e9, 13.75e9, 15.0e9)
PERIODS_MM = (4.8, 8.0)


def textbook_grounded_slab_gamma(f_hz, eps_r, tan_delta, h_m):
    """Gamma = (j*eta_s*tan(k_s*h) - eta0) / (j*eta_s*tan(k_s*h) + eta0), eta_s = eta0/sqrt(eps_r)
    (normal incidence; the brief's own closed form for the validation check), computed
    independently of ``LayeredGreenFunction`` (plain complex arithmetic, no TL/sigma machinery)."""
    eps_c = eps_r * (1.0 - 1j * tan_delta)
    k0 = 2.0 * np.pi * f_hz / 299792458.0
    k_s = k0 * np.sqrt(eps_c)
    eta_s = ETA0 / np.sqrt(eps_c)
    num = 1j * eta_s * np.tan(k_s * h_m) - ETA0
    den = 1j * eta_s * np.tan(k_s * h_m) + ETA0
    return num / den


@dataclass(frozen=True)
class SlabCheckPoint:
    frequency_hz: float
    eps_r: float
    tan_delta: float
    thickness_m: float
    gamma_mine: complex
    gamma_textbook: complex
    max_error: float


def bare_grounded_slab_checks(eps_r=2.2, thickness_m=1.575e-3, frequencies_hz=FREQUENCIES_HZ):
    """The bare (unpatterned) grounded-slab reflection vs the textbook closed form, lossless and
    lossy, across the band. Every ``max_error`` here should be at the solver's own floating-point
    noise floor (no fitting, no tolerance is "tuned" to pass)."""
    points = []
    for tan_delta in (0.0, 0.0009):
        for f in frequencies_hz:
            sub = Substrate("x", eps_r, tan_delta, thickness_m)
            green = LayeredGreenFunction(sub, f, backing="pec")
            planes = ReferencePlanes(0.0, None, thickness_m)  # planes on the substrate faces
            s = stack_s_parameters(green, planes)
            mine = s.component("S11_xx")
            textbook = textbook_grounded_slab_gamma(f, eps_r, tan_delta, thickness_m)
            points.append(SlabCheckPoint(f, eps_r, tan_delta, thickness_m, mine, textbook, abs(mine - textbook)))
    return tuple(points)


def vanishing_thickness_limit(eps_r=2.2, f_hz=12.5e9, thickness_m=1e-9):
    """Gamma -> -1 as h -> 0 (the ground plane short-circuits the pattern plane directly)."""
    sub = Substrate("x", eps_r, 0.0, thickness_m)
    green = LayeredGreenFunction(sub, f_hz, backing="pec")
    planes = ReferencePlanes(0.0, None, thickness_m)
    s = stack_s_parameters(green, planes)
    return complex(s.component("S11_xx"))


@dataclass(frozen=True)
class SurfaceWavePoint:
    period_mm: float
    frequency_hz: float
    k0: float
    k_substrate_re: float
    any_in_band: bool
    any_near_band: bool
    nearest_margin_over_k0: float


def surface_wave_margin_table(
    periods_mm=PERIODS_MM,
    frequencies_hz=FREQUENCIES_HZ,
    eps_r=2.2,
    tan_delta=0.0009,
    thickness_m=1.575e-3,
    floquet_order=32,
):
    """Sec. "Physics": for the default 4.8 mm and 8 mm periods at 10-15 GHz, whether any Floquet
    harmonic used in the impedance-matrix assembly falls in (or near) the grounded slab's TM0
    surface-wave interval, at normal incidence -- the margin is measured, never assumed."""
    points = []
    for p_mm in periods_mm:
        lattice = RectangularLattice(p_mm * 1e-3, p_mm * 1e-3)
        geometry = Geometry("rectangle", p_mm * 1e-3 * 0.5, p_mm * 1e-3 * 0.5)
        uc = UnitCell(lattice, geometry, Substrate("x", eps_r, tan_delta, thickness_m), Metal(), backing="pec")
        for f in frequencies_hz:
            green = LayeredGreenFunction(uc.substrate, f, backing="pec")
            r = surface_wave_proximity(uc, green, floquet_order)
            points.append(
                SurfaceWavePoint(
                    p_mm,
                    f,
                    r["k0"],
                    r["k_substrate_re"],
                    r["any_in_band"],
                    r["any_near_band"],
                    r["nearest_margin_over_k0"],
                )
            )
    return tuple(points)
