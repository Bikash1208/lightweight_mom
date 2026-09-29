"""Sec. 15-25: the heterogeneous-supercell local-periodicity-approximation (LPA) test suite.

Every neighbor geometry below is a fixed, explicit, deterministic list (Sec. 17: "do not choose
them randomly without a reproducible seed and record" -- this module simply never uses randomness,
which is stronger than a recorded seed). All values stay within the validated Phase 7 sweep range
(1.0-4.4 mm).
"""

from dataclasses import dataclass

from ..phase import principal_phase_deg, wrapped_difference_deg
from .supercell_geometry import build_supercell
from .supercell_solver import solve_isolated_cell, solve_supercell_center

CENTER_LW_MM = (2.8, 2.8)

# Sec. 17: weak contrast -- a uniform, small (0.2 mm) perturbation from the centre geometry.
WEAK_NEIGHBOR_MM = (3.0, 3.0)

# Sec. 17: moderate contrast -- 8 distinct values spanning ~1.8-3.8 mm (row-major, skipping centre).
MODERATE_NEIGHBORS_MM = (
    (1.8, 1.8), (2.2, 3.4), (3.8, 2.0),
    (3.0, 1.8),             (2.0, 3.8),
    (3.8, 3.8), (1.8, 2.6), (2.6, 3.0),
)

# Sec. 17: strong contrast -- the sweep's own extremes (1.0 mm / 4.4 mm) in every combination.
STRONG_NEIGHBORS_MM = (
    (1.0, 1.0), (4.4, 4.4), (1.0, 4.4),
    (4.4, 1.0),             (1.0, 1.0),
    (4.4, 4.4), (1.0, 4.4), (4.4, 1.0),
)

# Sec. 18: symmetric A-B-A / B-C-B / A-B-A arrangement.
SYMMETRIC_A_MM = (1.2, 1.2)
SYMMETRIC_B_MM = (3.8, 3.8)

# Sec. 19: the 5 x 5 outer ring (the 16 positions at Chebyshev distance 2 from centre), tiled from
# a fixed 16-value list (cycles MODERATE_NEIGHBORS_MM/STRONG_NEIGHBORS_MM as appropriate).
NOISE_FLOOR_RELATIVE = 1e-6   # Sec. 24: not an engineering choice -- measured solver self-
                              # consistency floor (the homogeneous-supercell exact-match test
                              # reproduces the isolated cell to ~1e-14 at this mesh/Floquet
                              # resolution; 1e-6 is a generous margin above that floor, not a
                              # judgement about what error is acceptable for synthesis).


def _mm_to_m(pairs):
    return tuple((l * 1e-3, w * 1e-3) for l, w in pairs)


def build_3x3_grid(kind, center_lw_mm=CENTER_LW_MM):
    """Sec. 17-18: one named 3x3 neighbor arrangement, row-major, centre at [1][1]."""
    if kind == "weak":
        neighbors = (WEAK_NEIGHBOR_MM,) * 8
    elif kind == "moderate":
        neighbors = MODERATE_NEIGHBORS_MM
    elif kind == "strong":
        neighbors = STRONG_NEIGHBORS_MM
    elif kind == "symmetric":
        # A B A / B C B / A B A, read row-major with the centre C excluded: A,B,A, B,B, A,B,A
        neighbors = (SYMMETRIC_A_MM, SYMMETRIC_B_MM, SYMMETRIC_A_MM,
                     SYMMETRIC_B_MM, SYMMETRIC_B_MM,
                     SYMMETRIC_A_MM, SYMMETRIC_B_MM, SYMMETRIC_A_MM)
    else:
        raise ValueError(f"unknown 3x3 kind {kind!r}")
    flat = list(neighbors)
    flat_mm = flat[:4] + [center_lw_mm] + flat[4:]
    grid_mm = [flat_mm[0:3], flat_mm[3:6], flat_mm[6:9]]
    return [[_mm_to_m([pair])[0] for pair in row] for row in grid_mm]


def extend_to_5x5(inner_3x3_m, outer_ring_mm):
    """Sec. 19, 22: places ``inner_3x3_m`` (metres) unchanged at the centre of a 5x5 grid (so the
    immediate-neighbor environment is *identical* to the corresponding 3x3 test -- the only
    difference is the newly-added outer ring, letting Sec. 22's locality comparison isolate
    exactly that effect) and fills the 16 outer positions from ``outer_ring_mm`` (mm, cycled)."""
    outer_m = _mm_to_m(outer_ring_mm)
    grid = [[None] * 5 for _ in range(5)]
    for i in range(3):
        for j in range(3):
            grid[i + 1][j + 1] = inner_3x3_m[i][j]
    k = 0
    for i in range(5):
        for j in range(5):
            if grid[i][j] is None:
                grid[i][j] = outer_m[k % len(outer_m)]
                k += 1
    return grid


@dataclass(frozen=True)
class LpaResult:
    name: str
    n: int
    isolated_px: complex
    embedded_px: complex
    complex_error: float
    magnitude_error: float
    phase_error_deg: float
    relative_complex_error: float
    n_unknowns_isolated: int
    n_unknowns_embedded: int


def compute_lpa_error(isolated_dipole, embedded_dipole):
    """Sec. 21: complex, magnitude and (wrapped) phase deviation between the embedded and
    isolated centre-cell dipole moments."""
    iso, emb = isolated_dipole.px, embedded_dipole.px
    complex_error = abs(emb - iso)
    magnitude_error = abs(abs(emb) - abs(iso))
    phase_error = float(abs(wrapped_difference_deg(
        float(principal_phase_deg(emb)), float(principal_phase_deg(iso)),
    )))
    relative = complex_error / abs(iso) if abs(iso) > 0 else float("inf")
    return complex_error, magnitude_error, phase_error, relative


def run_local_periodicity_suite(substrate, metal, period_m, frequency_hz,
                                 mesh_cells_per_subcell=8, floquet_factor=2,
                                 center_lw_mm=CENTER_LW_MM):
    """Sec. 25's required matrix: 3x3 weak/moderate/strong/symmetric, 5x5 moderate/strong."""
    center_lw_m = tuple(v * 1e-3 for v in center_lw_mm)
    isolated = solve_isolated_cell(
        period_m, period_m, center_lw_m[0], center_lw_m[1], substrate, metal, frequency_hz,
        mesh_cells_per_subcell, floquet_factor * mesh_cells_per_subcell,
    )

    cases = []
    for kind in ("weak", "moderate", "strong", "symmetric"):
        grid = build_3x3_grid(kind, center_lw_mm)
        cases.append((f"3x3_{kind}", 3, grid))

    moderate_3x3 = build_3x3_grid("moderate", center_lw_mm)
    strong_3x3 = build_3x3_grid("strong", center_lw_mm)
    cases.append(("5x5_moderate", 5, extend_to_5x5(moderate_3x3, MODERATE_NEIGHBORS_MM)))
    cases.append(("5x5_strong", 5, extend_to_5x5(strong_3x3, STRONG_NEIGHBORS_MM)))

    results = []
    for name, n, grid in cases:
        supercell = build_supercell(n, period_m, period_m, substrate, metal, grid)
        embedded = solve_supercell_center(
            supercell, center_lw_m, frequency_hz, mesh_cells_per_subcell,
            floquet_factor * mesh_cells_per_subcell * n,
        )
        complex_error, magnitude_error, phase_error, relative = compute_lpa_error(isolated.dipole, embedded.dipole)
        results.append(LpaResult(
            name=name, n=n, isolated_px=isolated.dipole.px, embedded_px=embedded.dipole.px,
            complex_error=complex_error, magnitude_error=magnitude_error, phase_error_deg=phase_error,
            relative_complex_error=relative, n_unknowns_isolated=isolated.n_unknowns,
            n_unknowns_embedded=embedded.n_unknowns,
        ))
    return isolated, results


# Sec. 23's INVALIDATED trigger: not a percentage chosen for engineering convenience, but a
# self-referential, dimensionless comparison against the quantity's OWN scale -- "the embedded
# response differs from the isolated one by at least half its own magnitude" is substantial by
# construction, independent of wherever an application-specific tolerance would eventually be
# drawn (Sec. 24's remedy: report the measured errors, and only call a deviation 'substantial'
# when it is substantial relative to the signal itself, never a fixed external percentage).
SUBSTANTIAL_DEVIATION_RELATIVE = 0.5


def classify_local_periodicity(results, noise_floor=NOISE_FLOOR_RELATIVE,
                                substantial=SUBSTANTIAL_DEVIATION_RELATIVE):
    """Sec. 23-24: no project document or cited scientific basis defines an acceptable LPA error
    tolerance for the intended holographic synthesis task, so this function never invents an
    engineering percentage (Sec. 24's explicit prohibition on '<10%=valid'-style thresholds).
    Two non-arbitrary criteria are used instead:

    - `"VALIDATED"`: every tested contrast level's deviation is within the solver's own measured
      numerical self-consistency floor (``noise_floor``) -- indistinguishable from the periodic
      assumption, a solver-precision fact, not a judgement call.
    - `"INVALIDATED"`: at least one tested, physically-realistic contrast level (within the
      validated 1.0-4.4 mm sweep range) shows a deviation at least as large as ``substantial``
      times the signal's own magnitude -- self-evidently a failure of the periodic assumption in
      that regime, whatever an eventual application tolerance turns out to be.
    - `"UNVALIDATED"`: neither of the above -- deviations are measurably nonzero (so `VALIDATED`
      cannot be claimed) but not self-evidently substantial either, and no sourced tolerance
      exists to adjudicate the remaining cases (Sec. 24's own remedy).

    Every measured error is reported in full regardless of which label is returned (Sec. 24:
    'report the measured errors first')."""
    if all(r.relative_complex_error <= noise_floor for r in results):
        return "VALIDATED"
    if any(r.relative_complex_error >= substantial for r in results):
        return "INVALIDATED"
    return "UNVALIDATED"
