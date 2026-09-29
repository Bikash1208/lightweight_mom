"""Phase 8C Sec. 31-32: numeric reciprocity and geometric-symmetry checks against the real,
already-solved Phase 7 data (`phase8c.data_source.load_full_database`) -- genuine measurements,
not assumptions. `MOM_FORMULATION.md`/`CONVENTIONS.md` already claim `Z` is complex-symmetric at
normal incidence (hence `S = S^T`, reciprocity); this module checks that claim numerically on real
solved data rather than repeating it as an assertion.
"""

from dataclasses import dataclass

import numpy as np

from ..database import family_block


@dataclass(frozen=True)
class ReciprocityCheck:
    n_cells_checked: int
    max_abs_error_co: float       # max |S21_xx - S12_xx|
    max_abs_error_cross: float    # max |S21_yx - S12_xy| (transpose relation)
    reciprocal: bool
    tolerance: float
    note: str


def check_reciprocity(database, frequency_index, tolerance=1e-6):
    """S = S^T (block form: S12 = S21^T) implies S21_xx = S12_xx and S21_yx = S12_xy -- checked
    directly on every successful cell's actual solved 4x4 matrix at ``frequency_index``."""
    errors_co, errors_cross = [], []
    for cell in database.cells:
        if cell.status != "ok":
            continue
        record = database.record(cell.index, frequency_index)
        s21 = family_block(record.s, "S21")
        s12 = family_block(record.s, "S12")
        errors_co.append(abs(s21[0, 0] - s12[0, 0]))
        errors_cross.append(abs(s21[1, 0] - s12[0, 1]))
    max_co, max_cross = float(np.max(errors_co)), float(np.max(errors_cross))
    reciprocal = max_co <= tolerance and max_cross <= tolerance
    note = (
        f"max|S21xx-S12xx|={max_co:.3e}, max|S21yx-S12xy|={max_cross:.3e} over "
        f"{len(errors_co)} cells -- consistent with `MOM_FORMULATION.md`'s claim that Z is "
        "complex-symmetric at normal incidence (S=S^T), within solver numerical precision."
        if reciprocal else
        f"max|S21xx-S12xx|={max_co:.3e}, max|S21yx-S12xy|={max_cross:.3e} EXCEEDS tolerance "
        f"{tolerance:.1e} -- reciprocity does NOT hold to the expected precision; investigate "
        "numerical error, polarization basis, port convention, or a material-model mismatch "
        "before using this data (Sec. 31)."
    )
    return ReciprocityCheck(
        n_cells_checked=len(errors_co), max_abs_error_co=max_co, max_abs_error_cross=max_cross,
        reciprocal=reciprocal, tolerance=tolerance, note=note,
    )


@dataclass(frozen=True)
class SymmetryCheck:
    n_pairs_checked: int
    max_abs_error: float
    symmetric: bool
    tolerance: float
    note: str


def check_lw_swap_symmetry(database, frequency_index, component_family="S11", tolerance=1e-3):
    """Sec. 32: a rectangular patch's x-polarized response at (L, W) should equal its
    y-polarized response at (W, L) under the x<->y mirror symmetry of the lattice/mesh (a real,
    testable physical symmetry -- checked only where both cells actually exist in the sweep, never
    assumed for a pair the library did not simulate). A modest tolerance is used because the mesh
    itself is not exactly x<->y symmetric for a non-square mesh count in general (`cells_x ==
    cells_y == 32` here, so it should hold well, but ordinary staircase differences are still
    possible for asymmetric (L, W) mesh alignment)."""
    by_lw = {}
    for cell in database.cells:
        if cell.status != "ok":
            continue
        key = (round(cell.parameters["length_m"] * 1e3, 9), round(cell.parameters["width_m"] * 1e3, 9))
        by_lw[key] = cell

    errors = []
    checked = set()
    for (l_mm, w_mm), cell in by_lw.items():
        swapped = by_lw.get((w_mm, l_mm))
        if swapped is None or (w_mm, l_mm) in checked:
            continue
        checked.add((l_mm, w_mm))
        r1 = database.record(cell.index, frequency_index)
        r2 = database.record(swapped.index, frequency_index)
        xx = family_block(r1.s, component_family)[0, 0]
        yy = family_block(r2.s, component_family)[1, 1]
        errors.append(abs(xx - yy))

    if not errors:
        return SymmetryCheck(0, float("nan"), False, tolerance, "no (L,W)/(W,L) pairs found in this sweep to check")
    max_err = float(np.max(errors))
    symmetric = max_err <= tolerance
    note = (
        f"{component_family}_xx(L,W) vs {component_family}_yy(W,L): max abs error "
        f"{max_err:.3e} over {len(errors)} pairs "
        f"({'within' if symmetric else 'EXCEEDS'} tolerance {tolerance:.1e})."
    )
    return SymmetryCheck(len(errors), max_err, symmetric, tolerance, note)
