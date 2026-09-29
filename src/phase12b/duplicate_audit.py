"""Phase 12B Task C: read-only duplicate-geometry audit. Does not rebuild the Phase 7 library or
re-run Phases 8C-11; every check here re-derives the *mesh-discretization* signature (never a
solve) either from the real, already-solved Phase 7 checkpoint or from Phases 9-11's own,
unmodified geometry constants.
"""

from dataclasses import dataclass

from .. import library_schema
from ..basis_functions import MeshSettings
from ..discretization import duplicate_geometry_groups, metal_mask_key
from ..geometry import Geometry, UnitCell
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..phase8c.data_source import DEFAULT_CHECKPOINT_JSON, load_full_database
from ..phase9 import local_periodicity as p9_lpa
from ..phase10 import environment_geometry as p10_eg
from ..phase11 import patterns as p11_patterns

# --------------------------------------------------------------------- Phase 7 library (32x32)

_LIB_CFG = library_schema.LibraryConfig.default()
LIBRARY_PERIOD_M = _LIB_CFG.period_m
LIBRARY_MESH = MeshSettings(**_LIB_CFG.mesh)
LIBRARY_SUBSTRATE = Substrate(
    name=_LIB_CFG.substrate["name"], eps_r=_LIB_CFG.substrate["eps_r"],
    tan_delta=_LIB_CFG.substrate["tan_delta"], thickness=_LIB_CFG.substrate["thickness_m"],
)
LIBRARY_LATTICE = RectangularLattice(LIBRARY_PERIOD_M, LIBRARY_PERIOD_M)

# --------------------------------------------------------------- Phases 9-11's own solve mesh
# Fixed throughout src/phase9/supercell_solver.py, src/phase10/*.py, src/phase11/reference_solver.py:
# mesh_cells_per_side / mesh_cells_per_subcell = 8 (dx = period/8 per sub-cell, regardless of
# supercell size n; floquet order does not affect the metal mask, so it is irrelevant here).
SUPERCELL_MESH_CELLS = 8
SUPERCELL_MESH = MeshSettings(SUPERCELL_MESH_CELLS, SUPERCELL_MESH_CELLS, 2 * SUPERCELL_MESH_CELLS)


def _unit_cell(l_mm, w_mm, period_m=LIBRARY_PERIOD_M, substrate=LIBRARY_SUBSTRATE):
    lattice = RectangularLattice(period_m, period_m) if period_m != LIBRARY_PERIOD_M else LIBRARY_LATTICE
    return UnitCell(lattice, Geometry("rectangle", l_mm * 1e-3, w_mm * 1e-3), substrate, Metal())


@dataclass(frozen=True)
class DuplicateReport:
    n_geometries: int
    n_distinct: int
    n_duplicate_groups: int
    groups: tuple


def phase7_library_duplicate_report():
    """Sec. 'Quantify the duplicates in the Phase 7 library': all 18x18=324 (L, W) geometries at
    the library's own locked 32x32 mesh, grouped by discretized metal mask (not by equal S
    values -- a stronger, purely-geometric check that holds for every frequency at once, since
    the mask does not depend on frequency)."""
    items = []
    for l in library_schema.L_VALUES_MM:
        for w in library_schema.W_VALUES_MM:
            items.append(((l, w), metal_mask_key(_unit_cell(l, w), LIBRARY_MESH)))
    groups = duplicate_geometry_groups(items)
    n_extra = sum(len(g) - 1 for g in groups)
    return DuplicateReport(len(items), len(items) - n_extra, len(groups), groups)


def phase7_delivered_lookup_duplicate_count(checkpoint_path=DEFAULT_CHECKPOINT_JSON, frequency_hz=12.5e9):
    """Cross-check against the brief's own claim ("169 distinct complex responses... at
    12.5 GHz"): counts distinct *solved* S11_xx complex values (rounded to solver precision) in
    the real checkpoint, independently of the geometric mask-based count above. Returns None if
    the checkpoint is not present."""
    import numpy as np

    try:
        db = load_full_database(checkpoint_path)
    except FileNotFoundError:
        return None
    freqs = db.frequencies
    idx = int(np.argmin(np.abs(freqs - frequency_hz)))
    values = set()
    for c in db.cells:
        if c.status != "ok":
            continue
        z = db.record(c.index, idx).component("S11_xx")
        values.add((round(z.real, 12), round(z.imag, 12)))
    return {"n_cells": sum(1 for c in db.cells if c.status == "ok"), "n_distinct_solved": len(values)}


# ------------------------------------------------------------------- Phases 9-11's own mesh

def phase9_strong_contrast_check():
    """The center (2.8, 2.8) mm vs. the 4 distinct strong-contrast corner values (1.0,1.0),
    (4.4,4.4), (1.0,4.4), (4.4,1.0) mm, at Phase 9's own 8-cell/side mesh: none should be
    mesh-identical to the center or to each other (they are meant to be a real, distinguishable
    strong perturbation)."""
    items = [("center", metal_mask_key(_unit_cell(*p9_lpa.CENTER_LW_MM), SUPERCELL_MESH))]
    for lw in sorted(set(p9_lpa.STRONG_NEIGHBORS_MM)):
        items.append((f"strong_{lw}", metal_mask_key(_unit_cell(*lw), SUPERCELL_MESH)))
    return duplicate_geometry_groups(items)


def phase10_center_duplicate_report():
    """Phase 10's 5 named centers (C1-C5), at Phase 9-11's own 8-cell/side mesh: which, if any,
    are mesh-identical to each other."""
    items = [(name, metal_mask_key(_unit_cell(*lw), SUPERCELL_MESH)) for name, lw in p10_eg.CENTER_CELLS_MM.items()]
    return duplicate_geometry_groups(items)


def phase10_holdout_leakage_check(holdout_center_names=("C5",)):
    """Sec. 'in particular any train/test leakage in Phase 10': does the held-out center's own
    isolated-cell mask collide with any *other* (training) center's mask, at the same mesh? A
    collision here would mean the blind test set is not actually geometrically distinct from
    training -- this checks the exact condition, not just 'are C1-C5 all pairwise distinct'."""
    groups = phase10_center_duplicate_report()
    holdout = set(holdout_center_names)
    leaks = [g for g in groups if (set(g) & holdout) and (set(g) - holdout)]
    return {"duplicate_groups": groups, "holdout_center_names": tuple(holdout_center_names), "leaks": tuple(leaks)}


def phase10_center_perturbation_duplicate_report(deltas_mm=(-0.6, -0.2, 0.0, 0.2, 0.6)):
    """For each center, the L-axis neighbor perturbations Phase 10 actually studies
    (weak=0.2mm, moderate=0.6mm, README/environment_geometry.py), W held at the center's own
    value: which perturbed points collide with another CENTER's perturbed points (not just within
    the same center's own family, which is not a leakage concern)."""
    items = []
    for name, (l, w) in p10_eg.CENTER_CELLS_MM.items():
        for d in deltas_mm:
            new_l = round(l + d, 2)
            if library_schema.L_VALUES_MM[0] <= new_l <= library_schema.L_VALUES_MM[-1]:
                items.append((f"{name}_L{d:+.1f}", metal_mask_key(_unit_cell(new_l, w), SUPERCELL_MESH)))
    groups = duplicate_geometry_groups(items)
    cross_center = tuple(g for g in groups if len({k.split("_L")[0] for k in g}) > 1)
    return {"all_groups": groups, "cross_center_groups": cross_center}


def phase11_pattern_endpoints_check():
    """Phase 11's checkerboard endpoints A_DEFAULT_MM=(1.4,1.4), B_DEFAULT_MM=(4.2,4.2): confirms
    they are not mesh-identical at the same mesh (they are meant to be the two, maximally
    different, checkerboard values)."""
    a = metal_mask_key(_unit_cell(*p11_patterns.A_DEFAULT_MM), SUPERCELL_MESH)
    b = metal_mask_key(_unit_cell(*p11_patterns.B_DEFAULT_MM), SUPERCELL_MESH)
    return a == b
