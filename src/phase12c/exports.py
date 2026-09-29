"""Phase 12C Task E: HFSS-ready exports for a small, deliberate set of grounded 8 mm cases.

Uses the Phase 6/9 export conventions unchanged (``cad_export.export_unit_cell_step``/
``export_supercell_step``, reference planes and time-convention note in the accompanying
``geometry.json``-style record) plus real MoM reference results (S-parameters for single cells,
dipole moments for supercells) for the SAME cases, so an eventual HFSS import has something to
compare against. No HFSS results are created or implied anywhere in this module.
"""

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .. import cad_export
from ..basis_functions import MeshSettings
from ..sparameters import compute_s_parameters
from . import design_point as dp
from .grounded_supercell import build_grounded_supercell
from .locality_study import L_FAMILY_MM, MESH_CELLS_PER_SUBCELL, floquet_for

# 5 single-cell cases spanning the measured phase curve (L_FAMILY_MM), including the near-resonance
# point (6.4 mm, Sec. B) -- a deliberate, small, documented set, not an exhaustive sweep.
SINGLE_CELL_L_MM = (1.6, 3.2, 4.8, 5.6, 6.4)


@dataclass(frozen=True)
class SingleCellExport:
    l_mm: float
    step_path: str
    s11_mag: float
    s11_phase_deg: float
    mesh_cells_per_side: int
    floquet_order: int


def export_single_cells(output_dir, l_values_mm=SINGLE_CELL_L_MM, mesh_cells_per_side=MESH_CELLS_PER_SUBCELL):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    floquet = floquet_for(mesh_cells_per_side)
    results = []
    for l_mm in l_values_mm:
        cell = dp.unit_cell_at(l_mm * 1e-3)
        mesh = MeshSettings(mesh_cells_per_side, mesh_cells_per_side, floquet)
        s = compute_s_parameters(cell, dp.DESIGN_FREQUENCY_HZ, mesh, dp.REFERENCE_PLANES,
                                  warn_grating=False, warn_surface_wave=False)
        step_path = output_dir / f"single_cell_L{l_mm:.1f}mm.step"
        cad_export.export_unit_cell_step(cell, step_path)
        results.append(SingleCellExport(l_mm, str(step_path), s.magnitude("S11_xx"), s.phase_deg("S11_xx"),
                                         mesh_cells_per_side, floquet))
    return tuple(results)


@dataclass(frozen=True)
class SupercellExport:
    name: str
    step_path: str
    grid_mm: tuple
    n: int
    n_metal_solids: int
    has_ground: bool


def export_supercells(output_dir, mesh_cells_per_subcell=MESH_CELLS_PER_SUBCELL):
    """Uniform (3x3, homogeneous control), smooth-gradient and strong-contrast 3x3 grounded
    supercells -- the brief's own required trio."""
    from .coupling_study import build_checkerboard_pattern_mm, build_gradient_pattern_mm

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    center_mm = L_FAMILY_MM[len(L_FAMILY_MM) // 2]
    lo_mm, hi_mm = L_FAMILY_MM[0], L_FAMILY_MM[-1]

    cases = {
        "uniform_3x3": [[(center_mm, center_mm)] * 3 for _ in range(3)],
        "smooth_gradient_3x3": build_gradient_pattern_mm(center_mm, lo_mm, hi_mm, 3),
        "strong_contrast_3x3": build_checkerboard_pattern_mm(center_mm, lo_mm, 3),
    }
    results = []
    for name, grid_mm in cases.items():
        grid_m = [[(l * 1e-3, w * 1e-3) for (l, w) in row] for row in grid_mm]
        supercell = build_grounded_supercell(3, dp.PERIOD_M, dp.PERIOD_M, dp.SUBSTRATE, dp.METAL, grid_m)
        step_path = output_dir / f"{name}.step"
        res = cad_export.export_supercell_step(supercell, step_path)
        results.append(SupercellExport(name, str(step_path), tuple(tuple(row) for row in grid_mm), 3,
                                        len(res.metals), res.ground is not None))
    return tuple(results)


def write_synthetic_hfss_reference_test(path):
    """Sec. E: "make sure the Phase 9 comparison tooling can ingest an HFSS CSV for these cases
    ... test it on a synthetic file only". Writes a clearly-labelled SYNTHETIC one-port CSV (no
    real HFSS run exists) purely to exercise ``phase9.hfss_reference.load_hfss_reference``'s
    one-port ingestion path end to end; the numbers below are made up for this test file and are
    never treated as, or compared against, a real reference anywhere in this project."""
    from ..phase9 import hfss_reference as hfss_mod

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        "cell_id,L_mm,W_mm,frequency_GHz,S11_real,S11_imag",
        "SYNTHETIC-G-01,6.0,6.0,12.5,-0.95,0.10",
        "SYNTHETIC-G-02,6.4,6.4,12.5,-0.70,0.60",
    ]
    path.write_text("\n".join(rows) + "\n")
    records = hfss_mod.load_hfss_reference(path)
    if records is None or any(r.s21 is not None for r in records):
        raise AssertionError("write_synthetic_hfss_reference_test: one-port ingestion did not round-trip correctly")
    return path, records


def write_export_manifest(single_cells, supercells, synthetic_hfss_path, output_path):
    manifest = {
        "single_cells": [asdict(r) for r in single_cells],
        "supercells": [asdict(r) for r in supercells],
        "synthetic_hfss_reference_test_csv": str(synthetic_hfss_path),
        "note": "No real HFSS results exist anywhere in this project; the synthetic CSV above is "
                "a test fixture for the one-port ingestion path only (Sec. E).",
    }
    Path(output_path).write_text(json.dumps(manifest, indent=2))
    return output_path
