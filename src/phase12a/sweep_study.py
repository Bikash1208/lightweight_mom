"""Phase 12A physics-sanity sweep: a grounded rectangular-PEC-patch sweep on the **same** Rogers
RT/Duroid 5880, 1.575 mm, 4.8 mm period, mesh and Floquet order Phase 7 used (``library_schema``,
unmodified) -- so the only thing that differs from the Phase 7 library is the ground plane
itself, making the "grounded vs ungrounded" comparison a fair one.

Not a production library (the brief's own words): only the diagonal L = W of Phase 7's 18x18
sweep, at a handful of frequencies, and no classification/database export -- direct
``compute_s_parameters`` calls, reported as measured numbers, plots and a written comparison.
"""

from dataclasses import dataclass

import numpy as np

from .. import library_schema
from ..basis_functions import MeshSettings
from ..excitation import ReferencePlanes
from ..geometry import Geometry, UnitCell
from ..materials import Metal, Substrate
from ..periodic_bc import RectangularLattice
from ..phase_extraction import phase_coverage
from ..sparameters import compute_s_parameters
from ..phase8c.data_source import DEFAULT_CHECKPOINT_JSON, load_full_database

# The same 18 L values (1.0-4.4 mm, 0.2 mm step) Phase 7 swept -- not re-invented (Sec.
# "Physics sanity" of the brief: "e.g. L = W over the Phase 7 range, ~10-20 values").
L_VALUES_MM = library_schema.L_VALUES_MM
SWEEP_FREQUENCIES_HZ = (10.0e9, 11.25e9, 12.5e9, 13.75e9, 15.0e9)
DESIGN_FREQUENCY_HZ = library_schema.LibraryConfig.default().f_design_hz   # 12.5 GHz, as Phase 7


def _library_substrate():
    cfg = library_schema.LibraryConfig.default()
    sub = cfg.substrate
    return Substrate(name=sub["name"], eps_r=sub["eps_r"], tan_delta=sub["tan_delta"], thickness=sub["thickness_m"])


def grounded_library_cell(length_width_m):
    """A square (L=W) rectangular-PEC-patch UnitCell on the exact Phase 7 period/substrate,
    grounded (``backing="pec"``)."""
    period_m = library_schema.LibraryConfig.default().period_m
    lattice = RectangularLattice(period_m, period_m)
    geometry = Geometry("rectangle", length_width_m, length_width_m)
    return UnitCell(lattice, geometry, _library_substrate(), Metal(), backing="pec")


def _mesh():
    m = library_schema.DEFAULT_MESH
    return MeshSettings(m["cells_x"], m["cells_y"], m["floquet_order"])


def _planes():
    rp = library_schema.DEFAULT_REFERENCE_PLANES
    return ReferencePlanes(rp["top_air_thickness"], None, _library_substrate().thickness)


@dataclass(frozen=True)
class SweepPoint:
    length_width_mm: float
    frequency_hz: float
    s11_xx_mag: float
    s11_xx_phase_deg: float
    s11_yx_mag: float           # cross-pol sanity check (should be ~0 for an unrotated square)


def run_grounded_patch_sweep(l_values_mm=L_VALUES_MM, frequencies_hz=SWEEP_FREQUENCIES_HZ):
    """Direct ``compute_s_parameters`` for every (L, frequency) pair; real MoM solves, no fitting,
    no lookup. Returns a tuple of ``SweepPoint``."""
    mesh, planes = _mesh(), _planes()
    points = []
    for l_mm in l_values_mm:
        cell = grounded_library_cell(l_mm * 1e-3)
        for f in frequencies_hz:
            s = compute_s_parameters(cell, f, mesh, planes, warn_grating=False)
            points.append(SweepPoint(
                l_mm, f, s.magnitude("S11_xx"), s.phase_deg("S11_xx"), s.magnitude("S11_yx"),
            ))
    return tuple(points)


def coverage_by_frequency(points):
    """{frequency_hz: phase_extraction.PhaseCoverage} from the sweep's own S11_xx phases --
    ``phase_coverage`` is Phase 7's own, unmodified circular-coverage function (Sec. "Physics
    sanity": "measure reflection-phase coverage, largest gap")."""
    freqs = sorted({p.frequency_hz for p in points})
    out = {}
    for f in freqs:
        table = [{"phase_deg": p.s11_xx_phase_deg} for p in points if p.frequency_hz == f]
        out[f] = phase_coverage(table)
    return out


# ------------------------------------------------------------- comparison with ungrounded Phase 7

def ungrounded_phase7_coverage(checkpoint_path=DEFAULT_CHECKPOINT_JSON, frequency_hz=DESIGN_FREQUENCY_HZ):
    """S11_xx and S21_xx circular phase coverage of the real, already-built, unmodified Phase 7
    library at ``frequency_hz`` -- read directly from the checkpoint (``Database``, Phase 5
    format), never re-solved, never re-derived from a doc. Returns
    ``{"status": "NOT_AVAILABLE_FROM_CURRENT_DATA"}`` if the checkpoint is not present (this
    project never fabricates missing data, CLAUDE.md)."""
    try:
        db = load_full_database(checkpoint_path)
    except FileNotFoundError:
        return {"status": "NOT_AVAILABLE_FROM_CURRENT_DATA"}
    freqs = db.frequencies
    idx = int(np.argmin(np.abs(freqs - frequency_hz)))
    if abs(freqs[idx] - frequency_hz) > 1e-6 * frequency_hz:
        raise ValueError(f"{frequency_hz} Hz is not an exact grid point of the Phase 7 library")
    result = {"status": "AVAILABLE", "frequency_hz": float(freqs[idx])}
    for name in ("S11_xx", "S21_xx"):
        table = [
            {"phase_deg": db.record(c.index, idx).phase_deg(name)}
            for c in db.cells if c.status == "ok"
        ]
        cov = phase_coverage(table)
        result[name] = {
            "coverage_deg": cov.coverage_deg, "largest_gap_deg": cov.largest_gap_deg,
            "n_points": cov.n_points,
        }
    return result


# ---------------------------------------------------------------------- convergence (Sec. "mesh/
# Floquet convergence study for one grounded patch")

@dataclass(frozen=True)
class ConvergencePoint:
    cells_per_side: int
    floquet_order: int
    s11_xx_mag: float
    s11_xx_phase_deg: float


def ungrounded_phase7_diagonal_table(checkpoint_path=DEFAULT_CHECKPOINT_JSON, frequency_hz=DESIGN_FREQUENCY_HZ):
    """[{"L_mm", "phase_deg"}] for the L=W diagonal of the real Phase 7 18x18 sweep at
    ``frequency_hz`` -- the exact same L values as ``run_grounded_patch_sweep``'s default, so the
    grounded-vs-ungrounded plot overlays directly. Returns ``[]`` if the checkpoint is absent
    (never fabricated)."""
    try:
        db = load_full_database(checkpoint_path)
    except FileNotFoundError:
        return []
    freqs = db.frequencies
    idx = int(np.argmin(np.abs(freqs - frequency_hz)))
    if abs(freqs[idx] - frequency_hz) > 1e-6 * frequency_hz:
        raise ValueError(f"{frequency_hz} Hz is not an exact grid point of the Phase 7 library")
    rows = []
    for c in db.cells:
        if c.status != "ok":
            continue
        l_mm, w_mm = c.parameters["length_m"] * 1e3, c.parameters["width_m"] * 1e3
        if abs(l_mm - w_mm) > 1e-9:
            continue   # diagonal (L = W) only, to match the grounded sweep
        rows.append({"L_mm": l_mm, "phase_deg": db.record(c.index, idx).phase_deg("S11_xx")})
    return rows


def convergence_study(length_width_m=2.8e-3, frequency_hz=DESIGN_FREQUENCY_HZ, sizes=(8, 16, 32, 64)):
    """S11_xx of one representative grounded patch (C3's own 2.8x2.8 mm, this project's own
    already-used center geometry -- phase9/phase10/phase11) at increasing mesh/Floquet
    resolution, ``floquet_order = 2 * cells_per_side`` each time (MOM_FORMULATION.md's own
    "2-4x cells" accuracy rule of thumb)."""
    period_m = library_schema.LibraryConfig.default().period_m
    lattice = RectangularLattice(period_m, period_m)
    geometry = Geometry("rectangle", length_width_m, length_width_m)
    cell = UnitCell(lattice, geometry, _library_substrate(), Metal(), backing="pec")
    planes = _planes()
    points = []
    for n in sizes:
        mesh = MeshSettings(n, n, 2 * n)
        s = compute_s_parameters(cell, frequency_hz, mesh, planes, warn_grating=False)
        points.append(ConvergencePoint(n, 2 * n, s.magnitude("S11_xx"), s.phase_deg("S11_xx")))
    return tuple(points)
