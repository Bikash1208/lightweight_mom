"""Phase 7: design-frequency extraction, the phase-vs-geometry lookup table, phase-collision
grouping and phase-coverage analysis -- all read from an already-solved `database.Database`
(Phase 5, unmodified), never recomputed or approximated.

## Phase source (req. 9)

The Phase 7 brief's own lookup-table schema (Sec. 10: `L_mm, W_mm, phase_deg, magnitude_db,
S11_real, S11_imag`) already names its response: **S11**, the co-polarized reflection
coefficient (`S11_xx`, x incident and received, `S_PARAMETERS.md` Sec. 1). This module honours
that literally. It is *not*, however, an independently physically re-derived choice, and this
module says so rather than pretending otherwise: the existing solver stack has **no ground
plane** (README "Assumptions" item 2, unchanged since Phase 1) -- unlike a typical grounded
reflectarray/holographic-antenna unit cell, where S21 is negligible and S11 alone carries the
full response, this project's rectangular patch on a bare substrate slab is a partially
reflective, partially transmissive element (Phase 6's case01 measured |S11_xx| around -13 to
-16 dB and |S21_xx| around -0.1 to -0.2 dB -- mostly transmissive, not mostly reflective). Which
of S11/S21 phase the eventual holographic generator needs depends on an architecture decision
(reflective vs. transmissive hologram, and whether a ground plane will be added later -- out of
scope for this phase and not invented here) that this project has not made. `resolve_phase_source`
therefore returns `component="S11_xx"` together with `status="SPECIFICATION_MANDATED"` (not
`"PHYSICALLY_VALIDATED"`) and an explicit note; `LIBRARY_REPORT.md` restates this as an open item.
`status="UNRESOLVED"` (req. 9's literal fallback) is reserved for a future caller that does not
pin phase_source via the table schema the way req. 10 does here.
"""

from dataclasses import dataclass

import numpy as np

from .library_schema import PHASE_COLLISION_BIN_DEG
from .validation_export import _db

PHASE_SOURCE_COMPONENT = "S11_xx"


@dataclass(frozen=True)
class PhaseSourceResolution:
    component: str
    status: str  # "SPECIFICATION_MANDATED" | "PHYSICALLY_VALIDATED" | "UNRESOLVED"
    note: str


def resolve_phase_source():
    """See module docstring. A function (not a constant) so the reasoning is inspectable and
    testable, and so a future caller can supply model information this module does not have."""
    return PhaseSourceResolution(
        component=PHASE_SOURCE_COMPONENT,
        status="SPECIFICATION_MANDATED",
        note=(
            "Req. 10's lookup-table schema names S11 explicitly; this project's existing solver "
            "has no ground plane, so S21 is not negligible here (unlike a typical grounded "
            "reflectarray element) and whether S11 or S21 phase is the one the eventual "
            "holographic generator needs is an architecture decision this project has not made. "
            "Both S11 and S21 (all four blocks) are stored in full in the library database "
            "regardless of this choice; only the phase-vs-geometry table and the coverage/"
            "collision/resonance analyses below are pinned to S11_xx."
        ),
    )


def design_frequency_index(database, f_design_hz, rtol=1e-9):
    """Exact index of ``f_design_hz`` in ``database``'s frequency grid (req. 8: an exact grid
    point, never interpolated). Raises ValueError if no stored frequency matches within
    floating-point tolerance -- this module never picks the nearest frequency silently."""
    freqs = database.frequencies
    matches = np.nonzero(np.abs(freqs - f_design_hz) <= rtol * max(abs(f_design_hz), 1.0))[0]
    if matches.size == 0:
        raise ValueError(
            f"design frequency {f_design_hz!r} Hz is not one of the {freqs.size} stored "
            f"frequencies ({freqs[0]!r}..{freqs[-1]!r} Hz); no interpolation is performed"
        )
    return int(matches[0])


def _record_metrics(record, component):
    z = record.component(component)
    mag = abs(z)
    return {
        "real": float(z.real),
        "imag": float(z.imag),
        "mag": float(mag),
        "mag_db": float(_db(np.array([mag]))[0]),
        "phase_deg": float(record.phase_deg(component)),
    }


@dataclass(frozen=True)
class DesignFrequencyRow:
    cell_id: str
    l_mm: float
    w_mm: float
    geometry_status: str
    solver_status: str
    s11: dict | None
    s21: dict | None
    s12: dict | None
    s22: dict | None


def extract_design_frequency(database, config):
    """One `DesignFrequencyRow` per geometry in ``config`` (req. 8): S11/S21/S12/S22 (real,
    imag, mag, mag_db, phase_deg) at the exact design-frequency grid point, for every one of the
    18x18 = 324 (L, W) combinations -- including geometries the solver marked invalid or that
    have no record at this frequency, whose S-blocks are `None` (never a substituted/approximate
    value, req. 28) but whose row is still present so the count stays 324 (req. 4)."""
    idx = design_frequency_index(database, config.f_design_hz)
    by_id = {c.cell_id: c for c in database.cells}
    rows = []
    for l_mm in config.l_values_mm:
        for w_mm in config.w_values_mm:
            # Match by stored geometry parameters (mm -> m), not by re-deriving a cell_id here,
            # so this function stays independent of _canonical.cell_id's hashing details.
            cell = next(
                (
                    c
                    for c in database.cells
                    if abs(c.parameters["length_m"] - l_mm * 1e-3) < 1e-15
                    and abs(c.parameters["width_m"] - w_mm * 1e-3) < 1e-15
                ),
                None,
            )
            if cell is None:
                raise ValueError(f"no cell for L={l_mm} mm, W={w_mm} mm in this database -- incomplete library")
            if cell.status != "ok":
                rows.append(DesignFrequencyRow(cell.cell_id, l_mm, w_mm, "INVALID", "FAILED", None, None, None, None))
                continue
            record = next((r for r in database.records_of(cell.index) if r.frequency_index == idx), None)
            if record is None:
                rows.append(DesignFrequencyRow(cell.cell_id, l_mm, w_mm, "VALID", "FAILED", None, None, None, None))
                continue
            rows.append(
                DesignFrequencyRow(
                    cell.cell_id,
                    l_mm,
                    w_mm,
                    "VALID",
                    "OK",
                    _record_metrics(record, "S11_xx"),
                    _record_metrics(record, "S21_xx"),
                    _record_metrics(record, "S12_xx"),
                    _record_metrics(record, "S22_xx"),
                )
            )
    return rows


# ---------------------------------------------------------------- phase-vs-geometry lookup table


def phase_vs_geometry_table(rows):
    """Req. 10's exact table: L_mm, W_mm, phase_deg, magnitude_db, S11_real, S11_imag -- for
    successful cells only (phase_deg/magnitude_db do not exist for a failed one; its L/W still
    appear in ``rows`` above with S-blocks `None`, which is where its absence is recorded)."""
    return [
        {
            "L_mm": r.l_mm,
            "W_mm": r.w_mm,
            "phase_deg": r.s11["phase_deg"],
            "magnitude_db": r.s11["mag_db"],
            "S11_real": r.s11["real"],
            "S11_imag": r.s11["imag"],
        }
        for r in rows
        if r.s11 is not None
    ]


# --------------------------------------------------------------------------- coverage / gaps


@dataclass(frozen=True)
class PhaseCoverage:
    n_points: int
    phase_min_deg: float  # naive min/max of the (unsorted) phase values themselves
    phase_max_deg: float
    naive_span_deg: float  # phase_max_deg - phase_min_deg (NOT circularly aware)
    largest_gap_deg: float  # largest empty arc on the circle between consecutive points
    coverage_deg: float  # 360 - largest_gap_deg: the circularly-aware "spanned" measure
    largest_gap_before_deg: float
    largest_gap_after_deg: float


def phase_coverage(table):
    """Req. 12, computed, not asserted: sorts the wrapped phase_deg values (each already on
    (-180, 180], `phase.principal_phase_deg`) around the circle and finds the largest empty arc
    between consecutive points (including the wrap-around arc from the largest value back to the
    smallest + 360). `coverage_deg = 360 - largest_gap_deg` is reported *alongside*, not instead
    of, the largest gap itself and the naive (non-circular) min/max/span, so a small library that
    happens to have points near +-180 does not get credited with a "span" it does not have."""
    if not table:
        raise ValueError("phase_coverage: no successful cells to analyze")
    phases = sorted(row["phase_deg"] for row in table)
    n = len(phases)
    if n == 1:
        return PhaseCoverage(1, phases[0], phases[0], 0.0, 360.0, 0.0, phases[0], phases[0])
    gaps = [phases[i + 1] - phases[i] for i in range(n - 1)]
    wrap_gap = (phases[0] + 360.0) - phases[-1]
    all_gaps = gaps + [wrap_gap]
    largest = max(all_gaps)
    gap_index = all_gaps.index(largest)
    before = phases[gap_index]
    after = phases[(gap_index + 1) % n]
    return PhaseCoverage(
        n_points=n,
        phase_min_deg=phases[0],
        phase_max_deg=phases[-1],
        naive_span_deg=phases[-1] - phases[0],
        largest_gap_deg=largest,
        coverage_deg=360.0 - largest,
        largest_gap_before_deg=before,
        largest_gap_after_deg=after,
    )


# ------------------------------------------------------------------------- collision analysis


def phase_collision_groups(table, bin_deg=PHASE_COLLISION_BIN_DEG):
    """Req. 11: bin the design-frequency S11 phase of every successful cell into ``bin_deg``-wide
    bins (default `library_schema.PHASE_COLLISION_BIN_DEG`, documented there) and return every
    bin with more than one candidate, each with its full geometry -- no winner is chosen."""
    if bin_deg <= 0.0:
        raise ValueError(f"bin_deg must be > 0, got {bin_deg!r}")
    bins = {}
    for row in table:
        bin_index = int(np.floor((row["phase_deg"] + 180.0) / bin_deg))
        bins.setdefault(bin_index, []).append(row)
    groups = []
    for bin_index, members in sorted(bins.items()):
        if len(members) < 2:
            continue
        lo = -180.0 + bin_index * bin_deg
        groups.append(
            {
                "phase_bin_deg": [lo, lo + bin_deg],
                "n_candidates": len(members),
                "candidates": sorted(members, key=lambda m: (m["L_mm"], m["W_mm"])),
            }
        )
    return groups
