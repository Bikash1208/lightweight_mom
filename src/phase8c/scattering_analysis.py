"""Phase 8C Sec. 17-21, 41: full-library complex-response analysis for every candidate control
quantity, at the design frequency and across the band. Builds on `phase_extraction.phase_coverage`
(Phase 7, unmodified) for circular statistics, but -- unlike Phase 7/8A/8B, which only ever
extracted S11xx -- extracts *any* linear component (S11xx, S21xx, S11yy, S21xy, ...) via
`database.component_value`, and tracks the full complex value, not phase alone (Sec. 19-20's
explicit "same phase, different magnitude are not the same state" requirement).
"""

from dataclasses import dataclass

import numpy as np

from ..phase_extraction import phase_coverage
from ..validation_export import _db

DEDUP_EPS = 1e-9  # bit-level dedup of complex values (absolute, since |S|<=1)
DEDUP_EPS_DEG = 1e-6  # same tolerance Phase 8B's distinct_states.py uses for phase-only dedup


@dataclass(frozen=True)
class ScatteringRow:
    l_mm: float
    w_mm: float
    real: float
    imag: float
    mag: float
    mag_db: float
    phase_deg: float


def extract_component_table(database, component_name, frequency_index):
    """One `ScatteringRow` per successful cell, for ``component_name`` (e.g. 'S11_xx', 'S21_yx')
    at ``frequency_index`` -- an exact stored frequency, never interpolated."""
    rows = []
    for cell in database.cells:
        if cell.status != "ok":
            continue
        record = database.record(cell.index, frequency_index)
        z = record.component(component_name)
        mag = abs(z)
        rows.append(
            ScatteringRow(
                l_mm=round(cell.parameters["length_m"] * 1e3, 9),
                w_mm=round(cell.parameters["width_m"] * 1e3, 9),
                real=float(z.real),
                imag=float(z.imag),
                mag=float(mag),
                mag_db=float(_db(np.array([mag]))[0]),
                phase_deg=float(record.phase_deg(component_name)),
            )
        )
    if not rows:
        raise ValueError(f"extract_component_table: no successful cells for {component_name!r}")
    return rows


@dataclass(frozen=True)
class ComplexStateAnalysis:
    n_geometries: int
    n_distinct_complex_states: int
    n_distinct_phase_states: int


def distinct_complex_and_phase_states(rows, dedup_eps=DEDUP_EPS, dedup_eps_deg=DEDUP_EPS_DEG):
    """Sec. 19, 41: counts distinct *complex* states (real+imag together) and distinct *phase*
    states separately -- two cells sharing a phase but not a magnitude must not collapse into one
    complex state (checked directly, not assumed). Grouped by an exact (rounded-to-tolerance) key
    rather than a 1-D sorted-neighbour pass, since a complex value has two independent components
    and a lexicographic sort can place two truly-close points far apart in sort order."""
    n = len(rows)
    n_decimals = max(0, -int(np.floor(np.log10(dedup_eps))))
    complex_keys = {(round(r.real, n_decimals), round(r.imag, n_decimals)) for r in rows}
    n_decimals_deg = max(0, -int(np.floor(np.log10(dedup_eps_deg))))
    phase_keys = {round(r.phase_deg, n_decimals_deg) for r in rows}
    return ComplexStateAnalysis(
        n_geometries=n,
        n_distinct_complex_states=len(complex_keys),
        n_distinct_phase_states=len(phase_keys),
    )


@dataclass(frozen=True)
class CandidateCoverage:
    component: str
    n_points: int
    coverage_deg: float
    largest_gap_deg: float
    magnitude_min: float
    magnitude_max: float
    magnitude_median: float
    magnitude_std: float
    phase_magnitude_correlation: float | None  # Pearson r of (phase_deg, mag) -- descriptive only


def candidate_coverage(rows, component_name):
    table = [{"phase_deg": r.phase_deg} for r in rows]
    cov = phase_coverage(table)
    mags = np.array([r.mag for r in rows])
    phases = np.array([r.phase_deg for r in rows])
    corr = None
    if len(rows) >= 2 and np.std(mags) > 0 and np.std(phases) > 0:
        corr = float(np.corrcoef(phases, mags)[0, 1])
    return CandidateCoverage(
        component=component_name,
        n_points=cov.n_points,
        coverage_deg=cov.coverage_deg,
        largest_gap_deg=cov.largest_gap_deg,
        magnitude_min=float(mags.min()),
        magnitude_max=float(mags.max()),
        magnitude_median=float(np.median(mags)),
        magnitude_std=float(np.std(mags)),
        phase_magnitude_correlation=corr,
    )


def coverage_vs_frequency(database, component_name):
    """Sec. 34: circular coverage of ``component_name`` across every stored frequency point."""
    freqs = database.frequencies
    results = []
    for i, f in enumerate(freqs):
        rows = extract_component_table(database, component_name, i)
        cov = phase_coverage([{"phase_deg": r.phase_deg} for r in rows])
        mags = [r.mag for r in rows]
        results.append(
            {
                "frequency_GHz": f / 1e9,
                "coverage_deg": cov.coverage_deg,
                "largest_gap_deg": cov.largest_gap_deg,
                "magnitude_min": min(mags),
                "magnitude_max": max(mags),
            }
        )
    return results
