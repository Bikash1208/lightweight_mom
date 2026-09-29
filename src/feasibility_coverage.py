"""Phase 8A Sec. 8-16: library phase-coverage analysis, read directly from the delivered
`unit_cell_library.json` (never re-solved, never interpolated across frequency) and built on top
of Phase 7's own circular-coverage arithmetic (`phase_extraction.phase_coverage`,
`phase_collision_groups`) so the numbers here are identical to, and regression-testable against,
`LIBRARY_REPORT.md`'s stored 24.84 deg / 335.16 deg.
"""

from dataclasses import dataclass

import numpy as np

from .library_classifier import detect_resonance, phase_usability
from .phase_extraction import phase_coverage


def library_table_at_frequency(library_doc, frequency_ghz, rtol=1e-6):
    """Sec. 8's exact-grid-point extraction, for every successfully-solved cell in
    ``library_doc`` (`unit_cell_library.write_library_json`'s schema) -- an exact frequency match
    only; raises rather than interpolating (Sec. 9, 26)."""
    rows = []
    for cell in library_doc["cells"]:
        if cell.get("solver_status") != "OK":
            continue
        freqs = cell["frequency_GHz"]
        matches = [i for i, f in enumerate(freqs) if abs(f - frequency_ghz) <= rtol * max(abs(frequency_ghz), 1.0)]
        if not matches:
            raise ValueError(
                f"frequency {frequency_ghz} GHz is not an exact grid point for cell "
                f"{cell['cell_id']}; no interpolation is performed"
            )
        e = cell["S11"][matches[0]]
        rows.append({
            "L_mm": cell["L_mm"], "W_mm": cell["W_mm"], "phase_deg": e["phase_deg"],
            "magnitude_db": e["mag_db"], "S11_real": e["real"], "S11_imag": e["imag"],
        })
    if not rows:
        raise ValueError("library_table_at_frequency: no successfully-solved cells in this library")
    return rows


@dataclass(frozen=True)
class PhaseDensity:
    n_valid_geometries: int
    n_distinct_phase_states: int
    phase_span_deg: float            # naive max-min (NOT circular); paired with coverage_deg separately
    largest_gap_deg: float
    smallest_gap_deg: float
    median_gap_deg: float
    phase_density_per_deg: float     # n_valid_geometries / 360: overall packing rate


def phase_density(table, dedup_eps_deg=1e-6):
    """Sec. 12: distribution of available phase states. ``dedup_eps_deg`` is the explicit,
    documented tolerance for treating two phase samples as "the same" state (bit-level dedup by
    default -- no two independently-simulated geometries are expected to land on the exact same
    phase, so this is a genuine-duplicate check, not a clustering heuristic)."""
    phases = sorted(row["phase_deg"] for row in table)
    n = len(phases)
    if n == 0:
        raise ValueError("phase_density: no successful cells")
    distinct = 1
    for i in range(1, n):
        if phases[i] - phases[i - 1] > dedup_eps_deg:
            distinct += 1
    if n < 2:
        gaps = [360.0]
    else:
        gaps = [phases[i + 1] - phases[i] for i in range(n - 1)]
        gaps.append((phases[0] + 360.0) - phases[-1])
    return PhaseDensity(
        n_valid_geometries=n, n_distinct_phase_states=distinct,
        phase_span_deg=phases[-1] - phases[0],
        largest_gap_deg=max(gaps), smallest_gap_deg=min(gaps), median_gap_deg=float(np.median(gaps)),
        phase_density_per_deg=n / 360.0,
    )


def phase_coverage_vs_frequency(library_doc):
    """Sec. 14: circular phase coverage at every one of the library's stored frequency points
    (not just f_design), so a narrow coverage at 12.5 GHz can be checked against the rest of the
    10-15 GHz band rather than assumed representative."""
    ok_cells = [c for c in library_doc["cells"] if c.get("solver_status") == "OK"]
    if not ok_cells:
        raise ValueError("phase_coverage_vs_frequency: no successful cells in this library")
    freq_ghz = ok_cells[0]["frequency_GHz"]
    for c in ok_cells:
        if c["frequency_GHz"] != freq_ghz:
            raise ValueError(
                f"cell {c['cell_id']} has a different frequency grid than the library's common "
                "grid; run feasibility_library_integrity.check_library_integrity first"
            )
    results = []
    for i, f_ghz in enumerate(freq_ghz):
        table = [
            {"L_mm": c["L_mm"], "W_mm": c["W_mm"], "phase_deg": c["S11"][i]["phase_deg"]}
            for c in ok_cells
        ]
        cov = phase_coverage(table)
        results.append({
            "frequency_GHz": f_ghz, "coverage_deg": cov.coverage_deg,
            "largest_gap_deg": cov.largest_gap_deg,
            "phase_min_deg": cov.phase_min_deg, "phase_max_deg": cov.phase_max_deg,
        })
    return results


def resonance_summary(library_doc):
    """Sec. 15: resonance classification recomputed straight from the stored |S11xx(f)| curve of
    every successful cell, via `library_classifier.detect_resonance` (Phase 7, unmodified) --
    same criterion, no new threshold. Returns (n_resonant, n_non_resonant, resonance_frequencies_hz)."""
    ok_cells = [c for c in library_doc["cells"] if c.get("solver_status") == "OK"]
    resonant_freqs = []
    n_resonant = n_non_resonant = 0
    for c in ok_cells:
        freqs_hz = [f * 1e9 for f in c["frequency_GHz"]]
        mag = [e["mag"] for e in c["S11"]]
        result = detect_resonance(freqs_hz, mag)
        if result.is_resonant:
            n_resonant += 1
            resonant_freqs.append(result.resonance_frequency_hz)
        else:
            n_non_resonant += 1
    return n_resonant, n_non_resonant, tuple(sorted(resonant_freqs))


def describe_coverage_vs_resonance(freq_coverage, n_resonant_cells, resonant_freqs_hz=()):
    """Sec. 16: descriptive only -- no causal claim ("resonance = required") is made regardless of
    what is found. If no resonant cells exist in-band (the Phase 7 library's actual, verified
    result), no correlation can be computed at all, and this function says so plainly rather than
    fabricating a statistic from zero data points."""
    if n_resonant_cells == 0:
        return {
            "resonant_cells_in_band": 0,
            "frequency_of_peak_coverage_GHz": max(freq_coverage, key=lambda r: r["coverage_deg"])["frequency_GHz"],
            "note": (
                "0 cells are classified RESONANT within the library's frequency band, so no "
                "resonance-vs-coverage correlation can be computed from this library. Any "
                "relationship between resonance and phase coverage remains untested by this data."
            ),
        }
    peak = max(freq_coverage, key=lambda r: r["coverage_deg"])
    return {
        "resonant_cells_in_band": n_resonant_cells,
        "resonance_frequencies_GHz": tuple(f / 1e9 for f in resonant_freqs_hz),
        "frequency_of_peak_coverage_GHz": peak["frequency_GHz"],
        "note": (
            "Resonant cells exist in-band; this module reports where the widest per-frequency "
            "phase coverage occurs alongside where resonances occur, but does not claim a causal "
            "relationship between the two without further analysis."
        ),
    }


def usable_geometry_count(table):
    """Sec. 12's "number of usable geometries": PHASE_USABLE per `library_classifier.phase_usability`
    (>= `library_schema.PHASE_USABLE_MIN_MAG_DB`), reused rather than re-derived."""
    return sum(1 for row in table if phase_usability(row["magnitude_db"]) == "PHASE_USABLE")
