"""Phase 7: deterministic, documented cell classification (req. 13-15).

Every field below is an *independent* flag, not one mutually-exclusive label: a cell can be both
RESONANT and PHASE_USABLE, or NON_RESONANT and GRATING_RISK, at once, and forcing a single status
enum would either lose that information or impose an arbitrary priority order the brief never
specifies. Every threshold used is a named constant in `library_schema.py`, documented there with
its rationale, and exercised at its boundary in `tests/test_library_classifier.py` -- never an
unexplained number invented in this module. Grating-order risk reuses
`discretization.is_specular_only`/`grating_onset_frequency` (Phase 4/5, unmodified) rather than a
new frequency/period computation.
"""

from dataclasses import dataclass

import numpy as np

from .discretization import grating_onset_frequency, is_specular_only
from .library_schema import PHASE_USABLE_MIN_MAG_DB, RESONANCE_MIN_PROMINENCE


@dataclass(frozen=True)
class ResonanceResult:
    is_resonant: bool
    resonance_frequency_hz: float | None
    resonance_metric: float | None  # peak prominence (see _prominence), the actual criterion
    peak_magnitude: float | None
    criterion: str  # human-readable, printed in reports/tests


def _prominence(mag, i):
    """Peak prominence of a local maximum at interior index ``i``: the peak value minus the
    higher of the two flanking minima to its left and right within the array (the standard
    topographic-prominence idea, computed directly rather than via an optional peak-finding
    dependency this lightweight project does not otherwise need)."""
    left_min = float(np.min(mag[: i + 1]))
    right_min = float(np.min(mag[i:]))
    return float(mag[i]) - max(left_min, right_min)


def detect_resonance(frequencies_hz, magnitude, min_prominence=RESONANCE_MIN_PROMINENCE):
    """Local maxima of ``magnitude`` (e.g. |S11_xx(f)|) at strictly interior frequency samples
    only -- a maximum sitting at the very edge of the modelled band may belong to a resonance
    whose true peak lies outside it, so it is not reported as a resonance found *inside* the
    band (req. 14). Returns the single largest-prominence resonance, or none found."""
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    magnitude = np.asarray(magnitude, dtype=float)
    n = magnitude.size
    criterion = (
        f"local maximum of |S11_xx(f)| at an interior frequency sample (excluding the two band "
        f"edges) with prominence (peak minus the higher flanking minimum) >= {min_prominence:g}"
    )
    if n < 3:
        return ResonanceResult(False, None, None, None, criterion)

    candidates = []
    for i in range(1, n - 1):
        if magnitude[i] > magnitude[i - 1] and magnitude[i] > magnitude[i + 1]:
            candidates.append((i, _prominence(magnitude, i)))
    if not candidates:
        return ResonanceResult(False, None, None, None, criterion)

    i_best, prominence = max(candidates, key=lambda c: c[1])
    is_resonant = bool(prominence >= min_prominence)
    return ResonanceResult(
        is_resonant=is_resonant,
        resonance_frequency_hz=float(frequencies_hz[i_best]),
        resonance_metric=float(prominence),
        peak_magnitude=float(magnitude[i_best]),
        criterion=criterion,
    )


def grating_risk(unit_cell, frequencies_hz):
    """GRATING_RISK iff a non-specular Floquet order propagates in air at any frequency in
    ``frequencies_hz`` (`discretization.is_specular_only`, unmodified) -- frequency- and
    period-aware, not a blanket label (req. 15). Since only the lattice period (fixed across the
    whole Phase 7 library) determines the onset, every cell shares one `onset_hz`, computed here
    from the actual `unit_cell`, not assumed."""
    onset_hz = grating_onset_frequency(unit_cell)
    at_risk = any(not is_specular_only(unit_cell, f) for f in frequencies_hz)
    return {"is_grating_risk": bool(at_risk), "grating_onset_hz": float(onset_hz)}


def phase_usability(mag_db):
    """PHASE_USABLE / PHASE_UNUSABLE from `library_schema.PHASE_USABLE_MIN_MAG_DB`: below that
    magnitude, the phase of an almost-zero complex number is numerically ill-conditioned, not a
    physically meaningful design quantity. `None` (no data -- a failed cell) is UNCLASSIFIED."""
    if mag_db is None:
        return "UNCLASSIFIED"
    return "PHASE_USABLE" if mag_db >= PHASE_USABLE_MIN_MAG_DB else "PHASE_UNUSABLE"


@dataclass(frozen=True)
class CellClassification:
    cell_id: str
    l_mm: float
    w_mm: float
    geometry_status: str  # VALID | INVALID
    solver_status: str  # OK | FAILED
    resonance: ResonanceResult | None
    phase_status: str  # PHASE_USABLE | PHASE_UNUSABLE | UNCLASSIFIED
    grating: dict | None


def classify_cell(cell_entry, database, design_row, min_prominence=RESONANCE_MIN_PROMINENCE):
    """One `CellClassification` combining every independent flag above for one geometry.
    ``design_row``: this cell's `phase_extraction.DesignFrequencyRow` (for phase usability)."""
    if cell_entry.status != "ok":
        return CellClassification(
            cell_entry.cell_id,
            design_row.l_mm,
            design_row.w_mm,
            "INVALID",
            "FAILED",
            None,
            "UNCLASSIFIED",
            None,
        )

    records = sorted(database.records_of(cell_entry.index), key=lambda r: r.frequency_index)
    expected_n = database.frequencies.size
    solver_status = "OK" if len(records) == expected_n else "FAILED"

    resonance = None
    grating = None
    if solver_status == "OK":
        freqs = np.array([r.frequency for r in records])
        mag = np.array([r.magnitude("S11_xx") for r in records])
        resonance = detect_resonance(freqs, mag, min_prominence)
        grating = grating_risk(database.unit_cell(cell_entry.index), freqs)

    phase_status = phase_usability(design_row.s11["mag_db"] if design_row.s11 is not None else None)

    return CellClassification(
        cell_entry.cell_id,
        design_row.l_mm,
        design_row.w_mm,
        "VALID",
        solver_status,
        resonance,
        phase_status,
        grating,
    )
