"""Task C: HFSS ingestion (Phase 12C Task E carried over -- the brief's own numbering calls it
"Task F" for that reason). Real HFSS data now exist under ``hfss_exports/phase12_grounded_8mm/``
(HFSS 2025 R1; see ``PROVENANCE`` below, recorded exactly as the brief states it, never invented).
Uses ``phase9.hfss_reference.read_touchstone`` (a new, general capability added there this phase)
-- this module is orchestration/comparison only.

Tolerance (defined and justified BEFORE looking at the comparison numbers, per the brief's own
instruction): the brief's own already-measured G1 result is a mean phase difference of -0.29 deg
and |S11| agreement to ~1e-5 -- two independent, real full-wave solvers (this project's MoM and
HFSS's FEM) agreeing this closely is itself the headline result. A generous, round, design-relevant
tolerance an order of magnitude above that measured agreement -- 1.0 deg phase, 0.05 dB magnitude --
is used to decide ``hfss_comparison_meets_tolerance``; it is not tuned to the data (it would also be
failed by, e.g., the brief's own cited "earlier 2-mode-port/copper-ground run" that was off by up to
-4.7 deg).
"""

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..phase import wrapped_difference_deg
from ..phase9.hfss_reference import read_touchstone

EXPORT_ROOT = Path(__file__).resolve().parent.parent.parent / "hfss_exports" / "phase12_grounded_8mm"
TOUCHSTONE_PATH = EXPORT_ROOT / "hfss_results" / "unit_HFSSDesign1.s36p"
EARLIER_RUN_TOUCHSTONE_PATH = EXPORT_ROOT / "hfss_results" / "unit_HFSSDesign1.s4p"
MOM_REFERENCE_CSV = EXPORT_ROOT / "single_cells" / "mom_reference_single_cells.csv"
MANIFEST_PATH = EXPORT_ROOT / "manifest.json"

# The one case with real HFSS Floquet-port data on disk (brief Sec. "Why" point 3): "More cases
# will be added as G<n>.s36p" -- G4 (the resonance case, L=6.0mm) is explicitly still pending.
CASE_WITH_HFSS_DATA = "G1"
ALL_SINGLE_CELL_IDS = ("G1", "G2", "G3", "G4", "G5", "G6")
PENDING_RESONANCE_CASE = "G4"

PHASE_TOLERANCE_DEG = 1.0
MAGNITUDE_TOLERANCE_DB = 0.05

PROVENANCE = {
    "hfss_version": "HFSS 2025 R1",
    "port": "Floquet port on the air-box top at z = +5 mm, no de-embed (= MoM reference plane)",
    "modes": "18 modes per Floquet port (up to 9.25 dB/mm attenuation at 15 GHz)",
    "backing": "PEC ground",
    "patch": "PEC, 35 micron thick (CAD only; solver treats it as zero-thickness PEC, README.md)",
    "convergence": "max delta-S 0.001",
    "sweep": "10-15 GHz, 401 points",
    "time_convention": "exp(+j*omega*t), same as this project's solver -- no conjugation",
    "setup_lesson": (
        "an earlier 2-mode-port / copper-ground run (unit_HFSSDesign1.s4p, kept on disk for this "
        "record) was off by up to -4.7 deg -- the 18-mode PEC-ground setup above is the one that "
        "reproduces the MoM reference this closely"
    ),
    "remaining_trend_hypothesis": (
        "a residual ~0.5 deg phase trend across the band (Sec. below) is plausibly caused by the "
        "35 micron patch thickness (the MoM reference is zero-thickness PEC) -- a hypothesis, not "
        "a finding; not tested here"
    ),
}


@dataclass(frozen=True)
class SingleCellHfssComparison:
    case_id: str
    l_mm: float
    frequency_ghz: np.ndarray
    mom_mag: np.ndarray
    mom_phase_deg: np.ndarray
    hfss_mode1_mag: np.ndarray
    hfss_mode1_phase_deg: np.ndarray
    hfss_mode2_mag: np.ndarray
    hfss_mode2_phase_deg: np.ndarray
    dphase_mode1_deg: np.ndarray
    dphase_mode2_deg: np.ndarray
    dmag_mode1_db: np.ndarray
    hfss_crosspol_mag: np.ndarray
    hfss_top_to_bottom_mag: np.ndarray

    @property
    def mean_dphase_deg(self):
        return float(np.mean(self.dphase_mode1_deg))

    @property
    def max_abs_dphase_deg(self):
        return float(np.max(np.abs(self.dphase_mode1_deg)))

    @property
    def max_dmag_db(self):
        return float(np.max(np.abs(self.dmag_mode1_db)))

    @property
    def max_crosspol_mag(self):
        return float(np.max(self.hfss_crosspol_mag))

    @property
    def max_top_to_bottom_mag(self):
        return float(np.max(self.hfss_top_to_bottom_mag))

    @property
    def meets_tolerance(self):
        return self.max_abs_dphase_deg <= PHASE_TOLERANCE_DEG and self.max_dmag_db <= MAGNITUDE_TOLERANCE_DB


def _load_mom_reference_rows(path=MOM_REFERENCE_CSV):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def compare_single_cell(
    case_id,
    touchstone_path=TOUCHSTONE_PATH,
    mom_csv_path=MOM_REFERENCE_CSV,
    top_port_name="FloquetPort1",
    bottom_port_name="FloquetPort2",
):
    """Real MoM-vs-HFSS comparison for one single-cell case: reads the Touchstone file, matches it
    against the MoM reference CSV's rows for ``case_id`` at the SAME (exact) frequencies, and
    reports mode1/mode2 self-reflection, cross-pol and top-to-bottom transmission -- reproducing
    (and independently verifying, ``tests/test_phase12d_hfss_ingest.py``) the already-measured
    ``G1_comparison.csv`` numbers the brief cites."""
    net = read_touchstone(touchstone_path)
    groups = net.port_mode_groups()
    top = groups[top_port_name]
    mode1_idx, mode2_idx = top[0] + 1, top[1] + 1  # 1-indexed, Touchstone convention
    bottom = groups.get(bottom_port_name)

    rows = [r for r in _load_mom_reference_rows(mom_csv_path) if r["case"] == case_id]
    if not rows:
        raise ValueError(f"{mom_csv_path}: no MoM reference rows for case {case_id!r}")
    l_mm = float(rows[0]["L_mm"])

    freq_ghz, mom_mag, mom_phase = [], [], []
    h1m, h1p, h2m, h2p, xpol, t2b = [], [], [], [], [], []
    for row in rows:
        f_hz = float(row["f_GHz"]) * 1e9
        idx = int(round((f_hz - net.frequency_hz[0]) / (net.frequency_hz[1] - net.frequency_hz[0])))
        if not (0 <= idx < len(net.frequency_hz)) or not np.isclose(net.frequency_hz[idx], f_hz, rtol=1e-9):
            raise ValueError(f"{touchstone_path}: no exact frequency match for {row['f_GHz']} GHz")
        freq_ghz.append(float(row["f_GHz"]))
        mom_mag.append(float(row["S11_xx_ref5mm_mag"]))
        mom_phase.append(float(row["S11_xx_ref5mm_phase_deg"]))
        s1 = net.s[idx, mode1_idx - 1, mode1_idx - 1]
        s2 = net.s[idx, mode2_idx - 1, mode2_idx - 1]
        h1m.append(abs(s1))
        h1p.append(float(np.degrees(np.angle(s1))))
        h2m.append(abs(s2))
        h2p.append(float(np.degrees(np.angle(s2))))
        xpol.append(abs(net.s[idx, mode2_idx - 1, mode1_idx - 1]))
        if bottom:
            t2b.append(abs(net.s[idx, bottom[0], mode1_idx - 1]))
        else:
            t2b.append(float("nan"))

    freq_ghz, mom_mag, mom_phase = map(np.asarray, (freq_ghz, mom_mag, mom_phase))
    h1m, h1p, h2m, h2p = map(np.asarray, (h1m, h1p, h2m, h2p))
    dphase1 = wrapped_difference_deg(h1p, mom_phase)
    dphase2 = wrapped_difference_deg(h2p, mom_phase)
    dmag_db = 20.0 * np.log10(h1m / mom_mag)
    return SingleCellHfssComparison(
        case_id,
        l_mm,
        freq_ghz,
        mom_mag,
        mom_phase,
        h1m,
        h1p,
        h2m,
        h2p,
        dphase1,
        dphase2,
        dmag_db,
        np.asarray(xpol),
        np.asarray(t2b),
    )


def reproduces_recorded_g1_comparison(comparison, recorded_csv_path=None, phase_tol_deg=1e-6, mag_tol=1e-9):
    """Sec. above: the brief requires this module's own, independently-built comparison to
    reproduce the already-measured G1 result it cites, not merely assert a new one. Compares
    against the pre-existing ``hfss_results/G1_comparison.csv`` this project's own maintainer
    supplied (real HFSS post-processing, not generated by this phase)."""
    path = recorded_csv_path or (TOUCHSTONE_PATH.parent / "G1_comparison.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != len(comparison.frequency_ghz):
        return False
    for row, dphase, dmag_db, xpol in zip(
        rows,
        comparison.dphase_mode1_deg,
        comparison.dmag_mode1_db,
        comparison.hfss_crosspol_mag,
    ):
        if abs(dphase - float(row["dphase_mode1_deg"])) > phase_tol_deg:
            return False
        recorded_dmag_db = 20.0 * np.log10(1.0 + float(row["dmag_mode1"]) / float(row["mom_mag"]))
        if abs(dmag_db - recorded_dmag_db) > 1e-3:  # csv stores linear dmag, not dB; coarser check
            return False
        if abs(xpol - float(row["hfss_crosspol_mag"])) > mag_tol:
            return False
    return True


def load_manifest(path=MANIFEST_PATH):
    return json.loads(Path(path).read_text())


# ---------------------------------------------------------------- geometric identity check

# Phase 12C's own earlier export (``results/phase12c_hfss_export/``) and the one actually sent to
# and returned from HFSS (``hfss_exports/phase12_grounded_8mm/``) overlap at exactly two L values
# (4.8 and 5.6 mm) -- both generated by this project's own ``cad_export.py``, but at different
# times, so this checks they are geometrically the SAME, not just assumed so (brief: "Check they're
# geometrically identical and say so; use the latter's IDs").
PHASE12C_EXPORT_ROOT = Path(__file__).resolve().parent.parent.parent / "results" / "phase12c_hfss_export"
GEOMETRIC_IDENTITY_PAIRS = (
    ("G2", PHASE12C_EXPORT_ROOT / "single_cell_L4.8mm.step", EXPORT_ROOT / "single_cells" / "G2_L4p8mm.step"),
    ("G3", PHASE12C_EXPORT_ROOT / "single_cell_L5.6mm.step", EXPORT_ROOT / "single_cells" / "G3_L5p6mm.step"),
)


@dataclass(frozen=True)
class GeometricIdentityCheck:
    case_id: str
    phase12c_path: str
    hfss_export_path: str
    metal_volume_relative_diff: float
    metal_bbox_matches: bool
    identical: bool


def check_geometric_identity(pairs=GEOMETRIC_IDENTITY_PAIRS, bbox_tol_mm=1e-6, volume_rtol=1e-9):
    """Re-imports both STEP files (independent of anything held in memory,
    ``cad_export.import_and_summarize``, unmodified) and compares their METAL solid's volume and
    bounding box. Requires the optional ``cadquery`` dependency; raises ``ImportError`` if absent
    (the caller decides whether to skip, matching this project's other cadquery-gated tests)."""
    from .. import cad_export

    results = []
    for case_id, path_a, path_b in pairs:
        a = cad_export.import_and_summarize(path_a)
        b = cad_export.import_and_summarize(path_b)
        vol_a, vol_b = a.metal.volume_mm3, b.metal.volume_mm3
        rel_diff = abs(vol_a - vol_b) / vol_a if vol_a else float("inf")
        bbox_matches = all(abs(x - y) <= bbox_tol_mm for x, y in zip(a.metal.bbox_mm, b.metal.bbox_mm))
        identical = bbox_matches and rel_diff <= volume_rtol
        results.append(GeometricIdentityCheck(case_id, str(path_a), str(path_b), rel_diff, bbox_matches, identical))
    return tuple(results)
