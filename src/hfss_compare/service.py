"""Upload-driven HFSS-vs-MoM comparison, plus a small on-disk ledger of past comparisons (feeds
the HFSS gate records the brief asks for, `docs/phases/PHASE13_WEBAPP.md`'s own gate). Reuses
`phase9.hfss_reference.read_touchstone` (structural validation) and
`phase12d.hfss_ingest.compare_single_cell` (the actual comparison) unmodified.
"""

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ..labels import HFSS_VALIDATED_OFF_RESONANCE, hfss_residual_at_resonance
from ..phase9.hfss_reference import read_touchstone
from ..phase12d.hfss_ingest import MOM_REFERENCE_CSV, compare_single_cell

DEFAULT_LEDGER_PATH = Path(__file__).resolve().parent.parent.parent / "results" / "hfss_compare_ledger.json"


@dataclass(frozen=True)
class ComparisonSummary:
    case_id: str
    touchstone_path: str
    n_ports: int
    frequency_count: int
    mean_dphase_deg: float
    max_abs_dphase_deg: float
    max_dmag_db: float
    meets_tolerance: bool
    labels: tuple


def validate_touchstone_upload(path):
    """Structural validation only (parses the file, doesn't compare it yet) -- raises
    `phase9.hfss_reference.TouchstoneError` with a clear message on anything malformed, before the
    app lets the user pick a case to compare against."""
    return read_touchstone(path)  # raises TouchstoneError on malformed input


def compare_uploaded_file(
    case_id,
    touchstone_path,
    mom_csv_path=MOM_REFERENCE_CSV,
    top_port_name="FloquetPort1",
    bottom_port_name="FloquetPort2",
):
    comparison = compare_single_cell(
        case_id,
        touchstone_path=touchstone_path,
        mom_csv_path=mom_csv_path,
        top_port_name=top_port_name,
        bottom_port_name=bottom_port_name,
    )
    labels = []
    if comparison.meets_tolerance:
        labels.append(HFSS_VALIDATED_OFF_RESONANCE)
    else:
        labels.append(hfss_residual_at_resonance(comparison.max_abs_dphase_deg))
    return ComparisonSummary(
        case_id=case_id,
        touchstone_path=str(touchstone_path),
        n_ports=None,
        frequency_count=len(comparison.frequency_ghz),
        mean_dphase_deg=comparison.mean_dphase_deg,
        max_abs_dphase_deg=comparison.max_abs_dphase_deg,
        max_dmag_db=comparison.max_dmag_db,
        meets_tolerance=comparison.meets_tolerance,
        labels=tuple(labels),
    ), comparison


def load_ledger(path=None):
    # `path` resolves against the module-level DEFAULT_LEDGER_PATH at CALL time, not def time --
    # so a caller (or a test, via monkeypatch) that overrides the module attribute instead of
    # passing `path` explicitly still gets the intended location, not a value frozen at import.
    path = Path(path) if path is not None else DEFAULT_LEDGER_PATH
    if not path.exists():
        return []
    return json.loads(path.read_text())


def record_comparison(summary: ComparisonSummary, path=None):
    """Appends `summary` to the on-disk ledger (creates it if absent)."""
    path = Path(path) if path is not None else DEFAULT_LEDGER_PATH
    ledger = load_ledger(path)
    entry = asdict(summary)
    entry["labels"] = list(entry["labels"])
    ledger.append(entry)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ledger, indent=2, default=str))
    return path
