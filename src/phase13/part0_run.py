"""Part 0 orchestration: runs every carry-over task from PHASE12D_ADDENDUM.md /
PHASE13_WEBAPP_PROMPT.md and assembles one JSON-serializable result. Reuses `convergence.py`,
`resonance_restatement.py`, `hfss_resonance_ingest.py` and `order_power_prediction.py` unmodified
-- this module only calls them and serializes what they return.
"""

import json
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path

import numpy as np

from ..phase12c import design_point as dp
from . import convergence as conv
from . import hfss_resonance_ingest as ing
from . import order_power_prediction as opp
from . import resonance_restatement as rr


def _to_jsonable(obj):
    if is_dataclass(obj) and not isinstance(obj, type):
        return {k: _to_jsonable(v) for k, v in asdict(obj).items()}
    if isinstance(obj, (np.ndarray,)):
        return [_to_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    return obj


@dataclass(frozen=True)
class Part0Result:
    convergence_g4: dict
    resonance_restatement: list
    step_size_cancellation_finding: str
    g1_comparison: dict
    g4_sheet_comparison: dict
    g6_comparison: dict
    g4_error_budget: list
    order_power_verdicts: list


def run_part0():
    """Runs the full Part 0 pipeline. Expensive (several real-MoM Richardson families); intended
    to be run once and its JSON output committed, exactly like every other phase's `run_phaseN_*`
    script."""
    unit_cell_g4 = dp.unit_cell_at(rr.RESONANT_L_MM * 1e-3)
    conv_result = conv.richardson_extrapolate_s11xx(unit_cell_g4, rr.RESONANT_FREQUENCY_HZ, dp.REFERENCE_PLANES)

    restatement = rr.restate_resonance_region_error(reference_fit=conv_result.phase)

    g1 = ing.compare_single_cell("G1")
    g4_sheet = ing.ingest_g4_sheet()
    g6 = ing.ingest_g6()
    g4_budget = ing.error_budget_for_g4()

    order_verdicts = opp.compare_predicted_vs_measured()

    return Part0Result(
        convergence_g4=_to_jsonable(conv_result),
        resonance_restatement=_to_jsonable(restatement),
        step_size_cancellation_finding=rr.STEP_SIZE_CANCELLATION_FINDING,
        g1_comparison=_to_jsonable(g1),
        g4_sheet_comparison=_to_jsonable(g4_sheet),
        g6_comparison=_to_jsonable(g6),
        g4_error_budget=_to_jsonable(g4_budget),
        order_power_verdicts=_to_jsonable(order_verdicts),
    )


def write_part0_json(result, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "part0_carry_overs.json"
    path.write_text(json.dumps(_to_jsonable(result), indent=2, sort_keys=True))
    return path
