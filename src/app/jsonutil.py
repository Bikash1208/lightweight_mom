"""Phase 13 Part 3: flattens solver return values (numpy scalars/arrays, complex numbers, frozen
dataclasses, Path objects) into plain JSON-safe Python -- the one place that logic lives, so every
router can just ``return to_jsonable(result)`` instead of repeating ad hoc ``.tolist()``/``float()``
conversions at each call site.
"""

import dataclasses
import math
from pathlib import Path

import numpy as np


def _finite_or_none(x):
    """`inf`/`-inf`/`nan` are valid Python floats but not valid JSON tokens Starlette's default
    encoder happens to accept (`json.dumps(..., allow_nan=True)` emits the bare, non-standard
    `Infinity`/`NaN`) -- browser `JSON.parse` rejects them outright. `None` is the standard-JSON
    way to carry "not a finite number" (e.g. `verify.supercell_check.CostEstimate`'s refused-case
    `estimated_runtime_s=float("inf")`)."""
    return x if math.isfinite(x) else None


def to_jsonable(value):
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return to_jsonable(value.tolist())
    if isinstance(value, np.floating):
        return _finite_or_none(float(value))
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, float):
        return _finite_or_none(value)
    if isinstance(value, complex):
        return {"real": _finite_or_none(value.real), "imag": _finite_or_none(value.imag)}
    if isinstance(value, Path):
        return str(value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return to_jsonable(dataclasses.asdict(value))
    return value
