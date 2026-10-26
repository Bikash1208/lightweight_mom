"""Phase 13 Part 3: `to_jsonable` must flatten every value type solver code actually returns."""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.app.jsonutil import to_jsonable


def test_numpy_scalars_and_arrays_become_plain_python():
    out = to_jsonable({"a": np.float64(1.5), "b": np.int64(3), "c": np.array([1.0, 2.0])})
    assert out == {"a": 1.5, "b": 3, "c": [1.0, 2.0]}
    json.dumps(out)  # must not raise


def test_complex_becomes_real_imag_dict():
    assert to_jsonable(1 + 2j) == {"real": 1.0, "imag": 2.0}


def test_path_becomes_string():
    p = Path(
        "/tmp/foo.json"
    )  # str(Path) uses the platform separator (backslashes on Windows); the contract is str(path)
    assert to_jsonable(p) == str(p)


def test_nested_dataclass_is_flattened():
    @dataclass(frozen=True)
    class Inner:
        value: np.ndarray

    @dataclass(frozen=True)
    class Outer:
        inner: Inner
        label: str

    out = to_jsonable(Outer(inner=Inner(value=np.array([1, 2, 3])), label="x"))
    assert out == {"inner": {"value": [1, 2, 3]}, "label": "x"}
    json.dumps(out)


def test_plain_values_pass_through_unchanged():
    assert to_jsonable("hello") == "hello"
    assert to_jsonable(True) is True
    assert to_jsonable(None) is None
