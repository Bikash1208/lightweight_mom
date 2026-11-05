"""fcc6691 regression: every PRE-EXISTING rotated-load / diode case is bit-identical to `main` (before the fix), except the two documented behaviour changes.

`tests/data/rotated_load_main_reference.json` = sha256 fingerprints of the exact complex128 results produced by the `main` checkout (generator:
`scripts/research/rotated_load_regression_cases.py`, run in a `main` worktree): load matrices of x / y / 0 / 90 deg diode loads, rotated loads with an explicit
`transverse_conducts` (False and True), `hfss_chain` along x, and the CP-2 cell (45 / 135 deg stubs, explicit transverse_conducts=True) S11 / S21 for the ideal and
MADP kits in three states.  Changed on purpose (default of a rotated load that does NOT pass `transverse_conducts`): `diode_load(..., 30.0)` and `hfss_chain(..., 30.0)`."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REF = json.loads((ROOT / "tests/data/rotated_load_main_reference.json").read_text(encoding="utf-8"))
CHANGED = {"load_matrix_d30_default", "load_matrix_hfss_chain_30_default"}


@pytest.fixture(scope="module")
def now():
    spec = importlib.util.spec_from_file_location("rotated_load_regression_cases", ROOT / "scripts/research/rotated_load_regression_cases.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.compute()


def test_reference_covers_the_cases(now):
    assert set(REF) == set(now) and len(REF) == 21


@pytest.mark.parametrize("key", sorted(set(REF) - CHANGED))
def test_preexisting_case_is_bit_identical_to_main(now, key):
    assert now[key]["sha256"] == REF[key]["sha256"] and now[key]["shape"] == REF[key]["shape"]


def test_documented_behaviour_change_only_for_default_rotated_loads(now):
    for key in CHANGED:
        assert now[key]["sha256"] != REF[key]["sha256"]  # the fix changes exactly these
    # the new default of a rotated load equals the explicit transverse_conducts=True load
    assert now["load_matrix_d30_default"]["sha256"] == now["load_matrix_d30_explicit_conducting"]["sha256"]
