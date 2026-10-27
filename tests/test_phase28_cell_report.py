"""Phase 28 item 3: the cell-candidate report is generated from the characterisation JSON files only; every number in it is a value of those files; the derived
quantities (convergence verdicts, tolerance changes, bandwidth, hardware) follow from the raw rows."""

import json
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_cell_characterise as ch  # noqa: E402
import phase28_cell_report as rep  # noqa: E402

CELLS = ROOT / "results" / "phase28" / "f10" / "cells"
REPORT = ROOT / "results" / "phase28" / "f10" / "CELL_CANDIDATES.md"
SEARCH = ROOT / "results" / "phase28" / "f10" / "search_summary.json"
LS = CELLS / "ls_threshold_J.json"
JSONS = sorted(p for p in CELLS.glob("cell_*.json") if not p.name.endswith(".raw.json"))

pytestmark = pytest.mark.skipif(
    not JSONS or not REPORT.exists(), reason="cell characterisation results not in the tree"
)

LABELS = [
    r"Phase \d+",
    r"item \d+",
    r"\d-bit",
    r"MADP-000907-14020",
    r"SMP1345-079LF",
    r"10\.0 GHz",
    r"\b45 degree",
    r"\b2 x 2",
]
TOKEN = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:e[+-]?\d+)?(?![\w])")


def _leaves(o, out):
    if isinstance(o, bool) or o is None:
        return
    if isinstance(o, (int, float)):
        out.append(float(o))
    elif isinstance(o, dict):
        for v in o.values():
            _leaves(v, out)
    elif isinstance(o, (list, tuple)):
        for v in o:
            _leaves(v, out)


def _matches(token, leaves):
    t = float(token)
    nd = len(token.split(".")[1]) if "." in token else 0
    scaled = [v for v in leaves if math.isfinite(v)]
    cand = scaled + [v * 1e3 for v in scaled] + [v * 100 for v in scaled]  # mm from m, percent from fractions
    return any(abs(round(v, nd) - t) < 1e-9 for v in cand)


def _docs():
    return [json.loads(p.read_text(encoding="utf-8")) for p in JSONS]


def test_report_is_what_the_generator_produces_from_the_jsons():
    search = json.loads(SEARCH.read_text(encoding="utf-8")) if SEARCH.exists() else None
    ls = json.loads(LS.read_text(encoding="utf-8")) if LS.exists() else None
    assert rep.build(_docs(), search, ls) == REPORT.read_text(encoding="utf-8")


def test_every_number_in_the_report_is_a_value_of_the_inputs():
    leaves = []
    for d in _docs():
        _leaves(d, leaves)
    if SEARCH.exists():
        _leaves(json.loads(SEARCH.read_text(encoding="utf-8")), leaves)
    if LS.exists():
        _leaves(json.loads(LS.read_text(encoding="utf-8")), leaves)
    text = REPORT.read_text(encoding="utf-8")
    for pat in LABELS:
        text = re.sub(pat, " ", text)
    bad = [t for t in TOKEN.findall(text) if not _matches(t, leaves)]
    assert not bad, sorted(set(bad))[:25]


def test_derived_block_follows_from_the_raw_rows():
    for d in _docs():
        again = ch.derive({k: v for k, v in d.items() if k != "derived"})
        assert again == d["derived"]
        assert (
            d["derived"]["hardware"]["diodes_total"]
            == d["derived"]["hardware"]["n_cells"] * d["derived"]["hardware"]["diodes_per_cell"]
        )
        # the convergence verdict is not 'converged' unless the numbers say so
        for v in d["derived"]["convergence"].values():
            assert v["converged"] == (v["max_relative_change"] <= v["tolerance"])


def test_tolerance_changes_are_relative_to_the_nominal_state_values():
    for d in _docs():
        m_last = d["derived"]["mesh_main"]
        nom = {
            r["state"]: r
            for r in d["rows"]
            if r["scale"] is None and r["f_ghz"] == 10.0 and r["mesh"] == m_last and r["floquet_factor"] == 1
        }
        for x in d["derived"]["tolerance"]:
            assert x["change_delta_ms"] == pytest.approx(
                (x["delta_ms"] - nom[x["state"]]["delta_ms"]) / abs(nom[x["state"]]["delta_ms"]), rel=1e-12
            )
