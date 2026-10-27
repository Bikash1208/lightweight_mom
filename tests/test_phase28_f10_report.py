"""Phase 28 item 2 at 10.0 GHz: the 10.0 GHz report and the frequency comparison are generated from their JSON only, and every number in them is a value of that JSON.
Skipped while the 10.0 GHz results are not in the tree."""

import json
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_compare_frequencies as cmpf  # noqa: E402
import phase28_feasibility_report as rep  # noqa: E402

F10 = ROOT / "results" / "phase28" / "f10"
JSON10 = F10 / "feasibility_results.json"
REPORT10 = F10 / "FEASIBILITY_REPORT.md"
COMP_JSON = F10 / "frequency_comparison.json"
COMP_MD = F10 / "FREQUENCY_COMPARISON.md"
JSON1054 = ROOT / "results" / "phase28" / "feasibility_results.json"

LABELS = [
    r"Phase \d+",
    r"item \d+",
    r"\d-bit",
    r"SMP1345-079LF(_1mA)?",
    r"MADP-000907-14020",
    r"10\.54 GHz",
    r"10\.0 GHz",
    r"\d{4}-\d{2}-\d{2}",
]
TOKEN = re.compile(r"(?<![\w.])-?\d+(?:\.\d+)?(?:e[+-]?\d+)?(?![\w])")

pytestmark = pytest.mark.skipif(not JSON10.exists(), reason="10.0 GHz results not in the tree")


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


def _matches(token, leaves, rel=None):
    t = float(token)
    if "e" in token:
        return any((v == 0) if t == 0 else (v != 0 and abs(v - t) <= 0.06 * abs(v)) for v in leaves)
    if rel is not None:  # %.4g formatting: four significant digits
        return any(abs(v - t) <= rel * max(abs(v), abs(t), 1e-300) for v in leaves if math.isfinite(v))
    nd = len(token.split(".")[1]) if "." in token else 0
    return any(abs(round(v, nd) - t) < 1e-9 for v in leaves if math.isfinite(v))


def _strip(text):
    for pat in LABELS:
        text = re.sub(pat, " ", text)
    return text


def test_f10_report_is_generated_from_its_json_and_every_number_is_in_it():
    doc = json.loads(JSON10.read_text(encoding="utf-8"))
    assert doc["constants"]["frequency_ghz"] == 10.0
    assert rep.build(doc) == REPORT10.read_text(encoding="utf-8")
    leaves = []
    _leaves(doc, leaves)
    bad = [t for t in TOKEN.findall(_strip(REPORT10.read_text(encoding="utf-8"))) if not _matches(t, leaves)]
    assert not bad, sorted(set(bad))[:20]


def test_f10_json_declares_its_frequency_table_and_the_launch_assumption():
    doc = json.loads(JSON10.read_text(encoding="utf-8"))
    assert "f10" in doc["table"]["glob"].replace("\\", "/") and doc["table"]["files"]
    assert doc["frequency_note"] and "launch fraction" in doc["frequency_note"]
    assert doc["anchor"]["applicable"] is False and doc["labels"][0] == "MODEL_ONLY"


@pytest.mark.skipif(not COMP_JSON.exists(), reason="comparison not generated")
def test_comparison_is_generated_from_the_two_jsons_and_its_numbers_are_in_its_json():
    comp = json.loads(COMP_JSON.read_text(encoding="utf-8"))
    again = cmpf.compare(
        json.loads(JSON1054.read_text(encoding="utf-8")), json.loads(JSON10.read_text(encoding="utf-8"))
    )
    assert comp == again
    assert cmpf.markdown(comp) == COMP_MD.read_text(encoding="utf-8")
    leaves = []
    _leaves(comp, leaves)
    bad = [t for t in TOKEN.findall(_strip(COMP_MD.read_text(encoding="utf-8"))) if not _matches(t, leaves, rel=6e-4)]
    assert not bad, sorted(set(bad))[:20]
