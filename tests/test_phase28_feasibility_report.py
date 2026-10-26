"""Phase 28 item 2: the report is generated from results/phase28/feasibility_results.json only; every number printed in it is a value of that file."""

import json
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase28_feasibility_report as rep  # noqa: E402

RESULTS = ROOT / "results" / "phase28" / "feasibility_results.json"
REPORT = ROOT / "results" / "phase28" / "FEASIBILITY_REPORT.md"

# Labels and names whose digits are part of a name, not a reported number (the ONLY text exempt from the check).
LABELS = [
    r"Phase \d+",
    r"item \d+",
    r"\d-bit",
    r"SMP1345-079LF(_1mA)?",
    r"MADP-000907-14020",
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
    if "e" in token:
        return any((v == 0) if t == 0 else (v != 0 and abs(v - t) <= 0.06 * abs(v)) for v in leaves)
    nd = len(token.split(".")[1]) if "." in token else 0
    return any(abs(round(v, nd) - t) < 1e-9 for v in leaves if math.isfinite(v))


def _strip_labels(text):
    for pat in LABELS:
        text = re.sub(pat, " ", text)
    return text


def test_every_number_in_the_report_is_a_value_of_the_json():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    leaves = []
    _leaves(doc, leaves)
    text = _strip_labels(REPORT.read_text(encoding="utf-8"))
    bad = [t for t in TOKEN.findall(text) if not _matches(t, leaves)]
    assert not bad, f"numbers in the report that are not values of the JSON: {sorted(set(bad))[:20]}"


def test_report_on_disk_is_exactly_what_the_generator_produces_from_the_json():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert rep.build(doc) == REPORT.read_text(encoding="utf-8")


def test_the_check_itself_catches_an_invented_number():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    leaves = []
    _leaves(doc, leaves)
    assert _matches(
        f"{doc['anchor']['rows'][0]['phase26']:.3f}", leaves
    )  # a real value, matched at the printed precision
    assert not _matches("123456.789", leaves)
    assert not _matches("0.0001234", leaves)


def test_anchor_rows_in_the_json_are_all_within_tolerance_and_the_report_says_so():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    assert doc["anchor"]["all_within_tolerance"] is True
    assert all(r["within_tolerance"] for r in doc["anchor"]["rows"])
    assert "`NOT_AVAILABLE_FROM_CURRENT_DATA" in REPORT.read_text(
        encoding="utf-8"
    )  # the 10 GHz row is labelled, not invented
    assert doc["labels"][0] == "MODEL_ONLY"


def test_an_invented_number_inserted_into_the_report_is_flagged():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    leaves = []
    _leaves(doc, leaves)
    text = _strip_labels(REPORT.read_text(encoding="utf-8") + "\nGain 31.37 dBi at 50 degrees.\n")
    bad = [t for t in TOKEN.findall(text) if not _matches(t, leaves)]
    assert "31.37" in bad


def test_design_point_table_cells_equal_the_json_cell_by_cell():
    """Not only 'some JSON value': every cell of the design-point table is the value of ITS pitch, variant and column."""
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    dp = doc["settings"]["design_point"]
    rows = {}
    for line in REPORT.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 11 and cells[1].startswith(("cp_", "lp_", "lpx_")):
            rows[(cells[0], cells[1])] = cells
    assert rows
    n = 0
    for c in doc["cases"]:
        if (c["aperture_r_mm"], c["theta0_deg"], c["phi0_deg"]) != (dp["aperture_r_mm"], dp["theta0_deg"], dp["phi0_deg"]):
            continue
        for name, v in c["variants"].items():
            cells = rows[(f"{c['pitch_mm']:.1f}", name)]
            expect = [
                f"{v['gain_headline_dbi']:.2f}",
                f"{v['quantisation_loss_db']:.2f}",
                f"{v['quantisation_loss_fixed_source_db']:.2f}",
                f"{v['quantisation_loss_renormalised_db']:.2f}",
                f"{v['radiated_power_change_db']:.2f}",
                f"{v['sll_db']:.1f}",
                f"{v['pointing_error_deg']:.2f}",
            ]
            assert cells[2:9] == expect, (c["pitch_mm"], name, cells[2:9], expect)
            n += 1
    assert n == len(rows)


def test_hardware_table_counts_equal_cells_times_bits_and_the_json():
    doc = json.loads(RESULTS.read_text(encoding="utf-8"))
    dp = doc["settings"]["design_point"]
    seen = 0
    for line in REPORT.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 16 and cells[1] in ("cp_1bit", "cp_2bit", "cp_3bit"):
            case = next(
                c
                for c in doc["cases"]
                if f"{c['pitch_mm']:.1f}" == cells[0] and (c["aperture_r_mm"], c["theta0_deg"], c["phi0_deg"]) == (dp["aperture_r_mm"], dp["theta0_deg"], dp["phi0_deg"])
            )
            b = case["bias"][cells[1]]
            bits = int(cells[1][3])
            assert int(cells[2]) == case["n_cells"] and int(cells[3]) == case["n_cells"] * bits == b["n_diodes"]
            assert float(cells[13]) == pytest.approx(b["dc_power_expected_w"], abs=0.051)
            seen += 1
    assert seen == 12
