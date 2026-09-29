"""Phase 8B Sec. 30-32: distinct-phase-state analysis, on synthetic data (exact control of
duplicates) and, when present, the real production library (reproducing 324 -> 169)."""

import pytest

from src.feasibility_schema import DEFAULT_LIBRARY_DIR
from src.phase8b.distinct_states import analyze_distinct_states, write_distinct_state_csv

_SYNTHETIC_TABLE = [
    {"L_mm": 1.0, "W_mm": 1.0, "phase_deg": 10.0},
    {"L_mm": 1.2, "W_mm": 1.0, "phase_deg": 10.0},          # exact duplicate of the row above
    {"L_mm": 1.4, "W_mm": 1.0, "phase_deg": 20.0},
    {"L_mm": 1.6, "W_mm": 1.0, "phase_deg": 20.0 + 1e-9},   # duplicate within DEDUP_EPS_DEG
    {"L_mm": 1.8, "W_mm": 1.0, "phase_deg": 30.0},
]


def test_analyze_distinct_states_counts_duplicates_correctly():
    result = analyze_distinct_states(_SYNTHETIC_TABLE)
    assert result.n_geometries == 5
    assert result.n_distinct_phase_states == 3   # {10, 20, 30}
    assert result.n_duplicate_records == 2


def test_geometric_diversity_and_phase_diversity_are_reported_separately():
    result = analyze_distinct_states(_SYNTHETIC_TABLE)
    assert result.n_geometries != result.n_distinct_phase_states
    assert result.n_geometries == 5   # never conflated with the smaller phase-state count


def test_every_row_gets_an_equivalence_group_and_duplicates_share_one():
    result = analyze_distinct_states(_SYNTHETIC_TABLE)
    groups_by_phase = {}
    for row in result.rows:
        groups_by_phase.setdefault(round(row.phase_deg), set()).add(row.equivalence_group)
    assert all(len(g) == 1 for g in groups_by_phase.values())


def test_write_distinct_state_csv(tmp_path):
    result = analyze_distinct_states(_SYNTHETIC_TABLE)
    path = write_distinct_state_csv(result, tmp_path / "out.csv")
    text = path.read_text()
    assert text.splitlines()[0] == "L_mm,W_mm,phase_deg,equivalence_group"
    assert len(text.splitlines()) == 6   # header + 5 rows


def test_analyze_distinct_states_rejects_empty_table():
    with pytest.raises(ValueError):
        analyze_distinct_states([])


@pytest.mark.skipif(
    not (DEFAULT_LIBRARY_DIR / "database" / "unit_cell_library.json").exists(),
    reason="Phase 7 production library not built",
)
def test_real_library_reproduces_324_geometries_169_distinct_states():
    from src.feasibility_coverage import library_table_at_frequency
    from src.feasibility_library_integrity import load_library_bundle

    bundle = load_library_bundle(DEFAULT_LIBRARY_DIR)
    table = library_table_at_frequency(bundle.library_doc, 12.5)
    result = analyze_distinct_states(table)
    assert result.n_geometries == 324
    assert result.n_distinct_phase_states == 169
