"""Design-frequency extraction (req. 8), the phase-vs-geometry table (req. 10), collision
grouping (req. 11) and phase-coverage analysis (req. 12) -- exact-grid-point behaviour, no
interpolation, and the coverage/gap arithmetic checked against hand-computed synthetic cases
before trusting it on real solver output."""

import warnings

import pytest

from src.library_builder import build_library_database
from src.library_schema import LibraryConfig
from src.phase_extraction import (
    DesignFrequencyRow,
    design_frequency_index,
    extract_design_frequency,
    phase_collision_groups,
    phase_coverage,
    phase_vs_geometry_table,
    resolve_phase_source,
)

warnings.filterwarnings("ignore", category=UserWarning)

_SUBSTRATE = {"name": "Rogers RT/Duroid 5880", "eps_r": 2.2, "tan_delta": 0.0009, "thickness_m": 1.575e-3}


def small_config():
    return LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=12.0e9,
        f_max_hz=13.0e9,
        f_step_hz=0.5e9,
        f_design_hz=12.5e9,
        l_values_mm=(1.0, 2.4, 4.2),
        w_values_mm=(1.0, 2.4, 4.2),
    )


@pytest.fixture(scope="module")
def built():
    cfg = small_config()
    return cfg, build_library_database(cfg, out_dir=None).database


# ------------------------------------------------------------------------ design frequency index


def test_design_frequency_index_is_exact_grid_point(built):
    cfg, database = built
    idx = design_frequency_index(database, cfg.f_design_hz)
    assert database.frequencies[idx] == pytest.approx(cfg.f_design_hz)


def test_design_frequency_index_rejects_a_non_grid_frequency(built):
    _, database = built
    with pytest.raises(ValueError, match="not one of"):
        design_frequency_index(database, 12.33e9)


# ------------------------------------------------------------------------------- extraction


def test_extract_design_frequency_has_one_row_per_geometry(built):
    cfg, database = built
    rows = extract_design_frequency(database, cfg)
    assert len(rows) == 9
    assert all(isinstance(r, DesignFrequencyRow) for r in rows)
    assert all(r.solver_status == "OK" for r in rows)
    assert all(r.s11 is not None and r.s21 is not None for r in rows)


def test_extract_design_frequency_values_match_the_database_directly(built):
    cfg, database = built
    rows = extract_design_frequency(database, cfg)
    idx = design_frequency_index(database, cfg.f_design_hz)
    row = next(r for r in rows if r.l_mm == 2.4 and r.w_mm == 2.4)
    cell = next(
        c
        for c in database.cells
        if abs(c.parameters["length_m"] - 2.4e-3) < 1e-12 and abs(c.parameters["width_m"] - 2.4e-3) < 1e-12
    )
    record = next(r for r in database.records_of(cell.index) if r.frequency_index == idx)
    assert row.s11["real"] == pytest.approx(record.component("S11_xx").real)
    assert row.s11["imag"] == pytest.approx(record.component("S11_xx").imag)
    assert row.s11["phase_deg"] == pytest.approx(record.phase_deg("S11_xx"))
    assert row.s21["mag"] == pytest.approx(record.magnitude("S21_xx"))


def test_phase_vs_geometry_table_has_exact_required_columns(built):
    cfg, database = built
    rows = extract_design_frequency(database, cfg)
    table = phase_vs_geometry_table(rows)
    assert len(table) == 9
    assert set(table[0]) == {"L_mm", "W_mm", "phase_deg", "magnitude_db", "S11_real", "S11_imag"}


def test_extract_design_frequency_raises_on_incomplete_database(built):
    """A config whose L/W grid the database does not fully cover must fail loudly, not silently
    return a partial table."""
    _, database = built
    mismatched_cfg = LibraryConfig(
        period_m=4.8e-3,
        substrate=_SUBSTRATE,
        f_min_hz=12.0e9,
        f_max_hz=13.0e9,
        f_step_hz=0.5e9,
        f_design_hz=12.5e9,
        l_values_mm=(1.0, 2.4, 3.6),
        w_values_mm=(1.0, 2.4, 3.6),  # 3.6 mm not in the built database
    )
    with pytest.raises(ValueError, match="no cell for"):
        extract_design_frequency(database, mismatched_cfg)


# --------------------------------------------------------------------------------- phase source


def test_resolve_phase_source_is_s11_xx_specification_mandated():
    res = resolve_phase_source()
    assert res.component == "S11_xx"
    assert res.status == "SPECIFICATION_MANDATED"
    assert res.note  # non-empty, explains the ground-plane caveat


# ------------------------------------------------------------------------------- coverage/gaps


def _table(phases_deg):
    return [
        {"L_mm": i, "W_mm": i, "phase_deg": p, "magnitude_db": -10.0, "S11_real": 0.0, "S11_imag": 0.0}
        for i, p in enumerate(phases_deg)
    ]


def test_coverage_four_evenly_spaced_points_has_90deg_gaps_and_270_coverage():
    cov = phase_coverage(_table([-135.0, -45.0, 45.0, 135.0]))
    assert cov.largest_gap_deg == pytest.approx(90.0)
    assert cov.coverage_deg == pytest.approx(270.0)


def test_coverage_single_point_has_full_gap_and_zero_coverage():
    cov = phase_coverage(_table([10.0]))
    assert cov.largest_gap_deg == pytest.approx(360.0)
    assert cov.coverage_deg == pytest.approx(0.0)


def test_coverage_wraparound_gap_detected_across_180_boundary():
    """Two points near +180/-180 should NOT be seen as far apart; the wrap-around gap from the
    largest value back to the smallest (+360) must be computed, not just consecutive diffs."""
    cov = phase_coverage(_table([-179.0, 0.0, 179.0]))
    # gaps: -179->0 = 179, 0->179 = 179, wrap 179->(-179+360)=181 -> 2
    assert cov.largest_gap_deg == pytest.approx(179.0)
    assert cov.coverage_deg == pytest.approx(181.0)


def test_coverage_does_not_claim_360_merely_because_points_exist():
    """A cautious library with only 3 clustered points must not report near-total coverage."""
    cov = phase_coverage(_table([1.0, 2.0, 3.0]))
    assert cov.coverage_deg < 10.0
    assert cov.largest_gap_deg > 350.0


def test_coverage_rejects_empty_table():
    with pytest.raises(ValueError):
        phase_coverage([])


# ------------------------------------------------------------------------------- collisions


def test_collision_groups_finds_clustered_phases():
    table = _table([10.0, 12.0, 100.0, -150.0])  # first two are within a 10-deg bin together
    groups = phase_collision_groups(table, bin_deg=10.0)
    assert len(groups) == 1
    assert groups[0]["n_candidates"] == 2
    assert {c["L_mm"] for c in groups[0]["candidates"]} == {0, 1}


def test_collision_groups_none_when_all_isolated():
    table = _table([0.0, 90.0, 180.0, -90.0])
    groups = phase_collision_groups(table, bin_deg=5.0)
    assert groups == []


def test_collision_groups_rejects_nonpositive_bin_width():
    with pytest.raises(ValueError):
        phase_collision_groups(_table([0.0]), bin_deg=0.0)


def test_collision_groups_does_not_choose_a_winner():
    """Every candidate in a colliding bin is kept, in a stable (sorted) order -- no automatic
    selection."""
    table = _table([5.0, 5.0, 5.0])
    groups = phase_collision_groups(table, bin_deg=10.0)
    assert len(groups) == 1
    assert groups[0]["n_candidates"] == 3
