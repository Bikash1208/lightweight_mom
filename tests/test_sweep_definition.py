"""Sweep definition: frequency grid, geometry-parameter axes, validation before any solve."""

import math
from pathlib import Path

import numpy as np
import pytest
from sweep_helpers import FREQ_8_12, config_dict, sweep_config

import src.runner as runner
from src.config import ConfigError, config_from_dict, load_config
from src.constants import C0
from src.sweep import (
    MAX_CELLS,
    PARAMETER_NAMES,
    FrequencyGrid,
    ParameterAxis,
    RefinementPolicy,
    SweepDefinition,
    SweepError,
    ValidationPolicy,
    base_parameters,
    inclusive_range,
    plan_sweep,
    unit_cell_from_context,
)

ROOT = Path(__file__).resolve().parent.parent


# ----------------------------------------------------------------- inclusive range / grid


def test_range_includes_both_end_points_with_exact_spacing():
    f = inclusive_range(8e9, 12e9, 0.5e9)
    assert f.size == 9 and f[0] == 8e9 and f[-1] == 12e9
    assert np.allclose(np.diff(f), 0.5e9, rtol=0, atol=1e-3)


@pytest.mark.parametrize(
    "start, stop, step, n",
    [(8e9, 12e9, 0.1e9, 41), (0.1, 1.0, 0.1, 10), (1e-3, 8e-3, 0.5e-3, 15), (0.0, 1.0, 0.25, 5), (-2.0, 2.0, 1.0, 5)],
)
def test_range_count_is_immune_to_float_rounding_and_last_point_is_stop(start, stop, step, n):
    v = inclusive_range(start, stop, step)
    assert v.size == n and v[-1] == stop and v[0] == start
    assert np.allclose(v, start + step * np.arange(n), rtol=0, atol=1e-9 * step)


def test_range_is_not_cumulative_so_it_cannot_drift():
    v = inclusive_range(0.0, 1.0, 1e-4)
    assert v.size == 10001 and v[5000] == pytest.approx(0.5, abs=1e-15)


def test_single_point_range():
    assert inclusive_range(5.0, 5.0, 1.0).tolist() == [5.0]


@pytest.mark.parametrize(
    "args, fragment",
    [
        ((8e9, 12e9, 0.3e9), "integer multiple"),
        ((8e9, 12e9, 0.0), "step must be > 0"),
        ((8e9, 12e9, -1e9), "step must be > 0"),
        ((12e9, 8e9, 1e9), "stop"),
        ((0.0, 1.0, 2.0), "integer multiple"),
        ((0.0, 1e6, 1.0), "exceeds the limit"),
        ((float("nan"), 1.0, 0.5), "finite"),
        ((0.0, float("inf"), 0.5), "finite"),
        ((True, 2.0, 1.0), "finite"),
        (("1", 2.0, 1.0), "finite"),
    ],
)
def test_range_rejects_bad_input(args, fragment):
    with pytest.raises(SweepError, match=fragment):
        inclusive_range(*args)


def test_frequency_grid_validates_and_counts():
    g = FrequencyGrid(8e9, 12e9, 0.5e9)
    assert g.count == 9 and g.frequencies().size == 9 and g.frequencies()[-1] == 12e9
    assert FrequencyGrid.single(10e9).frequencies().tolist() == [10e9]
    assert FrequencyGrid.single(10e9).count == 1
    for bad in ((0.0, 1e9, 1e8), (-1e9, 1e9, 1e8), (1e9, 2e9, 0.3e9), (float("nan"), 1e9, 1e8), (8e9, 12e9, "x")):
        with pytest.raises(SweepError):
            FrequencyGrid(*bad)


def test_frequency_grid_is_frozen():
    g = FrequencyGrid(8e9, 12e9, 1e9)
    with pytest.raises(Exception):
        g.start = 1.0


# ---------------------------------------------------------------- geometry parameter axes


@pytest.mark.parametrize(
    "name", ["period_x", "period_y", "rotation", "length", "type", "thickness", "eps_r", "frequency", "", "Length_m"]
)
def test_only_geometry_parameters_are_sweepable(name):
    with pytest.raises(SweepError, match="unknown geometry parameter"):
        ParameterAxis(name, (1e-3,))


def test_parameter_names_are_the_five_rectangle_fields_plus_the_four_annular_slot_fields():
    # PHASE15 Part E: PARAMETER_NAMES is the UNION of both geometry families' own sweepable
    # parameter names (ParameterAxis itself has no notion of "which geometry kind" -- see
    # sweep.py's own module docstring/comments); a single geometry's own `base_parameters()` only
    # ever returns ITS OWN family's subset, matched by `geometry_from_parameters`'s dispatch.
    assert set(PARAMETER_NAMES) == {
        "length_m",
        "width_m",
        "center_x_m",
        "center_y_m",
        "rotation_deg",
        "a_m",
        "b_m",
        "R_m",
        "r_m",
    }
    assert set(base_parameters(config_from_dict(config_dict()).unit_cell.geometry)) == {
        "length_m",
        "width_m",
        "center_x_m",
        "center_y_m",
        "rotation_deg",
    }


def test_annular_slot_geometry_round_trips_through_base_and_from_parameters():
    from src.geometry import AnnularSlotPatch
    from src.sweep import geometry_from_parameters

    geo = AnnularSlotPatch(a=4e-3, b=4e-3, R=1.5e-3, r=1.0e-3)
    params = base_parameters(geo)
    assert set(params) == {"a_m", "b_m", "R_m", "r_m", "center_x_m", "center_y_m", "rotation_deg"}
    rebuilt = geometry_from_parameters(AnnularSlotPatch.kind, params)
    assert isinstance(rebuilt, AnnularSlotPatch)
    assert rebuilt == geo


@pytest.mark.parametrize(
    "values, fragment",
    [
        ((), "no values"),
        ((1e-3, 1e-3), "duplicate"),
        ((1e-3, float("nan")), "finite"),
        ((1e-3, float("inf")), "finite"),
        ((True,), "finite"),
        (("3 mm",), "finite"),
        ((None,), "finite"),
    ],
)
def test_axis_rejects_bad_values(values, fragment):
    with pytest.raises(SweepError, match=fragment):
        ParameterAxis("length_m", values)


def test_axis_from_range_and_negative_zero_duplicate():
    assert ParameterAxis.from_range("length_m", 2e-3, 4e-3, 1e-3).values == pytest.approx((2e-3, 3e-3, 4e-3))
    with pytest.raises(SweepError, match="duplicate"):
        ParameterAxis("center_x_m", (0.0, -0.0))
    with pytest.raises(SweepError):
        ParameterAxis.from_range("length_m", 2e-3, 4e-3, 0.7e-3)


def test_definition_rejects_duplicate_axes_bad_policy_and_huge_products():
    ax = ParameterAxis("length_m", (3e-3, 4e-3))
    with pytest.raises(SweepError, match="duplicate sweep parameter"):
        SweepDefinition(axes=(ax, ParameterAxis("length_m", (5e-3,))))
    with pytest.raises(SweepError, match="on_invalid"):
        SweepDefinition(axes=(ax,), on_invalid="ignore")
    big = tuple(ParameterAxis(n, tuple(float(i) * 1e-6 for i in range(1, 1001))) for n in PARAMETER_NAMES[:3])
    assert math.prod(len(a.values) for a in big) > MAX_CELLS
    with pytest.raises(SweepError, match="exceeds the limit"):
        SweepDefinition(axes=big)


def test_policy_validation():
    assert RefinementPolicy().enabled is False
    for kwargs in (
        {"factor": 1},
        {"factor": 2.5},
        {"factor": True},
        {"tolerance": 0.0},
        {"tolerance": float("nan")},
        {"frequencies": "some"},
        {"enabled": 1},
    ):
        with pytest.raises(SweepError):
            RefinementPolicy(**kwargs)
    for kwargs in (
        {"min_cells_per_wavelength": 0.0},
        {"max_area_error": -0.1},
        {"strict": "yes"},
        {"min_cells_per_wavelength": float("inf")},
    ):
        with pytest.raises(SweepError):
            ValidationPolicy(**kwargs)


# ------------------------------------------------------------------------- planning cells


def axes(**named):
    return tuple(ParameterAxis(k, tuple(v)) for k, v in named.items())


def test_cell_expansion_is_cartesian_first_axis_slowest_with_unswept_parameters_from_base():
    cfg = sweep_config()
    d = SweepDefinition(axes=axes(length_m=(3e-3, 4e-3, 5e-3), rotation_deg=(0.0, 90.0)))
    plan = plan_sweep(cfg, d)
    assert len(plan.cells) == 6 and plan.n_solves == 6
    got = [(c.parameters["length_m"], c.parameters["rotation_deg"]) for c in plan.cells]
    assert got == [(3e-3, 0.0), (3e-3, 90.0), (4e-3, 0.0), (4e-3, 90.0), (5e-3, 0.0), (5e-3, 90.0)]
    assert [c.index for c in plan.cells] == list(range(6))
    for c in plan.cells:
        assert c.parameters["width_m"] == 2e-3 and c.parameters["center_x_m"] == 0.0 and c.status == "ok"


def test_axis_order_changes_cell_order_not_cell_content():
    cfg = sweep_config()
    a = plan_sweep(cfg, SweepDefinition(axes=axes(length_m=(3e-3, 4e-3), width_m=(1e-3, 2e-3))))
    b = plan_sweep(cfg, SweepDefinition(axes=axes(width_m=(1e-3, 2e-3), length_m=(3e-3, 4e-3))))
    assert {c.cell_id for c in a.cells} == {c.cell_id for c in b.cells}
    assert [c.cell_id for c in a.cells] != [c.cell_id for c in b.cells]


def test_no_axes_gives_exactly_the_base_cell_and_no_frequency_grid_uses_excitation_frequency():
    plan = plan_sweep(sweep_config(), SweepDefinition())
    assert len(plan.cells) == 1 and plan.frequencies.tolist() == [10e9]
    assert plan.cells[0].parameters == base_parameters(plan.config.unit_cell.geometry)


def test_cell_ids_are_deterministic_unique_and_sensitive_to_everything_that_matters():
    d = SweepDefinition(axes=axes(length_m=(3e-3, 4e-3, 5e-3)))
    ids = [c.cell_id for c in plan_sweep(sweep_config(), d).cells]
    assert ids == [c.cell_id for c in plan_sweep(sweep_config(), d).cells]
    assert len(set(ids)) == 3 and all(len(i) == 16 for i in ids)
    for change in (
        {"mesh.cells_x": 8},
        {"mesh.floquet_order": 33},
        {"substrate.eps_r": 2.3},
        {"substrate.tan_delta": 0.001},
        {"substrate.thickness": 1.6e-3},
        {"period_x": 9e-3},
        {"reference_planes.top_air_thickness": 1e-3},
        {"geometry.width": 1e-3},
    ):
        other = [c.cell_id for c in plan_sweep(sweep_config(**change), d).cells]
        assert not set(other) & set(ids), change


def test_unit_cells_are_built_from_parameters_and_rebuildable_from_context():
    cfg = sweep_config()
    plan = plan_sweep(cfg, SweepDefinition(axes=axes(length_m=(3e-3,), rotation_deg=(30.0,), center_x_m=(0.5e-3,))))
    cell = plan.cells[0]
    g = cell.unit_cell.geometry
    assert g.length == 3e-3 and g.width == 2e-3 and g.center_x == 0.5e-3
    assert g.rotation == pytest.approx(math.radians(30.0), rel=1e-15)
    again = unit_cell_from_context(plan.context, cell.parameters)
    assert (
        again.geometry == g and again.lattice == cell.unit_cell.lattice and again.substrate == cell.unit_cell.substrate
    )


def test_unswept_base_geometry_survives_the_degree_round_trip_to_1e15():
    cfg = sweep_config(**{"geometry.rotation_deg": 10.0})
    cell = plan_sweep(cfg, SweepDefinition()).cells[0]
    assert cell.unit_cell.geometry.rotation == pytest.approx(cfg.unit_cell.geometry.rotation, rel=1e-15)


# ------------------------------------------------------- geometry parameter validation

INVALID_CASES = [
    ("feature outside the cell", {"length_m": (9e-3,)}, {}, "extends outside"),
    ("rotated feature leaves the cell", {"rotation_deg": (10.0,)}, {"geometry.length": 8e-3}, "extends outside"),
    ("offset feature leaves the cell", {"center_x_m": (3e-3,)}, {}, "extends outside"),
    ("zero length", {"length_m": (0.0,)}, {}, "must be > 0"),
    ("negative width", {"width_m": (-1e-3,)}, {}, "must be > 0"),
    (
        "dipole wider than long",
        {"width_m": (4e-3,)},
        {"geometry.type": "dipole", "geometry.length": 3e-3},
        "length >= width",
    ),
    (
        "slot wider than long",
        {"width_m": (4e-3,)},
        {"geometry.type": "slot", "geometry.length": 3e-3},
        "length >= width",
    ),
    ("too small for the mesh", {"length_m": (0.2e-3,), "width_m": (0.2e-3,)}, {}, "no interior metal edges"),
]


@pytest.mark.parametrize("label, axis, overrides, fragment", INVALID_CASES, ids=[c[0] for c in INVALID_CASES])
def test_invalid_cells_are_reported_by_name_before_any_solve(monkeypatch, label, axis, overrides, fragment):
    calls = []
    monkeypatch.setattr(runner, "compute_s_parameters", lambda *a, **k: calls.append(1))
    cfg = sweep_config(**overrides)
    d = SweepDefinition(frequency=FrequencyGrid(9e9, 11e9, 1e9), axes=axes(**axis))
    with pytest.raises(SweepError, match=fragment) as exc:
        runner.run_sweep(cfg, definition=d)
    assert "cell 0" in str(exc.value) and "1 of 1" in str(exc.value)
    assert calls == []  # failed before the first solve


def test_error_message_lists_the_offending_parameters_and_counts_all_invalid_cells():
    cfg = sweep_config()
    d = SweepDefinition(axes=axes(length_m=(3e-3, 9e-3, 10e-3, 4e-3)))
    with pytest.raises(SweepError) as exc:
        plan_sweep(cfg, d)
    msg = str(exc.value)
    assert "2 of 4" in msg and "cell 1" in msg and "cell 2" in msg and "length_m=0.009" in msg
    assert "cell 0" not in msg and "cell 3" not in msg


def test_skip_mode_records_invalid_cells_with_reason_and_keeps_the_valid_ones():
    cfg = sweep_config()
    d = SweepDefinition(axes=axes(length_m=(3e-3, 9e-3, 4e-3)), on_invalid="skip")
    plan = plan_sweep(cfg, d)
    assert [c.status for c in plan.cells] == ["ok", "invalid", "ok"]
    bad = plan.cells[1]
    assert "extends outside" in bad.reason and bad.unit_cell is None and bad.checks is None
    assert len(plan.valid_cells) == 2 and plan.n_solves == 2
    assert bad.cell_id and bad.parameters["length_m"] == 9e-3


def test_boundary_cases_that_are_valid():
    """Touching the cell edge is allowed (Phase 1 rule); slot up to full width is a valid complement."""
    cfg = sweep_config()
    ok = plan_sweep(cfg, SweepDefinition(axes=axes(length_m=(8e-3,), width_m=(2e-3,))))
    assert ok.cells[0].status == "ok"
    slot = sweep_config(**{"geometry.type": "slot", "geometry.length": 4e-3, "geometry.width": 2e-3})
    assert plan_sweep(slot, SweepDefinition(axes=axes(length_m=(4e-3, 6e-3)))).cells[1].status == "ok"


def test_plan_without_definition_is_an_error():
    with pytest.raises(SweepError, match="no sweep definition"):
        plan_sweep(sweep_config())


# ------------------------------------------------------------- strict a priori validation


def test_strict_mode_fails_before_solving_when_the_mesh_is_too_coarse_for_the_wavelength():
    cfg = sweep_config()
    d = SweepDefinition(
        frequency=FrequencyGrid(8e9, 12e9, 4e9),
        axes=axes(length_m=(3e-3,)),
        validation=ValidationPolicy(min_cells_per_wavelength=60.0, strict=True),
    )
    with pytest.raises(SweepError, match="cells/wavelength"):
        plan_sweep(cfg, d)
    relaxed = SweepDefinition(
        frequency=d.frequency, axes=d.axes, validation=ValidationPolicy(min_cells_per_wavelength=60.0)
    )
    assert plan_sweep(cfg, relaxed).n_solves == 2  # 8 and 12 GHz


def test_strict_mode_rejects_grating_frequencies_and_misaligned_features():
    cfg = sweep_config()
    f_c = C0 / 8e-3
    grating = SweepDefinition(
        frequency=FrequencyGrid(f_c * 0.9, f_c * 1.1, f_c * 0.1), validation=ValidationPolicy(strict=True)
    )
    with pytest.raises(SweepError, match="grating order"):
        plan_sweep(cfg, grating)
    misaligned = SweepDefinition(axes=axes(length_m=(5.3e-3,)), validation=ValidationPolicy(strict=True))
    with pytest.raises(SweepError, match="metal area error"):
        plan_sweep(cfg, misaligned)
    assert plan_sweep(cfg, SweepDefinition(axes=axes(length_m=(5.3e-3,)))).cells[0].checks["area_ok"] is False


# ------------------------------------------------------------------- config file section


def test_config_without_sweep_section_has_none_and_shipped_config_is_unchanged():
    assert config_from_dict(config_dict()).sweep is None
    assert load_config(ROOT / "configs" / "config.yaml").sweep is None


def test_sweep_section_parses_units_ranges_lists_and_policy():
    cfg = sweep_config(
        {
            "frequency": FREQ_8_12,
            "parameters": {
                "length": {"start": "3 mm", "stop": "5 mm", "step": "1 mm"},
                "rotation_deg": [0, 45.5],
                "width": {"values": ["1 mm", 2.0e-3]},
            },
            "on_invalid": "skip",
            "validation": {
                "min_cells_per_wavelength": 12,
                "max_area_error": 0.1,
                "strict": True,
                "refinement": {"enabled": True, "factor": 3, "tolerance": 0.01, "frequencies": "all"},
            },
        }
    )
    sw = cfg.sweep
    assert (sw.frequency.start, sw.frequency.stop, sw.frequency.step) == (8e9, 12e9, 1e9)
    assert [a.name for a in sw.axes] == ["length_m", "rotation_deg", "width_m"]  # YAML order kept
    assert sw.axes[0].values == pytest.approx((3e-3, 4e-3, 5e-3)) and sw.axes[1].values == (0.0, 45.5)
    assert sw.axes[2].values == (1e-3, 2e-3)
    assert sw.on_invalid == "skip" and sw.validation.strict is True
    r = sw.validation.refinement
    assert (r.enabled, r.factor, r.tolerance, r.frequencies) == (True, 3, 0.01, "all")
    assert sw.validation.min_cells_per_wavelength == 12.0 and sw.validation.max_area_error == 0.1


def test_frequency_only_and_parameter_only_sweeps():
    f_only = sweep_config({"frequency": FREQ_8_12}).sweep
    assert f_only.axes == () and f_only.frequency.count == 5
    p_only = sweep_config({"parameters": {"length": {"values": ["3 mm"]}}}).sweep
    assert p_only.frequency is None
    assert sweep_config({}).sweep == SweepDefinition()


@pytest.mark.parametrize(
    "sweep, fragment",
    [
        ({"frequncy": FREQ_8_12}, "unknown key"),
        ({"frequency": {"start": "8 GHz", "stop": "12 GHz"}}, "missing required"),
        ({"frequency": {**FREQ_8_12, "count": 5}}, "unknown key"),
        ({"frequency": {"start": "8 GHz", "stop": "12 GHz", "step": "0.3 GHz"}}, "integer multiple"),
        ({"frequency": {"start": "8 furlong", "stop": "12 GHz", "step": "1 GHz"}}, "unknown unit"),
        ({"frequency": {"start": "12 GHz", "stop": "8 GHz", "step": "1 GHz"}}, "stop"),
        ({"parameters": {"period_x": [8e-3]}}, "unknown geometry parameter"),
        ({"parameters": {"type": ["dipole"]}}, "unknown geometry parameter"),
        ({"parameters": {"rotation": [0.0]}}, "unknown geometry parameter"),
        (
            {"parameters": {"length": {"values": [3e-3], "start": 1e-3, "stop": 2e-3, "step": 1e-3}}},
            "not both or neither",
        ),
        ({"parameters": {"length": {}}}, "not both or neither"),
        ({"parameters": {"length": {"start": 1e-3, "stop": 2e-3}}}, "missing required"),
        ({"parameters": {"length": {"start": "3", "stop": "5", "step": "1"}}}, ">= 1.0 m"),
        ({"parameters": {"length": {"values": "3 mm"}}}, "must be a list"),
        ({"parameters": {"length": [3e-3, 3e-3]}}, "duplicate"),
        ({"parameters": {"length": []}}, "no values"),
        ({"parameters": {"rotation_deg": ["45 deg"]}}, "finite number"),
        ({"parameters": {"rotation_deg": [True]}}, "finite number"),
        ({"parameters": [1, 2]}, "must be a mapping"),
        ({"on_invalid": "ignore"}, "on_invalid"),
        ({"validation": {"strict": "yes"}}, "boolean"),
        ({"validation": {"tolerance": 0.1}}, "unknown key"),
        ({"validation": {"refinement": {"factor": 1}}}, "factor"),
        ({"validation": {"refinement": {"frequencies": "min"}}}, "frequencies"),
        ({"validation": {"refinement": {"enabled": True, "tolerance": -1}}}, "tolerance"),
        ("frequency", "must be a mapping"),
    ],
)
def test_bad_sweep_sections_raise_config_error_with_a_useful_message(sweep, fragment):
    with pytest.raises(ConfigError, match=fragment):
        sweep_config(sweep)


def test_axis_values_outside_the_cell_are_not_a_config_error_but_a_planning_error():
    """Parsing checks the *form* of a value; the geometry rules are applied per cell by the plan."""
    cfg = sweep_config({"parameters": {"length": {"values": ["9 mm"]}}})
    with pytest.raises(SweepError, match="extends outside"):
        plan_sweep(cfg)


def test_example_sweep_config_loads_and_plans():
    cfg = load_config(ROOT / "examples" / "sweep_config.yaml")
    plan = plan_sweep(cfg)
    assert plan.n_solves == len(plan.valid_cells) * len(plan.frequencies) > 0
