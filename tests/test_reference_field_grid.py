"""Part 3: ReferenceFieldGrid construction, grid regularity, local-wavenumber estimation and the
TM0 sanity check, validated on the real example .fld file. The radial beta plateau this module
measures (~225 rad/m over r=10-140mm) independently reproduces the brief's OWN cited number,
confirming the 2-D unwrap + gradient approach works on this project's real data."""

from pathlib import Path

import numpy as np
import pytest

from src.dispersion.analytic_tm0 import tm0_beta
from src.fields import fld_parser as fp
from src.fields import reference_field_grid as rfg

EXAMPLE_PATH = Path(__file__).resolve().parent.parent / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
pytestmark = pytest.mark.skipif(not EXAMPLE_PATH.exists(), reason="example .fld file not present in this checkout")

EPS_R = 2.2
THICKNESS_M = 1.575e-3


@pytest.fixture(scope="module")
def example_grid():
    return fp.read_fld(EXAMPLE_PATH)


def test_grid_is_regular(example_grid):
    report = rfg.check_grid_regularity(example_grid)
    assert report.regular
    assert report.dx_m == pytest.approx(2e-3, abs=1e-9)
    assert report.dy_m == pytest.approx(2e-3, abs=1e-9)
    assert report.nx == 151 and report.ny == 151


def test_fld_to_reference_field_grid_phase_only(example_grid):
    grid = rfg.fld_to_reference_field_grid(example_grid)
    assert grid.amplitude_available is False
    assert grid.Ex is None
    assert grid.phase_deg.shape == (151, 151)
    assert grid.dx_m == pytest.approx(2e-3, abs=1e-9)
    assert grid.z_reference_m == pytest.approx(1.7e-3, abs=1e-9)
    # phase() falls back to phase_deg (equal_nan: the feed point is NaN in both)
    assert np.array_equal(grid.phase(), grid.phase_deg, equal_nan=True)


def test_rejects_non_phase_expression(example_grid):
    fake_mag = fp.FldGrid(
        header=fp.FldHeader(
            min_m=example_grid.header.min_m,
            max_m=example_grid.header.max_m,
            step_m=example_grid.header.step_m,
            coordinate_system="x",
            unit="mm",
            expression=fp.parse_expression("CmplxMag(ScalarZ(<Ex,Ey,Ez>))"),
        ),
        x_m=example_grid.x_m,
        y_m=example_grid.y_m,
        value=np.ones_like(example_grid.value),
        ordering=example_grid.ordering,
        path="fake",
    )
    with pytest.raises(ValueError, match="must be a CmplxPhase"):
        rfg.fld_to_reference_field_grid(fake_mag)


def test_fld_to_reference_field_grid_with_magnitude(example_grid):
    magnitude_fld = fp.FldGrid(
        header=example_grid.header,
        x_m=example_grid.x_m,
        y_m=example_grid.y_m,
        value=np.ones_like(example_grid.value),
        ordering=example_grid.ordering,
        path="synthetic_mag",
    )
    grid = rfg.fld_to_reference_field_grid(example_grid, magnitude_fld)
    assert grid.amplitude_available is True
    assert np.allclose(np.abs(grid.Ex[~np.isnan(grid.Ex)]), 1.0)


def test_radial_beta_plateau_matches_brief_cited_value(example_grid):
    """The headline empirical check: this module's own beta(r) estimate reproduces the brief's
    cited ~-225 rad/m (magnitude 225) plateau over r = 10-140mm, independently."""
    beta = rfg.local_wavenumber_field(example_grid)
    r = np.hypot(example_grid.x_m, example_grid.y_m)
    mask = (r > 0.010) & (r < 0.140)
    median_beta = float(np.median(beta[mask]))
    assert median_beta == pytest.approx(225.0, abs=5.0)


def test_radial_beta_profile_is_roughly_flat_in_the_plateau_region(example_grid):
    centers, medians = rfg.radial_beta_profile(example_grid, n_bins=20)
    in_plateau = (centers > 0.015) & (centers < 0.135)
    plateau_values = medians[in_plateau]
    plateau_values = plateau_values[~np.isnan(plateau_values)]
    assert plateau_values.std() / plateau_values.mean() < 0.15  # "constant" per the brief


def test_compare_beta_to_tm0_reports_a_real_relative_difference(example_grid):
    beta = rfg.local_wavenumber_field(example_grid)
    r = np.hypot(example_grid.x_m, example_grid.y_m)
    median_beta = float(np.median(beta[(r > 0.010) & (r < 0.140)]))
    comparison = rfg.compare_beta_to_tm0(median_beta, 10e9, EPS_R, THICKNESS_M)
    assert comparison.beta_tm0 == pytest.approx(tm0_beta(10e9, EPS_R, THICKNESS_M))
    assert comparison.n_eff_measured > 1.0
    # the brief flags this as a guided/surface wave, not the bare-slab TM0 itself -- expect a
    # real, nonzero difference (this is a sanity check, not an equality claim)
    assert comparison.relative_difference > 0.0


def test_cells_per_fringe_matches_the_briefs_worked_example():
    cpf = rfg.cells_per_fringe(225.0, period_m=8.0e-3)
    assert cpf == pytest.approx(3.5, abs=0.1)
    assert rfg.aliasing_warning(cpf) is not None


def test_aliasing_warning_absent_when_well_sampled():
    cpf = rfg.cells_per_fringe(225.0, period_m=1.0e-3)
    assert rfg.aliasing_warning(cpf) is None


def test_nan_at_feed_point_does_not_propagate_through_the_whole_row():
    grid = None
    grid = fp.read_fld(EXAMPLE_PATH)
    beta = rfg.local_wavenumber_field(grid)
    # far from the feed point, beta should be finite and near the plateau value
    assert np.isfinite(beta[0, 0])
    assert np.isfinite(beta[-1, -1])


def test_cell_centers_whole_cell_side_keeps_all_cells():
    """Phase 21 C1: 300.8 mm / 6.4 mm is 46.999... in floating point and used to floor to 46."""
    centers, n = rfg._cell_centers(-0.1504, 0.1504, 6.4e-3)
    assert n == 47 and len(centers) == 47
    assert centers[0] == pytest.approx(-0.1504 + 3.2e-3)
    # Also for an off-origin span and for the 300 mm / 8 mm docstring example (37, not 38).
    assert rfg._cell_centers(0.0, 0.3008, 6.4e-3)[1] == 47
    assert rfg._cell_centers(-0.15, 0.15, 8e-3)[1] == 37
