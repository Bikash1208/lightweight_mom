"""Grid options of `scripts/research/phase20_r4_build_library.py`: defaults reproduce the R4 constants exactly; the approved full grid is 41 x 217 points."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))

import phase20_r4_build_library as b  # noqa: E402

OLD_F = [8.0 + 0.25 * i for i in range(17)]
OLD_THETA = [7.5 * i for i in range(7)]
OLD_PHI = [30.0 * i for i in range(12)]


def _restore():
    b.GRID_F, b.GRID_THETA, b.GRID_PHI = list(OLD_F), list(OLD_THETA), list(OLD_PHI)


def test_default_constants_are_the_r4_grid():
    assert (b.GRID_F, b.GRID_THETA, b.GRID_PHI) == (OLD_F, OLD_THETA, OLD_PHI)
    assert len(b.grid_points()) == 73


def test_configure_grid_with_the_default_steps_changes_nothing():
    try:
        b.configure_grid()
        assert (b.GRID_F, b.GRID_THETA, b.GRID_PHI) == (OLD_F, OLD_THETA, OLD_PHI)
        b.configure_grid(f_step=0.25, theta_step=7.5, phi_step=30.0)
        assert (b.GRID_F, b.GRID_THETA, b.GRID_PHI) == (OLD_F, OLD_THETA, OLD_PHI)
    finally:
        _restore()


def test_full_approved_grid_is_41_by_217_points():
    try:
        b.configure_grid(f_step=0.1, theta_step=5.0, phi_step=15.0)
        assert len(b.GRID_F) == 41 and b.GRID_F[0] == 8.0 and abs(b.GRID_F[-1] - 12.0) < 1e-9
        assert b.GRID_THETA == [5.0 * i for i in range(10)]
        assert b.GRID_PHI == [15.0 * i for i in range(24)]
        assert len(b.grid_points()) == 1 + 9 * 24
    finally:
        _restore()


def test_option_parser_keeps_the_positional_arguments():
    pos, kw, dry = b.parse_grid_options(["x.py", "grid", "--f-step", "0.1", "40", "--dry-run", "reverse", "--theta-list", "0,5,10"])
    assert pos == ["x.py", "grid", "40", "reverse"] and dry is True
    assert kw == {"f_step": 0.1, "theta_list": [0.0, 5.0, 10.0]}
