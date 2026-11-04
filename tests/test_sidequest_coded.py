"""Side quest: grouped coded hologram model (scripts/research/sidequest_coded_model.py) — combinatorial / budget checks (no far-field evaluation)."""
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
warnings.simplefilter("ignore")
import sidequest_coded_model as M  # noqa: E402


@pytest.fixture(scope="module")
def k2():
    return M.Board([0, 180], "sector")


@pytest.fixture(scope="module")
def k4():
    return M.Board([0, 90, 180, 270], "sector")


def test_ideal_codes_static_fraction_is_chance_level(k2, k4):
    s2 = M.summarize(k2, k2.ideal_bits(np.zeros(2)))
    s4 = M.summarize(k4, k4.ideal_bits(np.zeros(4)))
    assert 0.44 < s2["static_frac"] < 0.54  # two uncorrelated 1-bit codes agree on ~50 % of the cells
    assert 0.08 < s4["static_frac"] < 0.22  # four: 2 / 16 = 12.5 %
    assert s2["lines_with_inverters"] == 1 and s2["lines_no_inverters"] == 2
    assert s4["lines_with_inverters"] <= 7 and s4["lines_no_inverters"] <= 14  # <= 2^K - 2 words, complementary pairs share a line


def test_greedy_respects_the_amplitude_budget_and_static_cells_are_identical(k2):
    eps = np.zeros(2)
    for budget in (0.25, 1.0):
        bits, s0, r, a0 = M.greedy(k2, eps, budget, 0.5, merge=True)
        assert np.all(r <= 1 - 10 ** (-budget / 20) + 1e-9)
        c = k2.contrib(eps)
        a_new = np.array([(c[k] * bits[k] * s0[k]).sum() for k in range(2)])
        assert np.all(a_new / a0 >= 10 ** (-budget / 20) - 1e-6)  # the coherent peak amplitude stays within the budget
        static = np.all(bits == bits[0], axis=0)
        assert np.all(bits[0][static] == bits[1][static])


def test_more_budget_never_reduces_the_static_fraction(k2):
    eps = np.zeros(2)
    f = [M.summarize(k2, M.greedy(k2, eps, b, 1.0, merge=False)[0])["static_frac"] for b in (0.0, 0.25, 0.5, 1.0)]
    assert f[0] <= f[1] <= f[2] <= f[3] and f[3] > f[0] + 0.2


def test_complement_classes():
    words = np.array([[1.0, -1.0], [-1.0, 1.0], [1.0, 1.0]])
    assert len(M.classes(words)) == 2  # {+-, -+} share a line, ++ is its own word
