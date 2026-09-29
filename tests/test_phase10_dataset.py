"""Sec. 18-21, 47: dataset harvesting/generation, symmetry-safe grouping, and the geometry-based
holdout."""

from src.materials import Metal, Substrate
from src.phase10 import dataset as ds_mod
from src.phase10 import environment_geometry as eg
from src.phase10.dataset import ALL_POSITIONS, EnvironmentSample
from src.phase10.neighborhood_study import run_full_neighborhood_study

SUB = Substrate("Rogers RT/Duroid 5880", 2.2, 0.0009, 1.575e-3)
MET = Metal()
P = 4.8e-3
F = 12.5e9


def _sample(sample_id, center, neighbors, center_name="C3"):
    return EnvironmentSample(sample_id, center_name, center, neighbors, "test", F, 1.0 + 0j, "test")


def test_harvest_structured_samples_from_a_real_small_study():
    study = run_full_neighborhood_study({"C3": (2.8, 2.8)}, P, SUB, MET, F, mesh_cells_per_subcell=8)
    samples = ds_mod.harvest_structured_samples(study, F)
    assert len(samples) > 0
    for s in samples:
        assert set(s.neighbors_lw_mm.keys()) == set(ALL_POSITIONS)
        assert s.source == "structured"


def test_generate_random_samples_produces_solved_full_9cell_environments():
    samples = ds_mod.generate_random_samples(eg.CENTER_CELLS_MM, 3, P, SUB, MET, F, mesh_cells_per_subcell=8, seed=42)
    assert len(samples) == 3
    for s in samples:
        assert s.source == "random"
        assert isinstance(s.s_complex, complex)
        for pos in ALL_POSITIONS:
            l, w = s.neighbors_lw_mm[pos]
            assert eg.in_domain(l, w)


def test_generate_random_samples_is_reproducible_with_a_fixed_seed():
    a = ds_mod.generate_random_samples(eg.CENTER_CELLS_MM, 2, P, SUB, MET, F, mesh_cells_per_subcell=8, seed=7)
    b = ds_mod.generate_random_samples(eg.CENTER_CELLS_MM, 2, P, SUB, MET, F, mesh_cells_per_subcell=8, seed=7)
    assert [s.center_lw_mm for s in a] == [s.center_lw_mm for s in b]
    assert [s.neighbors_lw_mm for s in a] == [s.neighbors_lw_mm for s in b]


def test_symmetry_canonical_key_groups_y_reflected_environments_together():
    center = (2.8, 2.8)
    neighbors_a = {"N": (3.0, 3.0), "E": (2.8, 2.8), "S": (1.0, 1.0), "W": (2.8, 2.8),
                   "NE": (2.8, 2.8), "NW": (2.8, 2.8), "SE": (2.8, 2.8), "SW": (2.8, 2.8)}
    neighbors_b = {"N": (1.0, 1.0), "E": (2.8, 2.8), "S": (3.0, 3.0), "W": (2.8, 2.8),   # N<->S swapped
                   "NE": (2.8, 2.8), "NW": (2.8, 2.8), "SE": (2.8, 2.8), "SW": (2.8, 2.8)}
    sa = _sample("a", center, neighbors_a)
    sb = _sample("b", center, neighbors_b)
    assert ds_mod.symmetry_canonical_key(sa) == ds_mod.symmetry_canonical_key(sb)


def test_symmetry_canonical_key_distinguishes_genuinely_different_environments():
    center = (2.8, 2.8)
    neighbors_a = {p: center for p in ALL_POSITIONS}
    neighbors_c = dict(neighbors_a)
    neighbors_c["N"] = (4.4, 4.4)
    sa = _sample("a", center, neighbors_a)
    sc = _sample("c", center, neighbors_c)
    assert ds_mod.symmetry_canonical_key(sa) != ds_mod.symmetry_canonical_key(sc)


def test_assign_splits_places_every_holdout_center_sample_in_blind_test():
    center = (2.8, 2.8)
    holdout_center = (4.2, 4.2)
    samples = [
        _sample("keep1", center, {p: center for p in ALL_POSITIONS}, "C3"),
        _sample("hold1", holdout_center, {p: holdout_center for p in ALL_POSITIONS}, "C5"),
        _sample("hold2", holdout_center, {p: holdout_center for p in ALL_POSITIONS}, "C5"),
    ]
    splits = ds_mod.assign_splits(samples, holdout_center_names=("C5",))
    assert splits["hold1"] == "blind_test"
    assert splits["hold2"] == "blind_test"


def test_assign_splits_never_separates_a_symmetry_equivalent_group():
    """N-only and S-only single-cardinal perturbations of the same magnitude are y-reflection
    images of each other (``MIRROR_NS``) and MUST share a canonical key -- and therefore a split,
    however ``assign_splits``'s shuffle falls."""
    center = (2.8, 2.8)
    n_pert = {"N": (3.0, 3.0), "E": (2.8, 2.8), "S": (2.8, 2.8), "W": (2.8, 2.8),
              "NE": (2.8, 2.8), "NW": (2.8, 2.8), "SE": (2.8, 2.8), "SW": (2.8, 2.8)}
    s_pert = {"N": (2.8, 2.8), "E": (2.8, 2.8), "S": (3.0, 3.0), "W": (2.8, 2.8),
              "NE": (2.8, 2.8), "NW": (2.8, 2.8), "SE": (2.8, 2.8), "SW": (2.8, 2.8)}
    assert ds_mod.symmetry_canonical_key(_sample("x", center, n_pert)) == \
        ds_mod.symmetry_canonical_key(_sample("y", center, s_pert))

    # Pad with enough unrelated, mutually distinct samples that the group shuffle is nontrivial.
    samples = [_sample("n_pert", center, n_pert, "C3"), _sample("s_pert", center, s_pert, "C3")]
    for i in range(30):
        distinct = {p: center for p in ALL_POSITIONS}
        distinct["N"] = (1.0 + 0.1 * i, 1.0 + 0.1 * i)
        samples.append(_sample(f"filler{i}", center, distinct, "C3"))
    splits = ds_mod.assign_splits(samples, holdout_center_names=())
    assert splits["n_pert"] == splits["s_pert"]
