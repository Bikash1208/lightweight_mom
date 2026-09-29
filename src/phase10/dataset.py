"""Sec. 18-21: the environment-response dataset. Every sample is a real, solved 3x3 heterogeneous
supercell (Sec. 13's full 9-cell representation -- justified directly by ``neighborhood_decision``'s
own measured finding that a cardinal-only, 5-cell model is NOT sufficient for this library, so the
dataset is built at the neighborhood size the evidence actually supports, not the cheaper one).

Two sources, both real MoM solves, neither fabricated:
  "structured" -- harvested directly from the already-computed `neighborhood_study` results (no
    extra solves): homogeneous, single-neighbor, opposite-pair, four-neighbor, asymmetric,
    diagonal-single, 9-cell-uniform and correlated-environment cases, for every center cell.
  "random"     -- additional full 9-cell environments (all 8 neighbor positions independently
    drawn) from a fixed-seed generator, to cover genuinely simultaneous multi-neighbor contrast
    that the structured cases (which vary at most 4 positions in controlled patterns) do not
    reach. Sec. 19's explicit coverage categories (homogeneous / weak / moderate / strong /
    symmetric / asymmetric / anisotropic / extremes / midpoints) are satisfied by the structured
    half; the random half adds combinatorial diversity within the same validated domain.
"""

from dataclasses import dataclass

import numpy as np

from . import environment_geometry as eg
from .grid_solver import solve_grid

ALL_POSITIONS = ("N", "E", "S", "W", "NE", "NW", "SE", "SW")
MIRROR_EW = {"N": "N", "S": "S", "E": "W", "W": "E", "NE": "NW", "NW": "NE", "SE": "SW", "SW": "SE"}
MIRROR_NS = {"N": "S", "S": "N", "E": "E", "W": "W", "NE": "SE", "SE": "NE", "NW": "SW", "SW": "NW"}
MIRROR_180 = {"N": "S", "S": "N", "E": "W", "W": "E", "NE": "SW", "SW": "NE", "NW": "SE", "SE": "NW"}

RANDOM_SEED = 1729  # fixed and recorded (Sec. 21/Rule 4's own spirit of reproducibility)


@dataclass(frozen=True)
class EnvironmentSample:
    sample_id: str
    center_name: str
    center_lw_mm: tuple
    neighbors_lw_mm: dict  # position -> (L_mm, W_mm), all 8 positions present
    category: str
    frequency_hz: float
    s_complex: complex
    source: str  # "structured" | "random"


def _grid_from_neighbors(center_lw_mm, neighbors_lw_mm):
    cardinal = {d: neighbors_lw_mm[d] for d in ("N", "E", "S", "W")}
    diagonal = {d: neighbors_lw_mm[d] for d in ("NE", "NW", "SE", "SW")}
    return eg.build_9cell_environment(center_lw_mm, cardinal_mm=cardinal, diagonal_mm=diagonal)


def harvest_structured_samples(study_results, frequency_hz):
    """Pulls every case already solved by ``neighborhood_study.run_full_neighborhood_study`` into
    dataset samples -- zero additional MoM solves."""
    samples = []
    for center_name, r in study_results.items():
        c = r.center_lw_mm
        # 5A/5B/5C/5D/5E (cardinal-only, diagonals == center) for every contrast label.
        for label, cases in r.five_cell_cases.items():
            for case_name, case in cases.items():
                neighbors = {d: c for d in ALL_POSITIONS}
                neighbors.update(case.cardinal_lw_mm)
                samples.append(
                    EnvironmentSample(
                        f"{center_name}_{label}_{case_name}",
                        center_name,
                        c,
                        neighbors,
                        f"structured_{case_name}_{label}",
                        frequency_hz,
                        case.px,
                        "structured",
                    )
                )
        # D1 diagonal-single cases.
        for label, cases in r.diagonal_cases.items():
            for case_name, case in cases.items():
                neighbors = {d: c for d in ALL_POSITIONS}
                neighbors.update(case.diagonal_lw_mm)
                samples.append(
                    EnvironmentSample(
                        f"{center_name}_{label}_{case_name}",
                        center_name,
                        c,
                        neighbors,
                        f"structured_{case_name}_{label}",
                        frequency_hz,
                        case.px,
                        "structured",
                    )
                )
        # 9-cell-uniform per contrast label.
        for label, px in r.nine_cell_true_px.items():
            contrast = r.contrast_values_mm[label]
            neighbors = {d: contrast for d in ALL_POSITIONS}
            samples.append(
                EnvironmentSample(
                    f"{center_name}_9cell_uniform_{label}",
                    center_name,
                    c,
                    neighbors,
                    f"structured_9cell_uniform_{label}",
                    frequency_hz,
                    px,
                    "structured",
                )
            )
        # correlated environment (2 patterns).
        low_high = r.correlated
        diag_a, card_b = (
            {d: low_high.a_lw_mm for d in ("NE", "NW", "SE", "SW")},
            {d: low_high.b_lw_mm for d in ("N", "E", "S", "W")},
        )
        neighbors1 = {**card_b, **diag_a}
        samples.append(
            EnvironmentSample(
                f"{center_name}_correlated_pattern1",
                center_name,
                c,
                neighbors1,
                "structured_correlated",
                frequency_hz,
                low_high.pattern_1_px,
                "structured",
            )
        )
        diag_b, card_a = (
            {d: low_high.b_lw_mm for d in ("NE", "NW", "SE", "SW")},
            {d: low_high.a_lw_mm for d in ("N", "E", "S", "W")},
        )
        neighbors2 = {**card_a, **diag_b}
        samples.append(
            EnvironmentSample(
                f"{center_name}_correlated_pattern2",
                center_name,
                c,
                neighbors2,
                "structured_correlated",
                frequency_hz,
                low_high.pattern_2_px,
                "structured",
            )
        )
    return samples


def generate_random_samples(
    centers_mm, n_samples, period_m, substrate, metal, frequency_hz, mesh_cells_per_subcell=8, seed=RANDOM_SEED
):
    """Sec. 19: fixed-seed (recorded, ``RANDOM_SEED``) combinatorial environments -- every one of
    the 8 neighbor positions independently drawn (uniformly, on the same validated 0.2 mm sweep
    grid), and the center drawn uniformly from the given center cells. Each is a genuine, freshly
    solved 3x3 heterogeneous-supercell MoM case (no fabricated values)."""
    rng = np.random.default_rng(seed)
    grid_values_mm = np.round(np.arange(eg.DOMAIN_L_MIN_MM, eg.DOMAIN_L_MAX_MM + 1e-9, eg.DOMAIN_STEP_MM), 6)
    center_names = list(centers_mm.keys())
    samples = []
    for k in range(n_samples):
        center_name = center_names[rng.integers(0, len(center_names))]
        center_lw_mm = centers_mm[center_name]
        neighbors_lw_mm = {
            d: (float(rng.choice(grid_values_mm)), float(rng.choice(grid_values_mm))) for d in ALL_POSITIONS
        }
        grid_mm = _grid_from_neighbors(center_lw_mm, neighbors_lw_mm)
        grid_m = eg.grid_to_m(grid_mm)
        sg = solve_grid(3, period_m, substrate, metal, grid_m, frequency_hz, mesh_cells_per_subcell)
        px = sg.dipole_at(1, 1).px
        samples.append(
            EnvironmentSample(
                f"random_{k:04d}",
                center_name,
                center_lw_mm,
                neighbors_lw_mm,
                "random_combinatorial",
                frequency_hz,
                px,
                "random",
            )
        )
    return samples


def _mirror_neighbors(neighbors_lw_mm, mapping):
    return {mapping[d]: neighbors_lw_mm[d] for d in ALL_POSITIONS}


def symmetry_canonical_key(sample: EnvironmentSample):
    """Rule 4: a sample must not be able to leak into both TRAIN and TEST via a physically
    equivalent mirrored arrangement. Conservative treatment (documented, not a numerical-equality
    claim): the 3x3 environment's y-reflection, x-reflection and 180-degree-rotation images are
    all mapped to ONE canonical key (lexicographically smallest), so any sample related to another
    by one of these reflections is grouped identically and can never be split across sets. This
    does not assert that mirrored configurations give numerically identical responses under this
    project's fixed x-polarized excitation convention -- it only prevents them from being treated
    as independent evidence for train/test isolation purposes."""
    variants = [dict(sample.neighbors_lw_mm)]
    for mapping in (MIRROR_EW, MIRROR_NS, MIRROR_180):
        variants.append(_mirror_neighbors(sample.neighbors_lw_mm, mapping))
    variant_keys = [tuple(sorted(v.items())) for v in variants]
    canonical = min(variant_keys)
    return (sample.center_lw_mm, canonical)


def assign_splits(samples, holdout_center_names=("C5",), val_fraction=0.15, test_fraction=0.15, seed=RANDOM_SEED):
    """Sec. 20-21. Two mechanisms, both applied:
    1. Sec. 21's geometry-based holdout: every sample whose center cell is in
       ``holdout_center_names`` is placed in the blind test set, unconditionally -- the model
       never sees this center cell's environment responses during training, testing genuine
       extrapolation to an unseen geometry region.
    2. Rule 4's symmetry-safe split: the remaining samples are grouped by
       ``symmetry_canonical_key`` and whole GROUPS (never individual samples) are assigned to
       train/validation/blind_test, using a fixed-seed shuffle of the group list, so a group
       cannot straddle two sets."""
    rng = np.random.default_rng(seed)
    assignment = {}
    holdout_ids, remaining = [], []
    for s in samples:
        if s.center_name in holdout_center_names:
            holdout_ids.append(s.sample_id)
        else:
            remaining.append(s)

    groups = {}
    for s in remaining:
        groups.setdefault(symmetry_canonical_key(s), []).append(s.sample_id)
    group_keys = list(groups.keys())
    rng.shuffle(group_keys)

    n_groups = len(group_keys)
    n_val = max(1, int(round(n_groups * val_fraction))) if n_groups > 0 else 0
    n_test = max(1, int(round(n_groups * test_fraction))) if n_groups > 0 else 0
    val_keys = set(group_keys[:n_val])
    test_keys = set(group_keys[n_val : n_val + n_test])

    for key in group_keys:
        split = "validation" if key in val_keys else ("blind_test" if key in test_keys else "train")
        for sid in groups[key]:
            assignment[sid] = split
    for sid in holdout_ids:
        assignment[sid] = "blind_test"
    return assignment
