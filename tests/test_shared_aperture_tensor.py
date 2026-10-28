"""Phase 29 angle-dependent tensor model (MODEL_ONLY): LUT provenance / symmetry and design covariance (translation by whole lattice pitches, 180 deg rotation)."""

import json
from pathlib import Path

import numpy as np
import pytest

from src.hologram import shared_aperture_tensor as st
from src.hologram import tensor_cp_design as t
from src.hologram.tensor_leaky_cp import SupercellTable

RES = Path(__file__).resolve().parents[1] / "results/phase29"
pytestmark = pytest.mark.skipif(not (RES / "angle_lut_v1.npz").exists(), reason="LUT not generated")


@pytest.fixture(scope="module")
def table():
    return SupercellTable.from_files(t.TABLE_GLOB)


@pytest.fixture(scope="module")
def lut(table):
    return st.AngleLUT(RES / "angle_lut_v1.npz", table)


def test_lut_provenance_and_symmetry():
    meta = json.loads((RES / "angle_lut_v1.json").read_text())
    assert (
        meta["version"] == "v1"
        and "MODEL_ONLY" in meta["labels"]
        and meta["kt_over_k0"] == 1.065
        and meta["nan_beta_count"] == 0
    )
    assert (
        meta["symmetry_checks"]["c4_relative_error_median"] < 1e-9
        and meta["symmetry_checks"]["c4_relative_error_max"] < 0.1
    )
    d = np.load(RES / "angle_lut_v1.npz")
    assert d["Y"].shape == (4, 12, 12, 2, 2) and np.isfinite(d["Y"]).all()


def test_lut_calibration_matches_the_table_strength_and_is_periodic(table, lut):
    strength = np.sqrt(np.abs(lut.dkk) ** 2 + np.abs(lut.dsk) ** 2).mean(axis=(1, 2)) * lut.c
    assert strength == pytest.approx(np.nanmean(np.sqrt(2.0 * table.alpha), axis=1), rel=1e-9)
    a, b, _ = lut.tables_at(np.array([0.3, 0.3 + np.pi, 0.3 - 2 * np.pi]))
    assert np.allclose(a[..., 0], a[..., 1]) and np.allclose(a[..., 0], a[..., 2])  # psi has period pi


def test_translation_and_rotation_covariance_of_the_design(table, lut):
    def run(ports, centre=(0.0, 0.0)):
        d = st.design_tensor(table, lut, 45.0, ports, [0.0] * len(ports), centre_mm=centre, n_iter=1, n_psi=120)
        return d, sorted(
            zip(
                np.round(d["x_m"] * 1e3 - centre[0], 3),
                np.round(d["y_m"] * 1e3 - centre[1], 3),
                d["slot"],
                np.round(d["axis_rad"], 6),
                strict=True,
            )
        )

    d0, c0 = run([(0.0, 0.0, 30.0, 225.0)])
    d1, c1 = run([(12.8, 6.4, 30.0, 225.0)], (12.8, 6.4))
    assert c0 == c1  # same relative cell map: the local incidence angle is taken from each cell's own feed
    assert d0["eta_conv"][0] == pytest.approx(d1["eta_conv"][0], rel=1e-3)
    d2, c2 = run(
        [(0.0, 0.0, 30.0, 45.0)]
    )  # 180 deg rotation of beam azimuth: the map is the 180 deg rotated one (slot axes unchanged: C2)
    assert d2["eta_conv"][0] == pytest.approx(d0["eta_conv"][0], rel=0.02) and c2 != []
