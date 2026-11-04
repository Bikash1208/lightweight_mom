"""Side quest: mushroom-type gate cell estimate (scripts/research/sidequest_mushroom_cell.py) — model sanity, no MoM run (uses the committed sheet table)."""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
warnings.simplefilter("ignore")
import sidequest_mushroom_cell as mc  # noqa: E402

RAW = json.loads((ROOT / "results/sidequest_hybrid/mushroom_theory_raw.json").read_text(encoding="utf-8"))["5.6"]


def test_lc_formula_is_sievenpiper():
    e = mc.lc_estimate(5.6e-3, 0.8e-3)
    assert abs(e["L_nH"] - 4e-7 * np.pi * 1.575e-3 * 1e9) < 1e-9
    assert abs(e["f0_GHz"] - 1 / (2 * np.pi * np.sqrt(e["L_nH"] * 1e-9 * e["C_pF"] * 1e-12)) / 1e9) < 1e-9


def test_on_surface_band_rises_with_phase_and_stays_below_f0():
    m = mc.Model(RAW["table"])
    f0 = m.f0("ON")[0]
    fs = [m.f_of_phase("ON", ph, None)[0] for ph in (120, 150, 180)]
    assert fs[0] < fs[1] < fs[2] < f0


def test_gap_covers_the_target_band_and_off_beta_is_near_the_hologram_average():
    pred = json.loads((ROOT / "hfss_exports/sidequest_mushroom_cell/predictions.json").read_text(encoding="utf-8"))
    g = pred["ON_stopband_GHz"]
    assert g["lower_edge_phase180_TM"] <= 10.24 and g["upper_edge_via_patch_resonance"] >= 10.84
    assert abs(pred["OFF_beta_over_k0_at_10p54"] - 1.065) < 0.03


def test_step_variants_have_the_named_parts():
    pytest.importorskip("cadquery")
    from src import cad_export as ce

    for st in ("ON", "OFF"):
        parts = ce.import_named_parts(str(ROOT / f"hfss_exports/sidequest_mushroom_cell/cell_{st}.step"))
        assert {"Substrate", "Ground", "Via", "Pad", "Patch", "Air", f"D_{st}"} <= set(parts)
        assert parts["Air"]["bbox_mm"][5] == pytest.approx(30.0) and parts["Via"]["n_solids"] == 1


def test_variant_predictions_are_ordered_and_files_exist():
    d = ROOT / "hfss_exports/sidequest_mushroom_variants"
    pred = json.loads((d / "predictions.json").read_text(encoding="utf-8"))
    f180 = {k: v["S0"]["M2"]["ON_180"] for k, v in pred["variants"].items()}
    assert f180["V1"] > f180["V2"] > f180["V3"]  # smaller patch -> higher resonance (calibrated M2)
    assert all(v["S0"]["M1"]["ON_180"] < v["S0"]["M2"]["ON_180"] for v in pred["variants"].values())  # M1 is the low-side alternative
    assert sorted(pred["files"]) == sorted(["V1_ON.step", "V1_OFF.step", "V2_ON.step", "V2_OFF.step", "V3_ON.step", "V4_ON.step"])
    for f in pred["files"]:
        assert (d / f).exists()


def test_variant_step_via_radius_and_patch_volume():
    pytest.importorskip("cadquery")
    from src import cad_export as ce

    d = ROOT / "hfss_exports/sidequest_mushroom_variants"
    p = ce.import_named_parts(str(d / "V4_ON.step"))
    assert p["Via"]["volume_mm3"] == pytest.approx(np.pi * 0.5**2 * 1.575, rel=1e-3)
    assert p["Patch"]["volume_mm3"] == pytest.approx((2.8**2 - 2.4**2) * 0.035, rel=1e-3)  # hole = pad 0.8 + 2 x 0.8 gap = 2.4
