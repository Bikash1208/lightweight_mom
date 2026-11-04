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
