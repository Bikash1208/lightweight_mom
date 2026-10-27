"""Phase 28 (cell redesign): what the substrate does to the loading a cell needs. For each board option: the mean sheet susceptance b_bar that puts the TM0 mode at beta/k0 = 1.053 (the item-2 mean), the bare-slab beta/k0 (b_bar = 0), and the leakage alpha the cell G's sheet tensor
(b_bar 4.46 mS, delta 2.29 mS, 10.0 GHz, extracted on RT/duroid 5880 1.575 mm) WOULD give on that board IF the sheet tensor were unchanged (an approximation: the tensor depends on the substrate; it must be re-extracted).
MODEL_ONLY; options listed, none chosen. python scripts/research/phase28_substrate_options.py"""
import json, math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scipy.optimize import brentq
from src.dispersion.leakage_model import leakage_alpha_first_order_tensor
from src.dispersion.sheet_model import find_beta_from_sheet_admittance
from src.frequency import free_space_wavenumber
from src.materials import Substrate

F = 10.0e9
K0 = free_space_wavenumber(F)
B_G, D_G = 4.457e-3, 2.293e-3  # G, minimum over states, mesh 96
OPTIONS = [("RT/duroid 5880", 2.2, 0.0009, 1.575e-3), ("RT/duroid 5880", 2.2, 0.0009, 0.787e-3), ("RT/duroid 5880", 2.2, 0.0009, 3.175e-3),
           ("RO4003C", 3.38, 0.0027, 0.813e-3), ("RO4003C", 3.38, 0.0027, 1.524e-3), ("RO3010", 10.2, 0.0022, 0.635e-3), ("RO3010", 10.2, 0.0022, 1.28e-3)]
rows = []
for name, er, td, h in OPTIONS:
    sub = Substrate(name, er, td, h)
    row = {"substrate": name, "eps_r": er, "tan_delta": td, "thickness_mm": h * 1e3}
    try:
        row["beta_over_k0_bare_slab"] = find_beta_from_sheet_admittance(1e-9j, sub, F).beta_rad_per_m / K0
    except Exception as e:
        row["beta_over_k0_bare_slab"] = None
    def f(b):
        return find_beta_from_sheet_admittance(1j * b * 1e-3, sub, F).beta_rad_per_m / K0 - 1.053
    try:
        row["b_bar_for_beta_1p053_ms"] = brentq(f, 0.05, 40.0, xtol=1e-9)
    except Exception:
        row["b_bar_for_beta_1p053_ms"] = None
    try:
        beta = find_beta_from_sheet_admittance(1j * B_G, sub, F).beta_rad_per_m
        row["G_beta_over_k0"] = beta / K0
        row["G_alpha_np_per_m_if_tensor_unchanged"] = float(leakage_alpha_first_order_tensor(beta, 1j * B_G, D_G, math.radians(30), sub, F))
    except Exception as e:
        row["G_beta_over_k0"] = None; row["G_alpha_np_per_m_if_tensor_unchanged"] = None
    rows.append(row)
out = {"labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"], "frequency_ghz": 10.0, "rows": rows}
Path(ROOT / "results/phase28/f10/two_diode/substrate_options.json").write_text(json.dumps(out, indent=1))
for r in rows: print(r)
