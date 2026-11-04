"""HFSS package: ONE mushroom-type gate unit cell (patch + via-diode), two STEP variants (ON: sheet D_ON, OFF: sheet D_OFF), eigenmode with a PEC lid. MODEL_ONLY predictions.
    python scripts/research/sidequest_mushroom_package.py  -> hfss_exports/sidequest_mushroom_cell/"""
import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import cadquery as cq
import numpy as np

import sidequest_mushroom_cell as mc
from src.frequency import free_space_wavenumber

OUT = ROOT / "hfss_exports/sidequest_mushroom_cell"
P, H, T = 6.4, 1.575, 0.035  # mm
W, PAD, GD, R_VIA, LID = 5.6, 0.8, 0.4, 0.2, 30.0  # patch side, pad side, diode gap, via radius, PEC lid height above z = 0 [mm]
HH = PAD / 2 + GD


def box(cx, cy, lx, ly, z0, z1):
    return cq.Workplane("XY").workplane(offset=z0).center(cx, cy).rect(lx, ly).extrude(z1 - z0)


def build(state):
    ring = [box(-(W / 2 + HH) / 2, 0, W / 2 - HH, W, 0, T), box((W / 2 + HH) / 2, 0, W / 2 - HH, W, 0, T),
            box(0, (W / 2 + HH) / 2, 2 * HH, W / 2 - HH, 0, T), box(0, -(W / 2 + HH) / 2, 2 * HH, W / 2 - HH, 0, T)]
    patch = ring[0]
    for r in ring[1:]:
        patch = patch.union(r)
    pad = box(0, 0, PAD, PAD, 0, T)
    via = cq.Workplane("XY").workplane(offset=-H).circle(R_VIA).extrude(H)
    sub = box(0, 0, P, P, -H, 0).cut(via)
    gnd = box(0, 0, P, P, -H - T, -H)
    air = box(0, 0, P, P, 0, LID).cut(patch).cut(pad)
    zs = T / 2
    x0, x1, y0, y1 = PAD / 2, HH, -PAD / 2, PAD / 2  # diode gap between pad edge and patch hole edge, current along +x
    sheet = cq.Workplane("XY").add(cq.Face.makeFromWires(cq.Wire.makePolygon([cq.Vector(x0, y0, zs), cq.Vector(x1, y0, zs), cq.Vector(x1, y1, zs), cq.Vector(x0, y1, zs)], close=True)))
    asm = cq.Assembly()
    for nm, o, c in (("Substrate", sub, (0.25, 0.55, 0.25)), ("Ground", gnd, (0.6, 0.6, 0.65)), ("Via", via, (1, 1, 0)), ("Pad", pad, (0.8, 0.55, 0.1)), ("Patch", patch, (0.9, 0.7, 0.2)), ("Air", air, (0.7, 0.8, 1.0))):
        asm.add(o, name=nm, color=cq.Color(*c, 1))
    asm.add(sheet, name="D_ON" if state == "ON" else "D_OFF", color=cq.Color(0.9, 0.1, 0.1, 1))
    asm.export(str(OUT / f"cell_{state}.step"), exportType="STEP")


def predictions():
    raw = json.loads((ROOT / "results/sidequest_hybrid/mushroom_theory_raw.json").read_text(encoding="utf-8"))["5.6"]
    m = mc.Model(raw["table"])
    lid = LID * 1e-3
    k0 = lambda f: free_space_wavenumber(f)  # noqa: E731
    table = []
    for ph in (0, 30, 60, 90, 120, 150, 180):
        ll = mc.C0 * np.radians(ph) / mc.P / (2 * np.pi) / 1e9 if ph else 0.0
        row = {"phase_deg": ph, "light_line_GHz": round(ll, 3)}
        for st in ("ON", "OFF"):
            fs = [x / 1e9 for x in m.f_of_phase(st, ph, lid, fmin=0.3e9, fmax=16e9, n=1500)]
            below = [round(x, 3) for x in fs if (ph == 0 or x < ll - 1e-3)]  # surface-wave (below the light line) crossings
            row[f"{st}_TM_surface_band_GHz"] = (below[-1] if below else None) if ph else None  # phase 0: the surface band starts at 0 Hz  # the highest crossing below the light line = the surface-wave band
            row[f"{st}_all_model_crossings_GHz"] = [round(x, 3) for x in fs][:5]
        table.append(row)
    f0 = m.f0("ON")[0] / 1e9
    f180 = m.f_of_phase("ON", 180, None)[0] / 1e9
    beta_off = [b / k0(10.54e9) for b in m.beta_at("OFF", 10.54e9, None) if b / k0(10.54e9) < 1.4][0]
    return {"label": "MODEL_ONLY (lumped wire-medium via model + MoM patch sheet; uncertainty ~ +-5 % in frequency, the HFSS eigenmode decides)",
            "cell_mm": {"lattice": P, "patch": W, "pad": PAD, "diode_gap": GD, "via_radius": R_VIA, "substrate": H, "lid_height_above_patch_plane": LID},
            "lc_theory": raw["lc"] | {"note": "Sievenpiper LC with w = 5.6, g = 0.8: C overestimates the MoM sheet capacitance (0.0975 pF) by 1.76x, f0 8.6 GHz vs 11.0 GHz; the design uses the MoM sheet"},
            "diode_values_10p54": {"D_ON_R_ohm": 2.44, "D_OFF_C_fF": 28.16},
            "ON_stopband_GHz": {"lower_edge_phase180_TM": round(f180, 3), "upper_edge_via_patch_resonance": round(f0, 3), "covers_10p54_pm_0p3": bool(f180 <= 10.24 and f0 >= 10.84)},
            "OFF_beta_over_k0_at_10p54": round(float(beta_off), 4), "hologram_average_beta_over_k0": 1.065,
            "eigenmode_table_with_lid": table,
            "reading_notes": ["TM surface band = the crossing just below the light line (f = c*phase/(2 pi a)); crossings at/above the light line are lid (parallel-plate) modes.",
                              "ON is lossy (R = 2.44 ohm): HFSS reports complex frequencies; compare Re(f).",
                              "OFF: the TM surface band leaves the 16 GHz window above ~120 deg (model), 'None' is not a failure.",
                              "Modes whose field is concentrated in the substrate between patch/pad and ground (large Ez at the via) are the surface-wave band; TE-like bands (Ey) are not in the table."]}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for st in ("ON", "OFF"):
        build(st)
    pred = predictions()
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1, default=float), encoding="utf-8", newline="\n")
    print(json.dumps({k: pred[k] for k in ("ON_stopband_GHz", "OFF_beta_over_k0_at_10p54")}, default=float))
    for r in pred["eigenmode_table_with_lid"]:
        print(r["phase_deg"], r["light_line_GHz"], r["ON_TM_surface_band_GHz"], r["OFF_TM_surface_band_GHz"])
