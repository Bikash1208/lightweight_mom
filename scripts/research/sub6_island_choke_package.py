"""SUB6 island cell: HFSS eigenmode package with bias vias, variant 2 (vias, no choke) and variant 3 (vias + analytic chokes), FR4 thickness from the choke table.  MODEL_ONLY.

Frame: the island package (slot along x, vias on x = 0 at y = -0.20 island / +2.35 disc) ROTATED by +90 deg about z, so the slot runs along y and ACROSS the slot = x = $ph_x
(the user's model: ACROSS = $ph_x, ALONG = $ph_y).  Stack (z, mm): patch 0..0.035, RT5880 -3.175..0, ground -3.210..-3.175, FR4 (t) below, bottom copper 35 um, AIR 30 mm below the
bottom copper with a RADIATION boundary on its bottom face, PEC lid on the top face of AIR (z = 150).  Lattice pairs on the four AIR side faces.

    python sub6_island_choke_package.py geometry [t_fr4_mm]     -> hfss_exports/sub6_island_choke/{ISL_vias2.step, ISL_choke3.step, geometry_checks.json}
"""
import json
import math
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sub6_island_choke as ch
import sub6_island_choke_layout as lay
import sub6_island_hfss_package as base

OUT = ROOT / "hfss_exports" / "sub6_island_choke"
RES = ROOT / "results" / "sub6_island_choke"
P, H, T, GT, LID = base.P, base.H, base.T, base.GT, base.LID
AIR_BELOW = 30.0
GRID = 0.001


def z_levels(t_fr4):
    z_gb = -H - GT
    z_fr4b = round(z_gb - t_fr4, 6)
    z_bot = round(z_fr4b - GT, 6)
    return dict(z_gb=z_gb, z_fr4b=z_fr4b, z_bot=z_bot, z_air=round(z_bot - AIR_BELOW, 6))


def g(v):
    return round(round(v / GRID) * GRID, 6)


def rot90(x, y):
    return (-y, x)  # +90 deg about z: old slot (x) -> new y; old across (y) -> new x


def choke_layout(t_fr4):
    """-> (design row, [choke dicts rebuilt with the table's parameters]) for the FR4 thickness."""
    rows = json.loads((RES / "choke_table.json").read_text())
    row = [r for r in rows if abs(r["h_mm"] - t_fr4) < 1e-9 and r["layout"]]
    row = sorted(row, key=lambda r: -r["zin_min_ohm"])[0]
    chokes = [lay.build(c["via"], c["dir"], c["flip"], row["w_mm"], row["spacing_mm"], row["l_mm"], row["r2_mm"], row["angle_deg"], c["n_runs"]) for c in row["layout"]["chokes"]]
    return row, chokes


def make(variant, t_fr4, chokes=None):
    """variant: 'vias' (2) | 'choke' (3)  -> {name: cq shape} (already rotated into the ACROSS = x frame)."""
    import cadquery as cq

    zl = z_levels(t_fr4)
    z_gb, z_fr4b, z_bot, z_air = zl["z_gb"], zl["z_fr4b"], zl["z_bot"], zl["z_air"]
    vr, pr, ar = base.VIA_R, base.PAD_R, base.ANTIPAD_R

    def box(x0, x1, y0, y1, za, zb):
        return cq.Workplane("XY").workplane(offset=za).center((x0 + x1) / 2, (y0 + y1) / 2).rect(x1 - x0, y1 - y0).extrude(zb - za)

    def cyl(x, y, r, za, zb):
        return cq.Workplane("XY").workplane(offset=za).center(x, y).circle(r).extrude(zb - za)

    patch = cq.Workplane("XY").circle(base.R).extrude(T).cut(box(-base.L / 2, base.L / 2, -base.W / 2, base.W / 2, -0.01, T + 0.01))
    island = box(-base.X_ISL, base.X_ISL, base.Y_ISL0, base.Y_ISL1, 0.0, T)
    ground = box(-P / 2, P / 2, -P / 2, P / 2, z_gb, -H)
    yi, yd = base.Y_VIA_I, base.Y_VIA_D

    def barrel_pad(y):
        return cyl(0.0, y, vr, z_bot, 0.0).union(cyl(0.0, y, pr, z_bot, z_fr4b))

    isl = island.union(barrel_pad(yi))
    dsc = patch.union(barrel_pad(yd))
    extra = {}
    if variant == "vias":  # the previous variant 2: 2.0 x 0.3 bias stubs (island towards +x, disc towards -x), no choke
        isl = isl.union(box(0.0, base.STUB_L, yi - base.STUB_W / 2, yi + base.STUB_W / 2, z_bot, z_fr4b))
        dsc = dsc.union(box(-base.STUB_L, 0.0, yd - base.STUB_W / 2, yd + base.STUB_W / 2, z_bot, z_fr4b))
    else:
        for name, c in (("CHOKE_ISLAND", chokes[0]), ("CHOKE_DISC", chokes[1])):
            solids = []
            w = c["w"]
            for k, poly in enumerate(c["parts"]):
                if k == 0:
                    continue  # pad already in barrel_pad
                if k == len(c["parts"]) - 1:  # radial stub: sector with an exact-centre arc
                    apex = (g(c["apex"][0]), g(c["apex"][1]))
                    pts = [(g(x), g(y)) for x, y in poly[1:]]
                    mid = pts[len(pts) // 2]
                    wp = cq.Workplane("XY").workplane(offset=z_bot).moveTo(*apex).lineTo(*pts[0]).threePointArc(mid, pts[-1]).close()
                    solids.append(wp.extrude(z_fr4b - z_bot).val())
                else:
                    xs = [g(x) for x, _ in poly]
                    ys = [g(y) for _, y in poly]
                    solids.append(box(min(xs), max(xs), min(ys), max(ys), z_bot, z_fr4b).val())
            extra[name] = cq.Compound.makeCompound(solids)
    ground = ground.cut(cyl(0.0, yi, ar, z_gb - 0.01, -H + 0.01)).cut(cyl(0.0, yd, ar, z_gb - 0.01, -H + 0.01))
    metal_solids = [dsc.val(), isl.val(), ground.val()]
    sub = box(-P / 2, P / 2, -P / 2, P / 2, -H, 0.0)
    air = box(-P / 2, P / 2, -P / 2, P / 2, z_air, LID)
    fr4 = box(-P / 2, P / 2, -P / 2, P / 2, z_fr4b, z_gb)
    pts = [cq.Vector(-base.ALONG / 2, base.Y_SH0, base.Z_SH), cq.Vector(base.ALONG / 2, base.Y_SH0, base.Z_SH), cq.Vector(base.ALONG / 2, base.Y_SH1, base.Z_SH), cq.Vector(-base.ALONG / 2, base.Y_SH1, base.Z_SH)]
    parts = {"METAL_PEC": cq.Compound.makeCompound(metal_solids), "SUBSTRATE": sub.val(), "SUBSTRATE_FR4": fr4.val(), "AIR": air.val(),
             "DIODE_SHEET": cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))}
    parts.update(extra)
    return {k: v.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 90.0) for k, v in parts.items()}


def intended(variant, t_fr4, chokes):
    zl = z_levels(t_fr4)
    xy = {0.0, P / 2, base.R, base.L / 2, base.W / 2, base.GC / 2, base.X_ISL, base.ALONG / 2, base.Y_ISL0, base.Y_ISL1, base.Y_SH0, base.Y_SH1, base.STUB_L,
          base.Y_VIA_I, base.Y_VIA_D, base.VIA_R, base.PAD_R, base.ANTIPAD_R, base.STUB_L + base.PAD_R}
    for yv in (base.Y_VIA_I, base.Y_VIA_D):
        xy |= {yv + base.PAD_R, yv - base.PAD_R, yv + base.STUB_W / 2, yv - base.STUB_W / 2}
    if chokes:
        for c in chokes:
            for poly in c["parts"][1:]:
                xy |= {g(v) for pt in poly for v in pt}
            xy |= {g(c["apex"][0]), g(c["apex"][1])}
    xy |= {-v for v in xy}
    z = {0.0, T, T / 2, -H, LID, -0.01, T + 0.01, zl["z_gb"], zl["z_fr4b"], zl["z_bot"], zl["z_air"], zl["z_gb"] - 0.01, -H + 0.01}
    z |= {-v for v in z}
    return dict(xy=np.array(sorted(round(v, 6) for v in xy)), z=np.array(sorted(round(v, 6) for v in z)), dir=np.array([-1.0, 0.0, 1.0]))


def geometry(t_fr4):
    import cadquery as cq
    from src.transmitarray.export import inexact_literals, snap_step_text, verify_snap
    from sub6_step_exact import audit_step_text, exactify_step_text

    OUT.mkdir(parents=True, exist_ok=True)
    row, chokes = choke_layout(t_fr4)
    cols = {"METAL_PEC": (0.9, 0.7, 0.2), "CHOKE_ISLAND": (0.9, 0.5, 0.2), "CHOKE_DISC": (0.9, 0.5, 0.2), "SUBSTRATE": (0.25, 0.55, 0.25), "SUBSTRATE_FR4": (0.2, 0.5, 0.2), "AIR": (0.7, 0.8, 1.0), "DIODE_SHEET": (0.9, 0.1, 0.1)}
    out = {}
    for variant, fn in (("vias", "ISL_vias2.step"), ("choke", "ISL_choke3.step")):
        parts = make(variant, t_fr4, chokes if variant == "choke" else None)
        asm = cq.Assembly()
        for nm, shp in parts.items():
            asm.add(shp, name=nm, color=cq.Color(*cols[nm], 1))
        path = OUT / fn
        asm.export(str(path), exportType="STEP")
        raw = path.read_text(encoding="utf-8")
        it = intended(variant, t_fr4, chokes if variant == "choke" else None)
        # rotated frame: the table is sign-symmetric and uses the same set for x and y, so it covers the rotation
        snapped, stats = snap_step_text(raw, it, description=f"sub6 island cell {variant}, FR4 {t_fr4} mm, mm, 0.001 mm grid")
        snapped, ex = exactify_step_text(snapped)
        path.write_text(snapped, encoding="utf-8", newline="\n")
        text = path.read_text(encoding="utf-8")
        bad = list(inexact_literals(text, it)) + audit_step_text(text)
        out[fn] = dict(exact_ok=len(bad) == 0, inexact_literals=len(bad), n_snapped=stats["n_snapped"], first_bad=[str(b) for b in bad[:5]])
        print(variant, out[fn], flush=True)
    return out


if __name__ == "__main__":
    if sys.argv[1] == "geometry":
        geometry(float(sys.argv[2]) if len(sys.argv) > 2 else 0.5)
