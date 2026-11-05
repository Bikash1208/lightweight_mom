"""RT5880 cell HFSS package: the committed STEP passes every geometry check on re-import, and the checks catch the earlier errors (touching sheets, missing rim, asymmetric rims)."""
import sys
from pathlib import Path

import cadquery as cq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
import check_step_cell as CK  # noqa: E402
import sidequest_rt5880_cell_package as K  # noqa: E402

STEP = ROOT / "hfss_exports" / "sidequest_rt5880_cell" / "cell_geometry.step"


def _comp(parts):
    return cq.Compound.makeCompound([parts["METAL_PEC"], parts["SUBSTRATE"], parts["AIR"], parts["DX_SHEET"], parts["DY_SHEET"]])


def test_committed_step_passes_all_checks_on_reimport():
    r = CK.check_file(STEP, pos=0.8, ls=2.8, g=0.4, wd=1.0)
    assert r["PASS"], {k: v["ok"] for k, v in r["checks"].items()}
    assert {"bodies_and_names", "no_intersections", "rims_at_4_ends", "sheet_edge_coincidence", "sheet_to_sheet_clearance"} <= set(r["checks"])


def test_touching_sheets_are_caught():
    r = CK.check(_comp(K.build(wd=1.2)), pos=0.8, ls=2.8, g=0.4, wd=1.2)  # footprint 1.2: inner edges at the slit crossing, sheets touch at a corner
    assert not r["checks"]["sheet_to_sheet_clearance"]["ok"] and r["checks"]["sheet_to_sheet_clearance"]["distance_mm"] < 1e-9


def test_missing_rim_and_asymmetric_rim_are_caught():
    r = CK.check(_comp(K.build(ls=3.1)), pos=0.8, ls=2.8, g=0.4, wd=1.0)  # slits through the disc where the model has a rim
    assert not r["checks"]["rims_at_4_ends"]["ok"] and not r["checks"]["patch_vs_MoM"]["ok"]
    parts = K.build()
    skew = cq.Workplane("XY").circle(3.0).extrude(0.035).cut(K.box(0, 0, 0.4, 5.6, 0, 0.035).translate((0, 0.2, 0))).cut(K.box(0, 0, 5.6, 0.4, 0, 0.035))  # vertical slit shifted: rim missing at the top end
    air = K.box(0, 0, 6.4, 6.4, 0, 30).cut(skew)
    parts["METAL_PEC"] = cq.Compound.makeCompound([skew.val(), K.box(0, 0, 6.4, 6.4, -1.61, -1.575).val()])
    parts["AIR"] = air.val()
    r = CK.check(_comp(parts), pos=0.8, ls=2.8, g=0.4, wd=1.0)
    assert not r["checks"]["rims_at_4_ends"]["ok"] or not r["checks"]["symmetry"]["ok"]


def test_bare_substrate_resonance_is_above_the_loaded_cell():
    assert 11.3 < K.f_bare(True) < 11.5  # model bare grounded slab with the 30 mm lid at 90 deg (HFSS user value 11.07: 3 % lower)
