"""Phase 29: STEP files of the two-feed HFSS check boards (MODEL_ONLY), from the SAME geometry data as the AEDT scripts (CELLS / FEEDS / board size parsed from `build_board*.py`):
`twofeed_board.step` (2 feeds) and `control_board.step` (1 feed), mm. Named bodies (CadQuery assembly): Substrate (solid, pin holes), Ground (zero-thickness sheet at z = 0 with a coax opening of radius r_out at each
feed), Metal_pattern (all slotted circles, zero-thickness sheets at z = 1.575, one compound), Pin_k (solid cylinders r 0.635 mm, z 0 .. 7.575 mm), Port_k (annulus sheet r 0.635 .. 2.05 mm at z = 0).
Zero-thickness sheets are the geometry the AEDT scripts build (Perfect E sheets); no air box. The STEP files stay gitignored (`hfss_exports/**/*.step`); copies go to C:\\hmsa_sync\\data\\phase29_twofeed\\.
Also writes the manual STEP-import setup into the package README (between markers)."""

import ast
import re
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "hfss_exports/phase29_twofeed"
SYNC = Path(r"C:\hmsa_sync\data\phase29_twofeed")
R_PIN, R_OUT, MONO_H, H, D = 0.635, 2.05, 7.575, 1.575, 6.0


def parse(script):
    txt = (PKG / script).read_text(encoding="utf-8")
    body = re.search(r"CELLS = \[[^\n]*\n(.*?)\n\]", txt, re.S).group(1)
    cells = ast.literal_eval("[" + body + "]")
    feeds = ast.literal_eval(re.search(r"FEEDS = (\[.*?\])", txt).group(1))
    wx = float(re.search(r'WX = "([\d.]+)mm"', txt).group(1))
    wy = float(re.search(r'WY = "([\d.]+)mm"', txt).group(1))
    return cells, feeds, wx, wy


def slot_pts(x, y, length, width, ang):
    a = np.radians(ang)
    c, s = np.cos(a), np.sin(a)
    return [
        (x + u * c - v * s, y + u * s + v * c)
        for u, v in (
            (-length / 2, -width / 2),
            (length / 2, -width / 2),
            (length / 2, width / 2),
            (-length / 2, width / 2),
        )
    ]


def build(cells, feeds, wx, wy, out):
    import cadquery as cq

    V = cq.Vector
    sub = cq.Workplane("XY").box(wx, wy, H, centered=(True, True, False))
    for fx, fy in feeds:
        sub = sub.cut(cq.Workplane("XY").center(fx, fy).circle(R_PIN).extrude(H))
    ground = cq.Face.makeFromWires(
        cq.Wire.makePolygon(
            [V(-wx / 2, -wy / 2, 0), V(wx / 2, -wy / 2, 0), V(wx / 2, wy / 2, 0), V(-wx / 2, wy / 2, 0)], close=True
        ),
        [cq.Wire.makeCircle(R_OUT, V(fx, fy, 0), V(0, 0, 1)) for fx, fy in feeds],
    )
    faces = []
    for x, y, length, width, ang in cells:
        outer = cq.Wire.makeCircle(D / 2, V(x, y, H), V(0, 0, 1))
        hole = cq.Wire.makePolygon([V(px, py, H) for px, py in slot_pts(x, y, length, width, ang)], close=True)
        faces.append(cq.Face.makeFromWires(outer, [hole]))
    pattern = cq.Compound.makeCompound(faces)
    assy = cq.Assembly(name="board")
    assy.add(sub.val(), name="Substrate")
    assy.add(ground, name="Ground")
    assy.add(pattern, name="Metal_pattern")
    for k, (fx, fy) in enumerate(feeds, start=1):
        assy.add(cq.Solid.makeCylinder(R_PIN, MONO_H, V(fx, fy, 0), V(0, 0, 1)), name=f"Pin_{k}")
    for k, (fx, fy) in enumerate(feeds, start=1):
        assy.add(
            cq.Face.makeFromWires(
                cq.Wire.makeCircle(R_OUT, V(fx, fy, 0), V(0, 0, 1)),
                [cq.Wire.makeCircle(R_PIN, V(fx, fy, 0), V(0, 0, 1))],
            ),
            name=f"Port_{k}",
        )
    assy.save(str(out), exportType="STEP")


def check(path, n_cells, n_feeds):
    """Re-read the STEP with OCP (XCAF): names of the bodies and counts of solids / faces per named body."""
    from OCP.IFSelect import IFSelect_RetDone
    from OCP.STEPCAFControl import STEPCAFControl_Reader
    from OCP.TCollection import TCollection_ExtendedString
    from OCP.TDataStd import TDataStd_Name
    from OCP.TDF import TDF_LabelSequence
    from OCP.TDocStd import TDocStd_Document
    from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    from OCP.XCAFDoc import XCAFDoc_DocumentTool

    doc = TDocStd_Document(TCollection_ExtendedString("doc"))
    rd = STEPCAFControl_Reader()
    rd.SetNameMode(True)
    assert rd.ReadFile(str(path)) == IFSelect_RetDone and rd.Transfer(doc)
    st = XCAFDoc_DocumentTool.ShapeTool_s(doc.Main())
    labs = TDF_LabelSequence()
    st.GetFreeShapes(labs)
    names = {}

    def walk(label):
        comps = TDF_LabelSequence()
        if st.IsAssembly_s(label):
            st.GetComponents_s(label, comps)
            for i in range(1, comps.Length() + 1):
                walk(comps.Value(i))
            return
        ref = label
        if st.IsReference_s(label):
            from OCP.TDF import TDF_Label

            ref = TDF_Label()
            st.GetReferredShape_s(label, ref)
        nm = TDataStd_Name()
        name = None
        for lb in (label, ref):
            if lb.FindAttribute(TDataStd_Name.GetID_s(), nm):
                name = nm.Get().ToExtString()
                if name and not name.startswith("="):
                    break
        shp = st.GetShape_s(ref)
        cnt = {}
        for key, kind in (("solids", TopAbs_SOLID), ("faces", TopAbs_FACE)):
            ex = TopExp_Explorer(shp, kind)
            n = 0
            while ex.More():
                n += 1
                ex.Next()
            cnt[key] = n
        names[name] = cnt

    for i in range(1, labs.Length() + 1):
        walk(labs.Value(i))
    ok = (
        names.get("Metal_pattern", {}).get("faces") == n_cells
        and sum(1 for k in names if k and k.startswith("Pin_")) == n_feeds
        and sum(1 for k in names if k and k.startswith("Port_")) == n_feeds
        and names.get("Substrate", {}).get("solids") == 1
        and names.get("Ground", {}).get("faces") == 1
    )
    return names, ok


SECTION = """<!-- STEP-SECTION-BEGIN -->
## Manual HFSS setup from the STEP files (no AEDT script)

Files (gitignored, copies in `C:\\hmsa_sync\\data\\phase29_twofeed\\`): `twofeed_board.step` (2 feeds), `control_board.step` (1 feed). mm. Bodies (names are in the STEP product tree): **Substrate** (solid, 1.575 mm, pin holes r 0.635),
**Ground** (zero-thickness sheet at z = 0, coax opening r 2.05 mm at each feed), **Metal_pattern** (all {n2} / {n1} slotted circles, D 6.0 mm with slot, zero-thickness sheets at z = 1.575), **Pin_1, Pin_2** (solid cylinders r 0.635 mm, z 0 .. 7.575 mm;
top 6.0 mm above the substrate), **Port_1, Port_2** (annulus sheets r 0.635 .. 2.05 mm at z = 0). Zero-thickness sheets are what the AEDT scripts build (no copper thickness). Import: Modeler > Import (STEP), heal off, units mm.
Feed positions (mm): control (0, 0); two-feed port 1 ({f1:.3f}, 0), port 2 ({f2:.3f}, 0).

| Body | Assign |
|---|---|
| Substrate | new material RT/duroid 5880: eps_r 2.2, tan delta 0.0009 |
| Ground, Metal_pattern | Boundaries > Perfect E (sheets) |
| Pin_k | material pec (solve inside off) |
| Port_k | Excitations > Lumped Port, 50 ohm, integration line from (x_k + 0.635, y_k, 0) to (x_k + 2.05, y_k, 0) mm (radially across the annulus, z = 0); Port_1, Port_2 named so |
| Air box (create it) | x, y: board + 15 mm each side, z from -15 mm to 1.575 + 30 mm; Boundaries > Radiation on the box (all faces) |

Solution setup: Driven Modal, 10.54 GHz single frequency, max delta S 0.02, <= 15 passes, lambda refinement on (an interpolating sweep 10.0-11.1 GHz is optional). Far field: Radiation > Insert Far Field Setup > Infinite Sphere, theta 0-180 deg,
phi -180..180 deg, step 1 deg. One solve serves every excitation: set the sources afterwards (Results > Edit Sources: power in W, phase in deg; port 1 alone, port 2 alone, sum = 0.5 W / 0 deg each, difference = 0.5 W / 0 deg and 0.5 W / 180 deg).

Export for each board: the S-matrix (S11, S22, S21; control: S11), and for each excitation the **Realized Gain RHCP and LHCP** (and Axial Ratio) as a function of theta, phi (CSV, 1 deg): control board: its single port; two-feed board: port 1 alone, port 2 alone, sum, difference.
Compare with the prediction tables above (`predictions.json` old model, `predictions_new_model.json` angle-dependent model): the quantity of interest is the SLL (highest RHCP lobe outside 1.5 HPBW of the peak) of the two-feed ports versus the control.
<!-- STEP-SECTION-END -->
"""


def main():
    SYNC.mkdir(parents=True, exist_ok=True)
    res = {}
    for script, fname in (
        ("build_board.py", "twofeed_board.step"),
        ("build_board_control_1feed.py", "control_board.step"),
    ):
        cells, feeds, wx, wy = parse(script)
        out = PKG / fname
        build(cells, feeds, wx, wy, out)
        names, ok = check(out, len(cells), len(feeds))
        shutil.copy2(out, SYNC / fname)
        res[fname] = {"cells": len(cells), "feeds": feeds, "bytes": out.stat().st_size, "names": names, "check_ok": ok}
        print(
            fname,
            f"{out.stat().st_size / 1e6:.2f} MB",
            "check_ok",
            ok,
            names if len(names) < 12 else {k: names[k] for k in sorted(names)},
        )
    c2, f2 = parse("build_board.py")[:2]
    c1 = parse("build_board_control_1feed.py")[0]
    sec = SECTION.format(n2=len(c2), n1=len(c1), f1=f2[0][0], f2=f2[1][0])
    rd = (PKG / "README.md").read_text(encoding="utf-8")
    if "<!-- STEP-SECTION-BEGIN -->" in rd:
        rd = re.sub(r"<!-- STEP-SECTION-BEGIN -->.*?<!-- STEP-SECTION-END -->\n", lambda m: sec, rd, flags=re.S)
    else:
        rd = rd.rstrip("\n") + "\n\n" + sec
    (PKG / "README.md").write_text(rd, encoding="utf-8")
    (SYNC / "README_STEP.txt").write_text(
        "phase29_twofeed STEP files (MODEL_ONLY). Bodies: Substrate, Ground (sheet z=0), Metal_pattern (sheets z=1.575), Pin_k (r 0.635, z 0..7.575), Port_k (annulus sheet z=0). See hfss_exports/phase29_twofeed/README.md in the repo (branch phase29-seven-hex).\n",
        encoding="utf-8",
    )
    return 0 if all(v["check_ok"] for v in res.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
