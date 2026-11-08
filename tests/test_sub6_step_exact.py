"""sub6 package-layer STEP clean-up: dimensions exact to 1e-9 mm after export (4.000 mm box, 1.00 x 0.35 mm sheet, the CELL1 notch cell)."""
import re
import sys
from pathlib import Path

import pytest

cq = pytest.importorskip("cadquery")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "research"))
import sub6_step_exact as X  # noqa: E402


def _points3(text):
    t = re.sub(r"\n\s+", "", text)
    out = []
    for m in re.finditer(r"CARTESIAN_POINT\('[^']*',\(([^()]*)\)\)", t):
        v = [float(s) for s in m.group(1).split(",")]
        if len(v) == 3:
            out.append(v)
    return out


def _export(parts, path):
    asm = cq.Assembly()
    for nm, shp in parts.items():
        asm.add(shp, name=nm, color=cq.Color(0.25, 0.55, 0.25, 1))
    asm.export(str(path), exportType="STEP")
    return Path(path).read_text(encoding="utf-8")


def test_box_4mm_and_sheet_exact(tmp_path):
    box = cq.Workplane("XY").box(4.0, 4.0, 4.0, centered=False).val()
    pts = [cq.Vector(0.0, 0.0, 0.0), cq.Vector(1.0, 0.0, 0.0), cq.Vector(1.0, 0.35, 0.0), cq.Vector(0.0, 0.35, 0.0)]
    sheet = cq.Face.makeFromWires(cq.Wire.makePolygon(pts, close=True))
    raw = _export({"BOX": box, "SHEET": sheet}, tmp_path / "t.step")
    clean, stats = X.exactify_step_text(raw)
    assert X.audit_step_text(clean) == []
    allowed = (0.0, 0.175, 0.35, 0.5, 1.0, 4.0)  # corners and face-centre placement points of the box and the sheet
    for p in _points3(clean):
        for c in p:
            assert min(abs(c - g) for g in allowed) < 1e-9, c
            assert c == round(c * 1000) / 1000  # exactly the double of a 3-decimal literal
    assert "3.9999999" not in clean and "0.3500000" not in clean
    assert stats["n_colours"] >= 1
    cols = re.findall(r"COLOUR_RGB\('',([^)]*)\)", re.sub(r"\n\s+", "", clean))
    assert all(abs(float(x) - round(float(x), 4)) == 0 for c in cols for x in c.split(","))


def test_noisy_literals_are_snapped_irrationals_kept():
    txt = (
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('x 7.9'),'2;1');\nENDSEC;\nDATA;\n"
        "#1 = CARTESIAN_POINT('',(3.9999999999,0.3500000001,0.));\n"
        "#2 = CARTESIAN_POINT('',(5.411952050785,-0.0350006));\n"
        "#3 = COLOUR_RGB('',0.250000000672,0.55000000909,0.250000000672);\n"
        "#4 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#5,'distance_accuracy_value','x');\n"
        "ENDSEC;\nEND-ISO-10303-21;\n"
    )
    assert len(X.audit_step_text(txt)) >= 4
    clean, _ = X.exactify_step_text(txt)
    assert X.audit_step_text(clean) == []
    lit = X.literals(clean)
    assert abs(lit[0] - 4.0) < 1e-9 and abs(lit[1] - 0.35) < 1e-9
    assert abs(lit[3] - 5.411952050785) < 1e-12 and abs(lit[4] - (-0.035)) < 1e-9
    assert "1.E-6" in clean  # uncertainty untouched


def test_sub6_cell1_notch_has_no_inexact_literals(tmp_path):
    import sub6_slotdisc_hfss_package as P

    parts = P.make("CELL1", True)
    raw = _export(parts, tmp_path / "c1.step")
    clean, stats = X.exactify_step_text(raw)
    assert X.audit_step_text(clean) == []
    # thickness 0.035 and the 1.00 mm notch edge are exactly on the grid
    vals = {round(c, 9) for p in _points3(clean) for c in p}
    assert 0.035 in vals and 0.5 in vals and -0.5 in vals
