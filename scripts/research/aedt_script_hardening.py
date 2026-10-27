"""Phase 28: shared hardening of the generated AEDT (IronPython) scripts after the first AEDT run of `phase28_loop_R4` (2026-10-27).

Two failures were found there and are fixed for every generated package:
1. `AssignMaster` accepts exactly ONE face. A side wall of the cell is two faces (substrate and air box), so the periodic walls are assigned per face: four master/slave pairs (substrate X, air X, substrate Y, air Y),
   each with its own u axis and an origin in the plane of ITS face (z = 0 for the substrate faces, z = H for the air faces); every `Faces:=` list holds one face id.
2. The fallback messages used a design name that was copied from another package. `note()` uses `oDesign.GetName()` and can never raise (a failed message is swallowed, so a fallback can never stop the script).
`harden(script)` rewrites a generated script accordingly: every `oDesktop.AddMessage(oProject.GetName(), <design>, level, text)` becomes `note(level, text)`, `note` is defined right after the design is set,
and the periodic-wall + Floquet-port block (from `def faces_at` to the analysis setup) is replaced by `PERIODIC_FLOQUET`.
"""

import re

NOTE_DEF = '''import sys
def note(level, text):
    # print-only, never raises: a failed message must not stop the script
    try:
        oDesktop.AddMessage(oProject.GetName(), oDesign.GetName(), level, text)
    except Exception:
        pass
    try:
        sys.stdout.write(str(text) + "\\n")
    except Exception:
        pass
'''

PERIODIC_FLOQUET = '''def one_face(body, x, y, z):
    # exactly one face id: AssignMaster / AssignSlave take one face each
    return int(oEditor.GetFaceByPosition(["NAME:FaceParameters", "BodyName:=", body, "XPosition:=", mm(x), "YPosition:=", mm(y), "ZPosition:=", mm(z)]))

# periodic walls: one master/slave pair PER FACE (substrate and air box separately), scan angle linked to the Floquet port. NOT yet run successfully in AEDT: if a pair fails the message names it,
# assign that pair by hand (README step 2).
for tag, body, z0, zc in (("s", "Substrate", 0.0, H / 2.0), ("a", "Air", H, H + LID / 2.0)):
    try:
        oBnd.AssignMaster(["NAME:MasterX" + tag, "Faces:=", [one_face(body, -P / 2, 0.0, zc)],
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(-P / 2), mm(-P / 2), mm(z0)], "UPos:=", [mm(-P / 2), mm(P / 2), mm(z0)]],
            "ReverseV:=", False])
        oBnd.AssignSlave(["NAME:SlaveX" + tag, "Faces:=", [one_face(body, P / 2, 0.0, zc)],
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(P / 2), mm(-P / 2), mm(z0)], "UPos:=", [mm(P / 2), mm(P / 2), mm(z0)]],
            "ReverseU:=", False, "Master:=", "MasterX" + tag, "UseScanAngle:=", True, "Phi:=", "scan_phi", "Theta:=", "scan_theta"])
    except Exception as exc:
        note(2, "Periodic X pair '" + tag + "' (" + body + ") failed (" + str(exc) + "): assign it by hand, README step 2.")
    try:
        oBnd.AssignMaster(["NAME:MasterY" + tag, "Faces:=", [one_face(body, 0.0, -P / 2, zc)],
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(-P / 2), mm(-P / 2), mm(z0)], "UPos:=", [mm(P / 2), mm(-P / 2), mm(z0)]],
            "ReverseV:=", False])
        oBnd.AssignSlave(["NAME:SlaveY" + tag, "Faces:=", [one_face(body, 0.0, P / 2, zc)],
            "CoordSysVector:=", ["NAME:Vectors", "Origin:=", [mm(-P / 2), mm(P / 2), mm(z0)], "UPos:=", [mm(P / 2), mm(P / 2), mm(z0)]],
            "ReverseU:=", False, "Master:=", "MasterY" + tag, "UseScanAngle:=", True, "Phi:=", "scan_phi", "Theta:=", "scan_theta"])
    except Exception as exc:
        note(2, "Periodic Y pair '" + tag + "' (" + body + ") failed (" + str(exc) + "): assign it by hand, README step 2.")
try:
    top = one_face("Air", 0.0, 0.0, H + LID)
    oBnd.AssignFloquetPort(["NAME:FloquetPort1", "Faces:=", [top], "NumModes:=", 2, "RenormalizeAllTerminals:=", True,
        "UseScanAngle:=", True, "Phi:=", "scan_phi", "Theta:=", "scan_theta", "DoDeembed:=", True, "DeembedDist:=", mm(LID),
        ["NAME:Modes", ["NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", False, "CharImp:=", "Zpi"],
                       ["NAME:Mode2", "ModeNum:=", 2, "UseIntLine:=", False, "CharImp:=", "Zpi"]],
        "ShowReporterFilter:=", False, "ReporterFilter:=", [True, True]])
except Exception as exc:
    note(2, "Floquet port assignment failed (" + str(exc) + "): assign it by hand on the top face of Air, 2 modes, de-embed " + mm(LID) + ", README step 2.")
'''

_REASSIGN = re.compile(r"(?P<i> *)oBnd\.DeleteBoundaries\((?P<a>.*?)\)\n(?P=i)oBnd\.AssignLumpedRLC\(args\)\n")
_ADD = re.compile(r'oDesktop\.AddMessage\(oProject\.GetName\(\), (?:"[^"]*"|oDesign\.GetName\(\)), ')


def _reassign(m):
    i, a = m["i"], m["a"]
    return (
        f"{i}try:\n{i}    oBnd.DeleteBoundaries({a})\n{i}    oBnd.AssignLumpedRLC(args)\n"
        f"{i}except Exception as exc2:\n{i}    note(2, 'reassigning the lumped RLC failed: ' + str(exc2))\n"
    )


def harden(text):
    """Rewrite a generated AEDT script as described in the module docstring (idempotent)."""
    if "def note(level, text):" in text:
        return text
    text = _ADD.sub("note(", text)
    text = _REASSIGN.sub(_reassign, text)
    if "def faces_at(x, y):" in text:
        a = text.index("def faces_at(x, y):")
        b = text.index('oAna = oDesign.GetModule("AnalysisSetup")')
        text = text[:a] + PERIODIC_FLOQUET + text[b:]
    for anchor in ('oEditor = oDesign.SetActiveEditor("3D Modeler")\n', "oDesign = oProject.GetActiveDesign()\n"):
        if anchor in text:
            text = text.replace(anchor, anchor + NOTE_DEF, 1)
            break
    else:
        raise ValueError("no design anchor found for note()")
    return text
