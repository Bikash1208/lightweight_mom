"""KiCad-python: import a Freerouting session into the routing-variant board and dump all B.Cu segments (locked previous ones included) to json.
    C:/kicad_tools/kicad/bin/python.exe sub6_kicad_stage_post.py <kicad_dir> <out.json>
"""
import json
import sys

import pcbnew

d, outj = sys.argv[1], sys.argv[2]
b = pcbnew.LoadBoard(d + "/sub6_board_route.kicad_pcb")
n0 = len(b.GetTracks())
print("SES import", pcbnew.ImportSpecctraSES(b, d + "/sub6_board_route.ses"), "tracks", n0, "->", len(b.GetTracks()), flush=True)
mm = pcbnew.ToMM
tr = []
for t in b.GetTracks():
    if t.GetClass() == "PCB_VIA":
        continue
    tr.append(dict(net=t.GetNetname(), layer=t.GetLayerName(), w=mm(t.GetWidth()), a=[mm(t.GetStart().x), -mm(t.GetStart().y)], b=[mm(t.GetEnd().x), -mm(t.GetEnd().y)]))
json.dump(tr, open(outj, "w"))
print("segments", len(tr), "vias", sum(1 for t in b.GetTracks() if t.GetClass() == "PCB_VIA"), flush=True)
