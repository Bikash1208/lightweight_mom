"""KiCad-python: import the Freerouting session into the full board, fill zones, save, dump tracks.
    C:/kicad_tools/kicad/bin/python.exe sub6_kicad_post.py <kicad_dir>
"""
import json
import sys

import pcbnew

d = sys.argv[1]
b = pcbnew.LoadBoard(d + "/sub6_board.kicad_pcb")
n0 = len(b.GetTracks())
ok = pcbnew.ImportSpecctraSES(b, d + "/sub6_board_route.ses")
print("SES import", ok, "tracks before", n0, "after", len(b.GetTracks()), flush=True)
import os
stub = d + "/stub.json"
if os.path.exists(stub):   # locked RETURN stub (pre-routed in the routing board, not part of the SES): add it to the full board
    for it in json.load(open(stub)):
        t = pcbnew.PCB_TRACK(b)
        t.SetStart(pcbnew.VECTOR2I(pcbnew.FromMM(it["a"][0]), pcbnew.FromMM(-it["a"][1])))
        t.SetEnd(pcbnew.VECTOR2I(pcbnew.FromMM(it["b"][0]), pcbnew.FromMM(-it["b"][1])))
        t.SetWidth(pcbnew.FromMM(it["w"]))
        t.SetLayer(b.GetLayerID(it["layer"]))
        t.SetNet(b.FindNet(it["net"]))
        b.Add(t)
    print("stub segments added:", len(json.load(open(stub))), flush=True)
pcbnew.ZONE_FILLER(b).Fill(b.Zones())
b.Save(d + "/sub6_board_routed.kicad_pcb")
mm = pcbnew.ToMM
tr = []
nvia_new = 0
for t in b.GetTracks():
    if t.GetClass() == "PCB_VIA":
        continue
    tr.append(dict(net=t.GetNetname(), layer=t.GetLayerName(), w=mm(t.GetWidth()), a=[mm(t.GetStart().x), -mm(t.GetStart().y)], b=[mm(t.GetEnd().x), -mm(t.GetEnd().y)]))
vias = [t for t in b.GetTracks() if t.GetClass() == "PCB_VIA"]
print("segments", len(tr), "vias", len(vias), "layers", sorted({x["layer"] for x in tr}), flush=True)
json.dump(tr, open(d + "/tracks.json", "w"))
