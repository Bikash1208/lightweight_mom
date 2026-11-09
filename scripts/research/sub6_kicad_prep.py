"""KiCad-python step 1: load board, fill zones, export Specctra DSN.  Run with C:/kicad_tools/kicad/bin/python.exe"""
import sys, pcbnew
d = sys.argv[1]
b = pcbnew.LoadBoard(d + "/sub6_board.kicad_pcb")
print("loaded: footprints", len(b.GetFootprints()), "tracks", len(b.GetTracks()), "zones", len(b.Zones()), "nets", b.GetNetCount(), flush=True)
filler = pcbnew.ZONE_FILLER(b)
filler.Fill(b.Zones())
for z in b.Zones():
    print("zone", z.GetZoneName(), z.GetLayerName(), "filled", z.IsFilled(), flush=True)
b.Save(d + "/sub6_board_filled.kicad_pcb")
ok = pcbnew.ExportSpecctraDSN(b, d + "/sub6_board.dsn")
print("DSN export", ok, flush=True)
