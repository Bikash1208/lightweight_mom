"""KiCad-python: export the routing variant to Specctra DSN.
 * only B.Cu is routable: F.Cu / In1.Cu / In2.Cu get a board-wide keep-out (this also blocks via placement everywhere);
 * Freerouting mis-handles the hexagonal outline (edges at 30/60 deg: pins inside the board become unreachable), so the DSN boundary is replaced by a large square and the
   outside of the hexagon (shrunk by EDGE_KEEP mm = copper-to-edge rule) is blocked on B.Cu by six keep-out quadrilaterals.
    C:/kicad_tools/kicad/bin/python.exe sub6_kicad_dsn.py <kicad_dir>
"""
import re
import sys

import pcbnew

EDGE_KEEP = 0.9   # mm, trace/copper edge to board outline (DRC min_copper_edge_clearance)
d = sys.argv[1]
b = pcbnew.LoadBoard(d + "/sub6_board_route.kicad_pcb")
print("loaded", len(b.GetFootprints()), len(b.GetTracks()), b.GetNetCount(), flush=True)
print("DSN", pcbnew.ExportSpecctraDSN(b, d + "/sub6_board_route.dsn"), flush=True)
t = open(d + "/sub6_board_route.dsn", encoding="utf-8").read()
m = re.search(r"\(boundary\s*\(path pcb 0\s+([^)]*)\)\s*\)", t)
coords = [float(v) for v in m.group(1).split()]
pts = list(zip(coords[0::2], coords[1::2]))
if pts[-1] == pts[0]:
    pts = pts[:-1]
ap = max(abs(x) for x, y in pts)  # apothem (um): vertices at x = +-apothem
s = (ap - EDGE_KEEP * 1000.0) / ap
hexs = [(x * s, y * s) for x, y in pts]
far = [(x * 6.0, y * 6.0) for x, y in pts]
big = 700000
poly = "%d %d  %d %d  %d %d  %d %d  %d %d" % (-big, -big, big, -big, big, big, -big, big, -big, -big)
t = t.replace(m.group(0), "(boundary\n      (path pcb 0  %s)\n    )" % poly)
quads = ""
for i in range(6):
    a, c = hexs[i], hexs[(i + 1) % 6]
    fa, fc = far[i], far[(i + 1) % 6]
    quads += '    (keepout "" (polygon B.Cu 0  %.1f %.1f  %.1f %.1f  %.1f %.1f  %.1f %.1f))\n' % (a[0], a[1], c[0], c[1], fc[0], fc[1], fa[0], fa[1])
big_ko = "".join('    (keepout "" (polygon %s 0  %s))\n' % (ly, " ".join("%.1f" % v for xy in far for v in xy)) for ly in ("F.Cu", "In1.Cu", "In2.Cu"))
anchor = "    (via "
t = t.replace(anchor, quads + big_ko + anchor, 1)
open(d + "/sub6_board_route.dsn", "w", encoding="utf-8", newline="\n").write(t)
print("outline keepouts: 6 quads, board-wide layer keepouts: 3, apothem(um)", ap, "shrunk by", EDGE_KEEP, "mm", flush=True)
