"""B.Cu PNGs (full board + connector zoom) from the routed KiCad board tracks.  MODEL_ONLY.  python sub6_final_13pad_png.py"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle

R = Path(__file__).resolve().parents[2] / "results" / "sub6_final_board"
D = R / "mesh96_13pad"
tr = json.load(open(D / "kicad" / "tracks.json"))
nets = json.load(open(D / "kicad" / "nets.json"))
rt = json.load(open(R / "mesh96" / "routing.json"))
cells = [c for c in csv.DictReader(open(R / "mesh96" / "fab" / "cell_table.csv")) if c["class"] == "SWITCHED"]


def draw(ax, s):
    ax.add_patch(Rectangle((-5.175, -5.175), 10.35, 10.35, fill=False, ec="m", lw=0.8, ls="--"))   # SMA keep-out
    ax.add_patch(Rectangle((-3.675, -3.675), 7.35, 7.35, fc="#f0c000", ec="none", alpha=0.7))       # SMA ground pad
    ax.add_patch(Circle((0, 0), 2.05, fc="white", ec="k", lw=0.3))                                  # pin opening
    import math as _m
    for k in range(12):
        a = 2 * _m.pi * k / 12 + _m.pi / 12
        ax.add_patch(Circle((2.7 * _m.cos(a), 2.7 * _m.sin(a)), 0.3, fc="#2ca02c", ec="none"))
    ax.add_patch(Polygon(rt["outline"], closed=True, fill=False, ec="k", lw=1))
    for t in tr:
        col = "#d62728" if t["net"] == "RETURN" else ("#1f77b4" if t["net"].startswith("S") else "#555555")
        ax.plot([t["a"][0], t["b"][0]], [t["a"][1], t["b"][1]], color=col, lw=t["w"] * s, solid_capstyle="round")
    for c in cells:
        for k in ("island", "disc"):
            ax.add_patch(Circle((float(c[k + "_via_x"]), float(c[k + "_via_y"])), 0.3, fc="#2ca02c", ec="none"))
    for p in nets["pads"]:
        ax.add_patch(Rectangle((p["x"] - 0.5, p["y"] - 0.5), 1, 1, fc="#ff7f0e", ec="none"))
    ax.set_aspect("equal")


fig, ax = plt.subplots(figsize=(13, 13))
draw(ax, 3.0)
ax.set_title("SUB6 final board B.Cu (13 pads, KiCad + Freerouting): blue S<string> drive 0.25 mm, red RETURN 0.5 mm, grey links 0.25 mm; MODEL_ONLY", fontsize=9)
fig.savefig(D / "fab" / "B_Cu_full.png", dpi=110, bbox_inches="tight")
plt.close(fig)
fig, ax = plt.subplots(figsize=(8, 8))
draw(ax, 25.0)
ax.set_xlim(-9, 9)
ax.set_ylim(-9, 9)
ax.set_title("SMA centre: B.Cu ground pad 7.35 x 7.35 (yellow), pin opening r 2.05, 12 blind GND vias (green), keep-out 10.35 x 10.35 (dashed)", fontsize=8)
fig.savefig(D / "fab" / "B_Cu_SMA_zoom.png", dpi=100, bbox_inches="tight")
plt.close(fig)
fig, ax = plt.subplots(figsize=(11, 11))
draw(ax, 9.0)
ax.set_xlim(-190, -140)
ax.set_ylim(-26, 26)
for p in nets["pads"]:
    ax.text(p["x"] - 1.2, p["y"], p["net"], ha="right", va="center", fontsize=8)
ax.set_title("Connector zoom: 13 pads, S01..S12 drive, RETURN = middle pad (7), pitch 2.54 mm", fontsize=9)
fig.savefig(D / "fab" / "B_Cu_connector_zoom.png", dpi=100, bbox_inches="tight")
print("png written")
