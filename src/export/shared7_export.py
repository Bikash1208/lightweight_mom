"""Phase 29 design point B3: AEDT scripts for the 7-port shared 150 mm hexagon (one hologram, 7 monopoles near the centre). Reuses the verified-style templates of
`compact7_export` (same calls as the Phase 26 board that ran in AEDT; wrapped, print-only fallbacks for everything unverified), with the board of the single
150 mm hexagon (270 x 310 mm) and the feeds of the GA result. MODEL_ONLY."""

from .compact7_export import _BUILD, D_MM, H_MM, MONO_H, R_OUT, R_PIN, _fill, aedt_weights_script  # noqa: F401  (re-exported: the weights script is identical)

BOARD_X, BOARD_Y = (
    270.0,
    310.0,
)  # hexagon (circumradius 150 mm, vertex at +/-y) bounding box 259.81 x 300 mm + 5 mm margin per side


def aedt_build_script(cells, feeds_mm, name="phase29_B3"):
    rows = ",\n    ".join(
        "({:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.4f})".format(
            c["x_mm"], c["y_mm"], c["slot_L_mm"], c["slot_W_mm"], c["slot_angle_deg"]
        )
        for c in cells
    )
    feeds = "[" + ", ".join("({:.3f}, {:.3f})".format(*f) for f in feeds_mm) + "]"
    s = _fill(
        _BUILD,
        title=f"shared-aperture 7-port tensor-CP board (150 mm hexagon, {len(cells)} cells)",
        NAME=name,
        H=H_MM,
        WX=f"{BOARD_X:g}",
        WY=f"{BOARD_Y:g}",
        FEEDS=feeds,
        D=D_MM,
        ROWS=rows,
        MONO_H=MONO_H,
        R_PIN=R_PIN,
        R_OUT=R_OUT,
    )
    s = s.replace(
        "the seven hexagons + 5 mm margin, centred on the cluster centre",
        "the 150 mm hexagon (vertex at +/-y) + 5 mm margin, centred on the hexagon centre",
    )
    s = s.replace(
        "ports 1..7: 1 = centre tile (broadside beam), 2..7 = outer tiles, beam azimuth 0, 60, ..., 300 deg at theta 30",
        "ports 1..7: 1 = centre feed (broadside beam), 2..7 = ring feeds, beam azimuth 0, 60, ..., 300 deg (tilt per README)",
    )
    return s
