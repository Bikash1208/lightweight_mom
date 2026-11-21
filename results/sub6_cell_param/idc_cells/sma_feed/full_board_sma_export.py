"""Export the user-supplied v3.1 R150 anchor/gate CSV as a full RF STEP.

The geometry and feed follow export_board_v2_br.py from the R100 crop package.
No HFSS session is opened here. Each state gets an independent STEP.
"""
import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import cadquery as cq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "research"))
from sub6_step_exact import audit_step_text, exactify_step_text

T = 0.035
H = 3.175
T_BIAS = 0.254
R_APERTURE = 150.0
R_ABS_IN = 156.0  # Clears the outermost disc (centre <=150, radius 5.45).
R_BOARD = 186.0  # 30 mm edge absorber follows the R100 crop board package.
R_DISC = 5.45
SLOT_L = 8.72
SLOT_W = 2.91
PITCH = 11.6
PIN_R = 0.635
PORT_R = 2.05
PIN_TOP = 8.925
SMA_PIN_R = 0.3048
SMA_ANTIPAD_R = 0.9017   # matches the r 0.9017 bore of the SMA flange
AIR_BELOW = 25.0
Z_GROUND_BOTTOM = -H - T
Z_BIAS_BOTTOM = Z_GROUND_BOTTOM - T_BIAS


def r3(value):
    return round(float(value), 3)


def slot_points(x, y, rotation_degrees):
    angle = math.radians(rotation_degrees)
    ca, sa = math.cos(angle), math.sin(angle)
    return [
        (r3(x + ca * u - sa * v), r3(y + sa * u + ca * v))
        for u, v in ((-SLOT_L / 2, -SLOT_W / 2), (SLOT_L / 2, -SLOT_W / 2),
                     (SLOT_L / 2, SLOT_W / 2), (-SLOT_L / 2, SLOT_W / 2))
    ]


def short_slot_points(x, y, rotation_degrees, u_min, u_max):
    """One closed slot beside a 0.5 mm centre copper strip, in the cell frame."""
    angle = math.radians(rotation_degrees)
    ca, sa = math.cos(angle), math.sin(angle)
    return [
        (r3(x + ca * u - sa * v), r3(y + sa * u + ca * v))
        for u, v in ((u_min, -SLOT_W / 2), (u_max, -SLOT_W / 2),
                     (u_max, SLOT_W / 2), (u_min, SLOT_W / 2))
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--state", choices=("BOTH", "A_ONLY", "B_ONLY", "A_ONLY_IDEAL",
                                             "B_ALL_OFF_IDEAL", "ALL_PLAIN", "B_ALL_SHORT"), default="BOTH")
    parser.add_argument("--feed", choices=("probe", "sma"), default="probe")
    parser.add_argument("--rotation-mode", choices=("input", "clean_swap", "scrambled"), default="input")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.csv_path.open(newline="")))
    assert len(rows) == 516, len(rows)
    assert sum(int(row["gate"]) for row in rows) == 16
    assert sum(row["beam"] == "A" for row in rows) == 258
    assert sum(row["beam"] == "B" for row in rows) == 258
    assert sum(row["gate_for"] == "A" for row in rows) == 8
    assert sum(row["gate_for"] == "B" for row in rows) == 8
    assert all(r3(float(row["x_mm"])) == r3(int(row["i"]) * PITCH) for row in rows)
    assert all(r3(float(row["y_mm"])) == r3(int(row["j"]) * PITCH) for row in rows)
    t0 = time.time()
    metal_parts = []
    via_xy = []
    layout_readback = []
    ideal_changed = 0
    shorted_volumes = []
    for index, row in enumerate(rows):
        x, y = float(row["x_mm"]), float(row["y_mm"])
        assert math.hypot(x, y) <= R_APERTURE + 1e-8
        assert math.hypot(x, y) + R_DISC < R_ABS_IN
        gate = bool(int(row["gate"]))
        shorted = args.state == "B_ALL_SHORT" and row["beam"] == "B" and not gate
        ideal_plain = ((args.state == "A_ONLY_IDEAL" and row["beam"] == "B"
                        and x <= 0 and float(row["r_mm"]) >= 29)
                       or (args.state == "B_ALL_OFF_IDEAL" and row["beam"] == "B")
                       or args.state == "ALL_PLAIN")
        if ideal_plain and not gate:
            ideal_changed += 1
        source_angle = float(row["rot_deg"])
        if args.rotation_mode == "clean_swap":
            angle = (math.degrees(math.atan2(y, x)) - source_angle) % 180
        elif args.rotation_mode == "scrambled":
            angle = (-source_angle) % 180
        else:
            angle = source_angle
        work = cq.Workplane("XY").center(x, y).circle(R_DISC)
        if shorted:
            for u_min, u_max in ((-SLOT_L / 2, -0.25), (0.25, SLOT_L / 2)):
                points = short_slot_points(x, y, angle, u_min, u_max)
                work = work.polyline([(px - x, py - y) for px, py in points]).close()
        elif not gate and not ideal_plain:
            points = slot_points(x, y, angle)
            work = work.polyline([(px - x, py - y) for px, py in points]).close()
        disc = work.extrude(T).val()
        if shorted:
            expected = (math.pi * R_DISC**2 - 2 * 4.11 * SLOT_W) * T
            assert abs(disc.Volume() - expected) < 0.002, (index, disc.Volume(), expected)
            shorted_volumes.append(disc.Volume())
        block = gate and ((args.state == "A_ONLY" and row["gate_for"] == "B") or
                          (args.state == "B_ONLY" and row["gate_for"] == "A"))
        if block:
            via = cq.Workplane("XY").workplane(offset=-H).center(x, y).rect(0.3, 0.3).extrude(H).val()
            disc = disc.fuse(via).clean()
            via_xy.append((x, y))
        metal_parts.append(disc)
        layout_readback.append({"i": int(row["i"]), "j": int(row["j"]), "xy_mm": [x, y],
                                "gate": gate, "beam": row["beam"], "gate_for": row["gate_for"], "block": block,
                                "ideal_plain_replacement": ideal_plain and not gate,
                                "plain_target": ideal_plain,
                                "centre_short_strip": shorted,
                                "rotation_deg": angle if not gate else None})
        if index % 100 == 99:
            print("DISCS", index + 1, "seconds", round(time.time() - t0, 1), flush=True)

    ground_work = cq.Workplane("XY").workplane(offset=Z_GROUND_BOTTOM).circle(R_BOARD)
    sma = args.feed == "sma"
    ground_work = ground_work.moveTo(0, 0).circle(SMA_ANTIPAD_R if sma else PORT_R)
    ground = ground_work.extrude(T).val()
    if sma:
        pin = cq.Workplane("XY").workplane(offset=-H).circle(SMA_PIN_R).extrude(PIN_TOP + H).val()   # only used to cut the substrate; the SMA STEP supplies the real pin
        metal_parts += [ground]
    else:
        pin = cq.Workplane("XY").workplane(offset=-H).circle(PIN_R).extrude(PIN_TOP + H).val()
        metal_parts += [ground, pin]
    metal = cq.Compound.makeCompound(metal_parts)
    sub = cq.Workplane("XY").workplane(offset=-H).circle(R_BOARD).extrude(H).val()
    if via_xy:
        cuts = [cq.Workplane("XY").workplane(offset=-H).center(x, y).rect(0.3, 0.3).extrude(H).val() for x, y in via_xy]
        sub = sub.cut(*cuts).clean()
    sub = sub.cut(pin).clean()
    bias = cq.Workplane("XY").workplane(offset=Z_BIAS_BOTTOM).circle(R_BOARD).extrude(T_BIAS).val()
    window = cq.Workplane("XY").workplane(offset=Z_BIAS_BOTTOM - 1).rect(12, 12).extrude(T_BIAS + 2).val()
    bias = bias.cut(window).clean()
    air_bot = (-AIR_BELOW) if args.feed == "sma" else (Z_BIAS_BOTTOM - 15)
    air = cq.Workplane("XY").workplane(offset=air_bot).rect(2 * (R_BOARD + 15), 2 * (R_BOARD + 15)).extrude(T + 40 - air_bot).val()
    asm = cq.Assembly(None, name="BOARD_R150_GATE")
    items = [("METAL_PEC", metal), ("SUBSTRATE", sub), ("BIAS_RT5880", bias), ("AIR", air)]
    if not sma:
        port = cq.Face.makeFromWires(cq.Wire.makeCircle(PORT_R, cq.Vector(0, 0, -H), cq.Vector(0, 0, 1)),
                                     [cq.Wire.makeCircle(PIN_R, cq.Vector(0, 0, -H), cq.Vector(0, 0, 1))])
        items.append(("PORT", port))
    for name, shape in items:
        asm.add(shape, name=name)
    tag = args.state + ("_SMA" if sma else "") + ("_" + args.rotation_mode.upper() if args.rotation_mode != "input" else "")
    step = args.output_dir / f"BOARD_R150_GATE_V31_{tag}.step"
    asm.export(str(step), exportType="STEP")
    snapped, snap_stats = exactify_step_text(step.read_text(encoding="utf-8"))
    step.write_text(snapped, encoding="utf-8", newline="\n")
    bad = audit_step_text(step.read_text(encoding="utf-8"))
    if args.state == "A_ONLY_IDEAL":
        assert ideal_changed == 234 and not via_xy
    elif args.state == "B_ALL_OFF_IDEAL":
        assert ideal_changed == 250 and not via_xy
    elif args.state == "ALL_PLAIN":
        assert ideal_changed == 500 and not via_xy
    elif args.state == "B_ALL_SHORT":
        assert len(shorted_volumes) == 250 and ideal_changed == 0 and not via_xy
    manifest = {"state": args.state, "rotation_mode": args.rotation_mode, "sites": len(rows),
                "gate_sites": 16, "blocked_vias": len(via_xy), "board_radius_mm": R_BOARD,
                "ideal_plain_replacements": ideal_changed,
                "plain_target_sites": sum(row["plain_target"] for row in layout_readback),
                "centre_short_strip_cells": len(shorted_volumes),
                "centre_short_strip_width_mm": 0.5 if shorted_volumes else None,
                "resulting_slot_lengths_mm": [4.11, 4.11] if shorted_volumes else None,
                "shorted_top_copper_volume_mm3_minmax": [min(shorted_volumes), max(shorted_volumes)] if shorted_volumes else None,
                "aperture_radius_mm": R_APERTURE, "absorber_ring_mm": [R_ABS_IN, R_BOARD],
                "substrate": {"material": "RT5880", "thickness_mm": H},
                "bias": {"material": "RT5880", "thickness_mm": T_BIAS},
                "ground": {"material": "PEC", "thickness_mm": T, "feed_opening_radius_mm": PORT_R},
                "feed": {"pin_radius_mm": PIN_R, "pin_z_mm": [-H, PIN_TOP], "port_z_mm": -H},
                "slot": {"length_mm": SLOT_L, "width_mm": SLOT_W, "disc_radius_mm": R_DISC},
                "step": {"path": str(step), "bytes": step.stat().st_size, "inexact_literals": len(bad)},
                "snap_stats": snap_stats, "rows": layout_readback,
                "dither_applied": False, "alpha_profile_applied": False,
                "elapsed_s": time.time() - t0}
    (args.output_dir / f"manifest_{tag}.json").write_text(json.dumps(manifest, indent=1))
    print("EXPORTED", step, step.stat().st_size, "inexact", len(bad), "seconds", round(time.time()-t0, 1), flush=True)
    if bad:
        raise RuntimeError(f"STEP audit found {len(bad)} inexact literals")


if __name__ == "__main__":
    main()
