"""Grounded-via BAR90 eigenmode cell. Always launches a new local AEDT session."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

import cell_param_pyaedt as cp
from opt4_geometry import Geometry, GEOMETRIES, check, LINE_W, CHANNEL_W

HERE = Path(__file__).resolve().parent


def configure() -> None:
    cp.set_profile("bar90")
    cp.PROFILE["file"] = "opt1_design_variables.csv"
    cp.DERIVED = ("island_left", "island_cx", "ground_bottom", "air_bottom")
    cp.ROW_TIMEOUT_S = 30 * 60


def values(g: Geometry, state: str = "OFF", direction: str = "across", phase: int = 90,
           lid_height_mm: float = 150) -> dict:
    v = cp.nominal()
    v.update(island_l=g.island_l,
             line_w=LINE_W if g.route != "G-ISL-noline" else 0,
             channel_w=CHANNEL_W if g.route == "G-DISC" else 0,
             air_top=lid_height_mm)
    v["Rs"], v["Xs"] = cp.state_values(state, v)
    v["$ph_x"], v["$ph_y"] = ((phase, 0) if direction == "across" else (0, phase))
    return v


def geometry_named(name: str) -> Geometry:
    for g in GEOMETRIES:
        if g.name == name:
            return g
    raise ValueError(f"Unknown geometry {name!r}")


def build(h, g: Geometry, out: Path, lid_height_mm: float = 150,
          top_boundary: str = "PEC") -> dict:
    dc = check(g)
    if dc["status"] != "PASS":
        raise RuntimeError(f"Geometry conflict: {dc['conflicts']}")
    v = values(g, lid_height_mm=lid_height_mm)
    cp.set_vars(h, v)
    for key, expression in {
        "island_left": "-slot_w/2+diode_gap",
        "island_cx": "island_left+island_w/2",
        "ground_bottom": "-rt_t-cu_t",
        "air_bottom": "ground_bottom",
    }.items():
        h[key] = expression
    h.modeler.model_units = "mm"
    mat = h.materials.add_material("RT5880_opt1")
    mat.permittivity = 2.2
    mat.dielectric_loss_tangent = 0.0009
    m = h.modeler
    box = lambda name, origin, size, material="pec": cp.require(
        m.create_box(origin, size, name=name, material=material), "create " + name)
    disc = cp.require(m.create_cylinder("Z", [0, 0, 0], "disc_r", "cu_t", name="DISC", material="pec"), "disc")
    slot = box("SLOT_CUT", ["-slot_w/2", "-slot_l/2", "-cu_t"],
               ["slot_w", "slot_l", "3*cu_t"], "vacuum")
    cp.require(m.subtract(disc.name, slot.name, keep_originals=False), "slot subtraction")
    if g.route == "G-DISC":
        channel = box("CHANNEL_CUT_TOP", ["island_cx-channel_w/2", "slot_l/2", "-cu_t"],
                      ["channel_w", "period/2-slot_l/2+0.1mm", "3*cu_t"], "vacuum")
        bottom_channel = box("CHANNEL_CUT_BOTTOM",
                             ["island_cx-channel_w/2", "-period/2-0.1mm", "-cu_t"],
                             ["channel_w", "period/2-slot_l/2+0.1mm", "3*cu_t"], "vacuum")
    else:
        channel = None
    if channel:
        cp.require(m.subtract(disc.name, channel.name, keep_originals=False), "channel subtraction")
        if g.route == "G-DISC":
            cp.require(m.subtract(disc.name, bottom_channel.name, keep_originals=False),
                       "bottom channel subtraction")
    island = box("ISLAND", ["island_left", "-island_l/2", 0], ["island_w", "island_l", "cu_t"])
    if g.route == "G-DISC":
        line = box("OUTBOUND_LINE", ["island_cx-line_w/2", "island_l/2", 0],
                   ["line_w", "period/2-island_l/2", "cu_t"])
        receiver = box("RETURN_LINE", ["island_cx-line_w/2", "-period/2", 0],
                       ["line_w", "period/2-island_l/2", "cu_t"])
    elif g.route == "G-ISL-a":
        line = box("OUTBOUND_LINE", ["-line_w/2", "5mm", 0],
                   ["line_w", "period/2-5mm", "cu_t"])
        receiver = box("RECEIVER_TAB", ["-line_w/2", "-period/2", 0],
                       ["line_w", "period/2-5mm", "cu_t"])
    elif g.route == "G-ISL-b":
        line = box("OUTBOUND_LINE", ["5mm", "-line_w/2", 0],
                   ["period/2-5mm", "line_w", "cu_t"])
        receiver = box("RECEIVER_TAB", ["-period/2", "-line_w/2", 0],
                       ["period/2-5mm", "line_w", "cu_t"])
    else:
        line = receiver = None
    if line:
        if g.route == "G-DISC":
            cp.require(m.unite([island.name, line.name, receiver.name]),
                       "continuous island bias bus union")
        else:
            cp.require(m.unite([disc.name, line.name, receiver.name]), "disc periodic line union")
    ground = box("SOLID_GROUND", ["-period/2", "-period/2", "ground_bottom"],
                 ["period", "period", "cu_t"])
    substrate = box("RT5880", ["-period/2", "-period/2", "-rt_t"],
                    ["period", "period", "rt_t"], "RT5880_opt1")
    via_x = "-3mm" if g.route == "G-DISC" else "island_cx"
    via = cp.require(m.create_cylinder("Z", [via_x, "0mm", "-rt_t"],
                                       "0.15mm", "rt_t+cu_t", name="GROUND_VIA",
                                       material="pec"), "ground via")
    pad = cp.require(m.create_cylinder("Z", [via_x, "0mm", 0],
                                       "0.3mm", "cu_t", name="VIA_TOP_PAD",
                                       material="pec"), "via top pad")
    cp.require(m.subtract(substrate.name, via.name, keep_originals=True), "drill RT5880")
    target = disc if g.route == "G-DISC" else island
    cp.require(m.unite([target.name, pad.name, via.name, ground.name]), "via short to solid ground")
    sheet = cp.require(m.create_rectangle("XY", ["-slot_w/2", "-diode_w/2", "cu_t/2"],
                                          ["diode_gap", "diode_w"], name="DIODE_SHEET", material="vacuum"), "diode sheet")
    cp.require(h.assign_impedance_to_sheet(sheet.name, name="BAR90_02LRH",
                                           resistance="Rs", reactance="Xs"), "diode impedance")
    air = box("AIR", ["-period/2", "-period/2", "air_bottom"],
              ["period", "period", "air_top-air_bottom"], "vacuum")
    zmid = (v["air_top"] - v["rt_t"] - v["cu_t"]) / 2
    face = lambda pos: cp.require(m.get_faceid_from_position(pos, assignment="AIR"), "AIR face")
    cp.require(h.assign_lattice_pair([face([-v["period"]/2, 0, zmid]), face([v["period"]/2, 0, zmid])],
                                     phase_delay="InputPhaseDelay", phase_delay_param1="$ph_x", name="Lattice_X"), "x lattice")
    cp.require(h.assign_lattice_pair([face([0, -v["period"]/2, zmid]), face([0, v["period"]/2, zmid])],
                                     phase_delay="InputPhaseDelay", phase_delay_param1="$ph_y", name="Lattice_Y"), "y lattice")
    top_face = face([0, 0, v["air_top"]])
    if top_boundary == "PEC":
        cp.require(h.assign_perfect_e([top_face], name="PEC_LID"), "PEC lid")
    elif top_boundary == "IMPEDANCE_377":
        cp.require(h.assign_impedance_to_sheet([top_face], name="TOP_IMPEDANCE_377",
                                               resistance="377ohm", reactance="0ohm"), "377-ohm top")
    else:
        raise ValueError(f"Unknown top boundary {top_boundary!r}")
    cp.require(h.assign_perfect_e([face([0, 0, -v["rt_t"]-v["cu_t"]])], name="PEC_GROUND_BOTTOM"), "PEC ground bottom")
    setup = cp.require(h.create_setup(cp.SETUP), "eigenmode setup")
    setup.props["MinimumFrequency"] = "4.5GHz"
    setup.props["NumModes"] = 8
    setup.props["MaxDeltaFreq"] = 0.1
    setup.props["MaximumPasses"] = 15
    setup.props["MinimumPasses"] = 3
    setup.props["ConvergeOnRealFreq"] = True
    cp.require(setup.update(), "setup update")
    validation, ok = h.validate_full_design(output_dir=str(out))
    summary = {"status": "BUILT" if ok else "VALIDATION_FAIL", "aedt": h.desktop_class.aedt_version_id,
               "geometry": dc, "objects": m.object_names, "boundaries": [b.name for b in h.boundaries],
               "stack": {"top_cu": True, "rt5880_mm": v["rt_t"], "solid_ground": True,
                         "fr4": False, "vias": True, "via_drill_mm": 0.3, "via_pad_mm": 0.6,
                         "antipads": False, "stubs": False, "air_below_ground": False},
               "floor_boundary": "SOLID_PEC_GROUND", "mesh_operations": [], "validation": validation}
    summary["top_boundary"] = {"type": top_boundary, "height_mm": lid_height_mm}
    if top_boundary == "IMPEDANCE_377":
        boundary = next(b for b in h.boundaries if b.name == "TOP_IMPEDANCE_377")
        summary["top_boundary"]["Resistance"] = boundary.props["Resistance"]
        summary["top_boundary"]["Reactance"] = boundary.props["Reactance"]
        cp.require(boundary.props["Resistance"] == "377ohm" and
                   boundary.props["Reactance"] == "0ohm", "top impedance readback")
    cp.write_json(out / "build_summary.json", summary)
    cp.require(ok, "native design validation")
    cp.require(h.save_project(), "save project")
    return summary


def solve_csv(h, g: Geometry, csv_path: Path, out: Path, cores: int, project: Path,
              lid_height_mm: float = 150, project_variant: str | None = None) -> None:
    cp.PROJECT_VARIANT = project_variant or g.name
    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    failed = []
    try:
        for row in rows:
            run_id = row["*"]
            v = values(g, row["state"], row["direction"], int(row["phase_deg"]), lid_height_mm)
            if (out / f"{run_id}.json").exists():
                continue
            errors = []
            for attempt in (1, 2):
                try:
                    cp.solve_row(h, out, run_id, v, cores)
                    break
                except Exception as exc:
                    errors.append(f"attempt {attempt}: {type(exc).__name__}: {str(exc)[:300]}")
                    print(f"ROW_FAILED {run_id}: {errors[-1]}", flush=True)
                    h.release_desktop()
                    h = cp.connect(project)
            else:
                failed.append({"run_id": run_id, "attempts": 2, "errors": errors,
                               "state": row["state"], "direction": row["direction"],
                               "phase_deg": int(row["phase_deg"]), "geometry": g.name})
                cp.write_json(out / "FAILED_ROWS.json", failed)
                print(f"ROW_SKIPPED_AS_FAILED {run_id}", flush=True)
    finally:
        h.release_desktop()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", required=True)
    ap.add_argument("--project", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--lid-height-mm", type=float, choices=(150, 300), default=150)
    ap.add_argument("--top-boundary", choices=("PEC", "IMPEDANCE_377"), default="PEC")
    ap.add_argument("--project-variant", default=None)
    args = ap.parse_args()
    configure()
    g = geometry_named(args.geometry)
    if check(g)["status"] != "PASS":
        raise RuntimeError(f"Geometry conflict: {check(g)['conflicts']}")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.build and args.project.exists():
        raise RuntimeError("Refusing to overwrite project")
    if not args.build and not args.project.exists():
        raise RuntimeError("Built project missing")
    args.project.parent.mkdir(parents=True, exist_ok=True)
    h = cp.connect(args.project)
    if args.build:
        try:
            print(json.dumps(build(h, g, args.output, args.lid_height_mm, args.top_boundary), indent=2), flush=True)
        finally:
            h.release_desktop()
    elif args.csv:
        solve_csv(h, g, args.csv, args.output, args.cores, args.project,
                  args.lid_height_mm, args.project_variant)
    else:
        h.release_desktop()
        raise RuntimeError("Choose --build or --csv")


if __name__ == "__main__":
    main()
