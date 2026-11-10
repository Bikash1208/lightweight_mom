"""Native HFSS periodic eigenmode study; no MoM ON-state calculations.

Documented API references and limitations: API_REFERENCES.md and RUN.md.
All results remain NOT_HFSS_VALIDATED until nominal VERIFY passes.
"""
from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
SETUP = "EigenCell"
AEDT_VERSION = os.environ.get("SUB6_AEDT_VERSION", "2025.1")   # set on a machine with another install (and report it); default = reference version
TARGETS = {("OFF", "along"): (5.771, 147), ("OFF", "across"): (6.202, 180),
           ("ON", "along"): (4.835, 115), ("ON", "across"): (4.905, 79)}

# Diode profiles. bar50 = user's calibration cell (VERIFY/artefact model check only);
# bar90 = Infineon BAR90-02LRH (TSLP-2-7) for SCREEN. Sheet values are at the profile's own W/L.
PROFILES = {
    "bar50": {"file": "design_variables_bar50.csv", "boundary": "BAR50_02V", "gap": 1.00, "w": 0.35,
              "states": {"OFF": (5.424, -79.320), "ON": (1.05, 7.65)}},
    "bar90": {"file": "design_variables_bar90.csv", "boundary": "BAR90_02LRH", "gap": 0.30, "w": 0.60,
              "states": {"OFF": (22.23, -335.4), "ON": (1.60, 29.2)}},
}
PROFILE = None
PROJECT_VARIANT = "full"   # set from the project build_summary.json in main()
ROW_TIMEOUT_S = 0   # set from --row-timeout-min


def bar90_sheet(c_pf=0.15, rp=3000.0, r_on=0.8, ls_nh=0.4, wl=2.0, f_hz=5.8e9, on=False):
    """BAR90-02LRH 5.8 GHz sheet impedance: (C || Rp) + Ls, or r_on + Ls, scaled by W/L."""
    w = 2 * math.pi * f_hz
    xl = w * ls_nh * 1e-9
    if on:
        return r_on * wl, xl * wl
    y = complex(1 / rp, w * c_pf * 1e-12)
    z = 1 / y + 1j * xl
    return z.real * wl, z.imag * wl


# Sensitivity states (nominal geometry only); values derived from the device model.
SENS_STATES = {"OFF_C013": bar90_sheet(c_pf=0.13), "OFF_C018": bar90_sheet(c_pf=0.18),
               "OFF_RP2K": bar90_sheet(rp=2000.0), "OFF_RP5K": bar90_sheet(rp=5000.0)}


def set_profile(name):
    global PROFILE
    PROFILE = dict(PROFILES[name], name=name)


def states():
    return PROFILE["states"]


def sheet_scale(v):
    """Rescale sheet impedance if W/L differs from the profile's nominal W/L."""
    return (v["diode_w"] / v["diode_gap"]) / (PROFILE["w"] / PROFILE["gap"])


def state_values(state, v=None):
    rs, xs = {**states(), **(SENS_STATES if PROFILE["name"] == "bar90" else {})}[state]
    k = sheet_scale(v) if v else 1.0
    return rs * k, xs * k


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def variables():
    with (HERE / PROFILE["file"]).open(newline="", encoding="utf-8-sig") as f:
        return {r["name"]: r for r in csv.DictReader(f)}


def nominal():
    return {k: float(v["nominal"]) for k, v in variables().items()}


def fabrication(v):
    """Analytic 2D clearances, including circle/rectangle/slot and bottom traces.

    Via is solid PEC in the authoritative calibration, not a plated hollow tube.
    Air is excluded from solid-overlap checks by design.
    """
    p, r = v["period"], v["disc_r"]
    sl, sw = v["slot_l"], v["slot_w"]
    il, iw, off = v["island_l"], v["island_w"], v["island_offset"]
    left = -sw / 2 + v["diode_gap"]
    vr, pr, ar = v["via_drill"] / 2, v["pad_dia"] / 2, v["antipad_dia"] / 2
    ix, iy = v["island_via_x"], v["island_via_y"]
    dx, dy = v["disc_via_x"], v["disc_via_y"]
    # Distance of outside point to the slot rectangle (disc via must stay outside).
    dist_slot = math.hypot(max(abs(dx) - sw / 2, 0), max(abs(dy) - sl / 2, 0))
    # Longest diode feature along y: the sheet, or the copper pad when the profile defines one.
    ext = max(v["diode_w"], v.get("pad_l", 0.0))
    # Slot corners must remain inside the disc to preserve the closed slot topology.
    clear = {
        "annular_ring": pr - vr,
        "ground_to_barrel": ar - vr,
        "disc_via_to_disc_rim": r - math.hypot(dx, dy) - vr,
        "disc_via_to_slot": dist_slot - vr,
        "island_via_to_island": min(ix-left, left+iw-ix, iy-(off-il/2), off+il/2-iy)-vr,
        "island_to_slot_right": sw / 2 - (left + iw),
        "island_to_slot_y": (sl-il)/2-abs(off),
        "diode_end_y": il/2-abs(v["diode_pos"])-ext/2,
        "diode_to_slot_end": sl/2-abs(off+v["diode_pos"])-ext/2,
        "slot_corner_to_disc_rim": r-math.hypot(sl/2, sw/2),
        "disc_to_periodic_disc": p-2*r,
        "bottom_pad_separation": math.hypot(ix-dx, iy-dy)-2*pr,
        "antipad_separation": math.hypot(ix-dx, iy-dy)-2*ar,
    }
    # Each finite rectangular stub plus circular pad must fit and stay apart.
    for name, x, y in [("island", ix, iy), ("disc", dx, dy)]:
        clear[name+"_bottom_to_cell_edge"] = min(p/2-abs(x)-max(pr,v["stub_w"]/2),
              p/2-(y+v["stub_l"]), y-pr+p/2)
    clear["stub_to_stub_x"] = abs(dx-ix)-v["stub_w"]
    # Diode sheet must not intersect the island via at its island-side edge.
    sx0, sx1 = -sw/2, left
    sy = off + v["diode_pos"]
    clear["sheet_to_island_via"] = math.hypot(max(sx0-ix, 0, ix-sx1),
        max(abs(iy-sy)-v["diode_w"]/2, 0))-vr
    if "pad_w" in v:
        # Pads sit on the island and the patch edge; both must be wholly inside their metal.
        clear["pad_inside_island_x"] = iw - v["pad_w"]
        clear["patch_pad_to_disc_rim"] = r - math.hypot(sw/2 + v["pad_w"], abs(off+v["diode_pos"]) + v["pad_l"]/2)
    failures = {k: a for k, a in clear.items() if a < .15-1e-9}
    # Sheet must lie within the pad footprint: gap equals pad gap, width no longer than the pad.
    if "pad_l" in v and v["diode_w"] > v["pad_l"] + 1e-9:
        failures["sheet_wider_than_pad"] = v["diode_w"]
    return {"status": "PASS" if not failures else "FAB_FAIL", "clearances_mm": clear,
            "failures": failures}


FAB_WAIVERS = {}   # DISABLED 2026-11-10: the user dropped the only waived case (anti-pad 0.4 mm); nothing is waivable now


def fab_gate(v, run_id=""):
    g = fabrication(v)
    fails = dict(g["failures"]); waived = {}
    for tag, keys in FAB_WAIVERS.items():
        if tag in run_id:
            for k in keys:
                if k in fails:
                    waived[k] = fails.pop(k)
    return {"status": "PASS" if not fails else "FAB_FAIL", "clearances_mm": g["clearances_mm"], "failures": fails, "waived": waived}


def require(value, context):
    if value is False or value is None:
        raise RuntimeError(context)
    return value


def set_vars(h, v):
    spec = variables()
    for k, val in v.items():
        h[k] = f"{val:.12g}{spec[k]['unit']}"


def connect(project, new=False):
    from ansys.aedt.core import Hfss
    return Hfss(project=str(Path(project).resolve()), design="CellParam", solution_type="Eigenmode",
                version=AEDT_VERSION, non_graphical=True, new_desktop=True, close_on_exit=False)


def build(h, out, floor="pec", air_below_mm=30.0, floor_z_mm=None, mesh_ops=True, variant="full"):
    v = nominal()
    v["air_below"] = float(air_below_mm)
    gate = fabrication(v)
    require(gate["status"] == "PASS", "Nominal fabrication conflict")
    if gate["status"] != "PASS":
        raise RuntimeError(str(gate))
    set_vars(h, v)
    derived = {"island_left": "-slot_w/2+diode_gap", "island_cx": "island_left+island_w/2",
               "ground_bottom": "-rt_t-cu_t", "fr4_bottom": "ground_bottom-fr4_t",
               "metal_bottom": "fr4_bottom-cu_t",
               "air_bottom": (f"{float(floor_z_mm):.12g}mm" if floor_z_mm is not None else "fr4_bottom-air_below")}
    for k, expression in derived.items():
        h[k] = expression
    h.modeler.model_units = "mm"
    for name, eps, loss in [("RT5880_cell", 2.2, .0009), ("FR4_epoxy_cell", 4.4, .02)]:
        mat = h.materials.add_material(name)
        mat.permittivity = eps
        mat.dielectric_loss_tangent = loss
    m = h.modeler
    def box(name, origin, sizes, material="pec"):
        return require(m.create_box(origin, sizes, name=name, material=material), "create "+name)
    def cyl(name, x, y, z, radius, height, material="pec"):
        return require(m.create_cylinder("Z", [x,y,z], radius, height, name=name, material=material), "create "+name)
    patch = cyl("PATCH", 0, 0, 0, "disc_r", "cu_t")
    cut = box("SlotCut", ["-slot_w/2","-slot_l/2","-cu_t"], ["slot_w","slot_l","3*cu_t"], "vacuum")
    require(m.subtract(patch.name, cut.name, keep_originals=False), "slot subtraction")
    island = box("ISLAND", ["island_left","island_offset-island_l/2",0], ["island_w","island_l","cu_t"])
    ground = box("GROUND", ["-period/2","-period/2","ground_bottom"], ["period","period","cu_t"])
    rt = box("SUBSTRATE", ["-period/2","-period/2","-rt_t"], ["period","period","rt_t"], "RT5880_cell")
    fr = box("SUBSTRATE_FR4", ["-period/2","-period/2","fr4_bottom"], ["period","period","fr4_t"], "FR4_epoxy_cell")
    metal = [patch.name, island.name, ground.name]
    fine = []
    # Diagnostic variants (mode identification only): "novia" = no barrels/pads/stubs/anti-pads/holes,
    # "nostub" = vias, pads, anti-pads kept but both bottom stubs removed.
    for tag in ([] if variant == "novia" else ["island", "disc"]):
        x, y = tag+"_via_x", tag+"_via_y"
        barrel = cyl(tag+"_BARREL", x,y,"metal_bottom","via_drill/2","-metal_bottom")
        pad = cyl(tag+"_PAD", x,y,"metal_bottom","pad_dia/2","cu_t")
        metal += [barrel.name, pad.name]
        if variant == "full":
            stub = box(tag+"_STUB", [x+"-stub_w/2",y,"metal_bottom"], ["stub_w","stub_l","cu_t"])
            metal.append(stub.name)
        anti = cyl(tag+"_ANTIPAD",x,y,"ground_bottom-cu_t","antipad_dia/2","3*cu_t","vacuum")
        require(m.subtract(ground.name,anti.name,keep_originals=False), "ground anti-pad")
        for body, z, height in [(rt, "-rt_t-cu_t", "rt_t+2*cu_t"),
                                (fr, "fr4_bottom-cu_t", "fr4_t+2*cu_t")]:
            hole = cyl(tag+"_HOLE_"+body.name,x,y,z,"via_drill/2",height,"vacuum")
            require(m.subtract(body.name,hole.name,keep_originals=False), "substrate hole")
    if variant == "trace":
        # Realistic periodic bias trace (user request): replaces the stubs. Piece A leaves the island via straight to y=disc_via_y, then
        # along -x to the -x lattice face; piece B comes in from the +x lattice face along y=disc_via_y to this cell's disc via.
        # Lattice periodicity makes A of one cell and B of its -x neighbour one continuous series-string trace. Direction/length are
        # an ASSUMPTION (routed board data not available) - flagged to the user.
        for tname, org, size in [
                ("TRACE_A_V", ["island_via_x-trace_w/2", "disc_via_y-trace_w/2", "metal_bottom"],
                              ["trace_w", "island_via_y-disc_via_y+trace_w/2", "cu_t"]),
                ("TRACE_A_H", ["-period/2", "disc_via_y-trace_w/2", "metal_bottom"],
                              ["island_via_x+trace_w/2+period/2", "trace_w", "cu_t"]),
                ("TRACE_B", ["disc_via_x", "disc_via_y-trace_w/2", "metal_bottom"],
                            ["period/2-disc_via_x", "trace_w", "cu_t"])]:
            metal.append(box(tname, org, size).name)
    # One named metal object with three disconnected conducting components, as STEP standard.
    require(m.unite(metal), "METAL_PEC unite")
    m[metal[0]].name = "METAL_PEC"
    m["METAL_PEC"].solve_inside = False
    sheet = require(m.create_rectangle("XY", ["-slot_w/2","island_offset+diode_pos-diode_w/2","cu_t/2"],
                 ["diode_gap","diode_w"], name="DIODE_SHEET", material="vacuum"), "diode sheet")
    require(h.assign_impedance_to_sheet(sheet.name, name=PROFILE["boundary"], resistance="Rs", reactance="Xs"), "impedance boundary")
    air = box("AIR", ["-period/2","-period/2","air_bottom"],
              ["period","period","air_top-air_bottom"], "vacuum")
    # Use nominal face centres, never unstable guessed face indices.
    zmid = (v["air_top"]-v["rt_t"]-2*v["cu_t"]-v["fr4_t"]-v["air_below"])/2
    face = lambda pos: require(m.get_faceid_from_position(pos, assignment="AIR"), "AIR face")
    px = [face([-v["period"]/2,0,zmid]), face([v["period"]/2,0,zmid])]
    py = [face([0,-v["period"]/2,zmid]), face([0,v["period"]/2,zmid])]
    require(h.assign_lattice_pair(px, phase_delay="InputPhaseDelay", phase_delay_param1="$ph_x",name="Lattice_X"), "x lattice")
    require(h.assign_lattice_pair(py, phase_delay="InputPhaseDelay", phase_delay_param1="$ph_y",name="Lattice_Y"), "y lattice")
    require(h.assign_perfect_e([face([0,0,v["air_top"]])],name="PEC_LID"), "PEC top")
    bottom = float(floor_z_mm) if floor_z_mm is not None else -v["rt_t"]-v["cu_t"]-v["fr4_t"]-v["air_below"]
    require(h.save_project(), "save partial geometry before bottom boundary")
    bottom_face = face([0,0,bottom])
    if floor == "pec":
        # Calibration setup: deliberately unassigned. HFSS's outer default is PEC.
        floor_boundary = "UNASSIGNED_DEFAULT_PEC"
    elif floor == "impedance377":
        require(h.assign_impedance_to_sheet([bottom_face], name="ABSORBING_FLOOR_377",
                    resistance="377ohm", reactance="0ohm"), "377-ohm floor")
        floor_boundary = "IMPEDANCE_377_OHM_PER_SQUARE"
    else:
        raise ValueError(f"Unknown floor {floor!r}")
    # No implicit element-count limit; dimensions remain AEDT variable expressions.
    applied_mesh_ops = []
    for names, size, inside, label in ([] if not mesh_ops else [([rt.name],"mesh_rt",True,"RT_mesh"),
                ([fr.name],"mesh_fr4",True,"FR4_mesh"),
                # The united metal contains the top pattern, via barrels, pads and stubs.
                # Apply the calibration's 0.1-mm fine limit to the whole disconnected compound;
                # this is conservative for top metal and avoids unreliable curved-face probing.
                (["METAL_PEC"],"mesh_fine",False,"Metal_mesh"),
                ([sheet.name],"mesh_fine",False,"Diode_mesh")]):
        require(h.mesh.assign_length_mesh(names,inside_selection=inside,maximum_length=size,
                       maximum_elements=None,name=label), "mesh "+label)
        applied_mesh_ops.append(label)
    setup = require(h.create_setup(SETUP), "eigenmode setup")
    setup.props["MinimumFrequency"] = "4.5GHz"
    setup.props["NumModes"] = 8
    setup.props["MaxDeltaFreq"] = .1
    setup.props["MaximumPasses"] = 15
    setup.props["MinimumPasses"] = 3
    setup.props["ConvergeOnRealFreq"] = True
    require(setup.update(), "setup update")
    info, ok = h.validate_full_design(output_dir=str(out))
    summary = {"status": "BUILT" if ok else "VALIDATION_FAIL", "aedt": h.desktop_class.aedt_version_id,
        "objects": m.object_names, "boundaries": [b.name for b in h.boundaries],
        "air_faces": len(air.faces), "air_bottom_mm": bottom, "floor_boundary":floor_boundary,
        "frame": "slot y; across x; along y",
        "mesh_operations": applied_mesh_ops, "variant": variant, "diode_profile": PROFILE["name"], "sheet_mm": [v["diode_gap"], v["diode_w"]],
        "sheet_nominal_rs_xs": [v["Rs"], v["Xs"]], "fabrication": gate, "validation": info}
    write_json(out/"build_summary.json",summary)
    print(json.dumps(summary, indent=2),flush=True)
    require(ok, "Native geometry validation failed")
    require(h.save_project(), "save project")
    return summary


def read_rows(path):
    with path.open(newline="",encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    spec = variables()
    if not rows or "*" not in rows[0]:
        raise ValueError("Expected * run ID and AEDT variable columns")
    seen = set()
    for row in rows:
        if row["*"] in seen:
            raise ValueError("Duplicate run ID")
        seen.add(row["*"])
        unknown = set(row)-set(spec)-{"*"}
        if unknown:
            raise ValueError(f"Unknown parameters: {unknown}")
        v = nominal()
        for k, val in row.items():
            if k != "*" and val != "":
                v[k] = float(val)
                if not float(spec[k]["min"])-1e-9 <= v[k] <= float(spec[k]["max"])+1e-9:
                    raise ValueError(f"{k} outside declared range in {row['*']}")
        g = fab_gate(v, row["*"])
        if g["status"] != "PASS":
            raise ValueError(f"FAB_FAIL {row['*']}: {g['failures']}")
        yield row["*"], v


def case_metadata(v):
    allowed = list(states()) + (list(SENS_STATES) if PROFILE["name"] == "bar90" else [])
    state = next((s for s in allowed if all(abs(a-b)<1e-6 for a,b in zip((v["Rs"],v["Xs"]), state_values(s,v)))), None)
    if state is None:
        raise ValueError("Unknown diode state")
    x,y = v["$ph_x"],v["$ph_y"]
    if (x == 0) == (y == 0):
        raise ValueError("Exactly one nonzero phase required")
    return state, "across" if x else "along", max(x,y)


UNIT_SI = {"mm": 1e-3, "deg": math.pi/180, "ohm": 1.0, "": 1.0}
DERIVED = ("island_left", "island_cx", "ground_bottom", "fr4_bottom", "metal_bottom", "air_bottom")


def values_match(want, got, unit):
    """AEDT returns SI (mm->m, deg->rad); also accept the raw number in the variable's own unit."""
    if not isinstance(got, (int, float)) or isinstance(got, bool) or not math.isfinite(got):
        return False
    return any(math.isclose(got, c, rel_tol=1e-9, abs_tol=1e-12) for c in (want*UNIT_SI.get(unit, 1.0), want))


def read_back(h, v):
    """Read every intended variable back from AEDT. Returns (readback_si, mismatches)."""
    spec = variables(); rb = {}; bad = {}
    for k, want in v.items():
        got = h.get_evaluated_value(k)
        rb[k] = got if isinstance(got, (int, float)) and not isinstance(got, bool) else None
        if not values_match(want, got, spec[k]["unit"]):
            bad[k] = {"intended": want, "read_back": got}
    for k in DERIVED:                       # derived geometry expressions, kept for the record
        got = h.get_evaluated_value(k)
        rb[k] = got if isinstance(got, (int, float)) and not isinstance(got, bool) else None
    # Floor depth consistency when the AIR bottom is parametric (artefact projects): fr4_bottom - air_below.
    try:
        param = "air_below" in str(h["air_bottom"])
    except Exception:
        param = False
    if param:
        want_mm = -v["rt_t"] - v["cu_t"] - v["fr4_t"] - v["air_below"]
        if not values_match(want_mm, rb["air_bottom"], "mm"):
            bad["air_bottom"] = {"intended": want_mm, "read_back": rb["air_bottom"]}
    return rb, bad


def readback_ok(result):
    """Re-check a stored result row. Returns (ok, reason). Rows without read-back are NOT accepted."""
    rb = result.get("readback_si")
    if not rb:
        return False, "no read-back stored (solved before the read-back rule)"
    spec = variables()
    for k, want in result["variables"].items():
        if not values_match(want, rb.get(k), spec[k]["unit"]):
            return False, f"{k}: intended {want}, read back {rb.get(k)}"
    return True, ""


def solve_row(h, out, run_id, v, cores):
    state,direction,phase = case_metadata(v)
    set_vars(h,v)
    rb, bad = read_back(h, v)
    if bad:
        raise RuntimeError(f"READBACK_MISMATCH {run_id}: {json.dumps(bad)}")
    # DeleteFullVariation is invalid before the first solution exists in AEDT 2025 R1.
    if getattr(h, "_sub6_has_solved", False):
        require(h.cleanup_solution(), "Clean Up Solutions")
    print(f"SOLVE {run_id}: {state} {direction} {phase}deg; requested cores={cores}",flush=True)
    require(h.save_project(), "save before solve")
    start=time.time()
    require(h.analyze_setup(SETUP,cores=cores,use_auto_settings=True,blocking=False), "HFSS solve failed to start; see AEDT messages/license diagnostics")
    # Hard per-row limit (covers meshing + adaptive passes). Poll instead of blocking so it can be enforced.
    time.sleep(20)
    while h.desktop_class.are_there_simulations_running:
        if ROW_TIMEOUT_S and time.time()-start>ROW_TIMEOUT_S:
            try:
                h.stop_simulations(clean_stop=False)
            except Exception as exc:
                print("stop_simulations raised:",exc,flush=True)
            raise RuntimeError(f"ROW_TIMEOUT: {run_id} exceeded {ROW_TIMEOUT_S/60:.0f} min after {time.time()-start:.0f} s; solve stopped")
        time.sleep(15)
    h._sub6_has_solved = True
    # Official PyAEDT Eigenmode example: query names instead of guessing EigenMode(i) syntax.
    qs = h.post.available_report_quantities(quantities_category="Eigen Q")
    fs = h.post.available_report_quantities(quantities_category="Eigen Modes")
    if len(qs) != len(fs) or len(fs) < 8:
        raise RuntimeError(f"Incomplete mode export: frequencies={fs}; Q={qs}")
    modes=[]
    for i,(fq,qq) in enumerate(zip(fs,qs),1):
        freq=h.post.get_solution_data(expressions=fq,setup_sweep_name=SETUP+" : LastAdaptive",report_category="Eigenmode")
        q=h.post.get_solution_data(expressions=qq,setup_sweep_name=SETUP+" : LastAdaptive",report_category="Eigenmode")
        fv=float(freq.get_expression_data()[1][0])
        qv=float(q.get_expression_data()[1][0])
        unit=freq.units_data.get(fq, "")
        scale={"Hz":1e-9,"kHz":1e-6,"MHz":1e-3,"GHz":1}.get(unit)
        resolved_by="units_data"
        if scale is None:
            # AEDT 2025 R1 returns an empty unit for Eigen Modes. The setup fixes the minimum
            # eigenfrequency at 4.5 GHz, and SI prefixes are 1000x apart, so at most one prefix can put the
            # raw value in 4-100 GHz. Accept only a unique bracket; record it. Otherwise refuse.
            fits=[sc for sc in (1e-9,1e-6,1e-3,1.0) if 4.0<=fv*sc<=100.0]
            if len(fits)!=1:
                raise RuntimeError(f"Unknown reported frequency units {unit!r} and raw value {fv!r} is not unambiguous")
            scale=fits[0]; resolved_by="setup_min_frequency_bracket(4.5GHz..100GHz)"
        if not math.isfinite(fv) or not math.isfinite(qv) or qv<=0:
            raise RuntimeError("Failed physics check: invalid eigenfrequency or Q")
        modes.append({"mode":i,"f_ghz":fv*scale,"q":qv,"frequency_expression":fq,"q_expression":qq,
                      "raw_frequency":fv,"raw_unit":unit,"unit_resolved_by":resolved_by})
    h.export_profile(SETUP,output_file=str(out/(run_id+".prof")))
    result={"run_id":run_id,"state":state,"direction":direction,"phase_deg":phase,
            "variables":v,"project_variant":PROJECT_VARIANT,"readback_si":rb,"modes":modes,"seconds":time.time()-start,
            "label":"HFSS_SIMULATED / NOT_MEASURED", "branch_identity":"REQUIRES_EZ_REVIEW"}
    write_json(out/(run_id+".json"),result)
    return result


def export_modes(path, results):
    fields=["run_id","state","direction","phase_deg","mode","f_ghz","q","branch_identity"]
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in results:
            for m in r["modes"]:
                w.writerow({k: (m[k] if k in m else r[k]) for k in fields})


def verify(h,out,cores):
    results=[]; table=[]
    for (state,direction),(tf,tq) in TARGETS.items():
        v=nominal();v["Rs"],v["Xs"]=state_values(state,v)
        v["$ph_x"],v["$ph_y"]=(90,0) if direction=="across" else (0,90)
        r=solve_row(h,out,"VERIFY_"+state+"_"+direction,v,cores);results.append(r)
        # Match by both frequency and Q; do not silently substitute the high-Q FR4 slab branch.
        candidates=[m for m in r["modes"] if abs(m["f_ghz"]/tf-1)<=.01 and abs(m["q"]/tq-1)<=.30]
        best=min(r["modes"],key=lambda m:(abs(m["f_ghz"]/tf-1)/.01)**2+(abs(m["q"]/tq-1)/.30)**2)
        row={"state":state,"direction":direction,"target_f_ghz":tf,"target_q":tq,
             "f_ghz":best["f_ghz"],"q":best["q"],"mode":best["mode"],
             "f_error_pct":100*abs(best["f_ghz"]/tf-1),"q_error_pct":100*abs(best["q"]/tq-1),
             "status":"PASS" if len(candidates)==1 else "FAIL" if not candidates else "AMBIGUOUS_MODE_HARD_STOP"}
        table.append(row)
        print("VERIFY_ROW "+json.dumps(row),flush=True)
    gate={"status":"PASS" if all(r["status"]=="PASS" for r in table) else "HARD_STOP_VERIFY",
          "aedt":h.desktop_class.aedt_version_id,"reference_version":"2025.1","reference":"USER_REPORTED_HFSS",
          "table":table,"note":"Frequency/Q replication gate only; field branch identity still requires Ez review."}
    write_json(out/"verify.json",gate); export_modes(out/"VERIFY_modes.csv",results)
    print("VERIFY_RESULT "+json.dumps(gate),flush=True)
    if gate["status"]!="PASS":
        raise RuntimeError("HARD_STOP_VERIFY: no screening allowed")
    return gate


def artefact_cases(h,out,cores,label,air_below):
    results=[]
    for state in ("OFF","ON"):
        for direction in ("across","along"):
            v=nominal();v["Rs"],v["Xs"]=state_values(state,v)
            # solve_row() re-applies every variable, so the floor depth under test MUST be set here;
            # otherwise the nominal 30 mm silently overwrites a 60 mm build (bug found 2026-11-10 21:15).
            v["air_below"]=float(air_below)
            v["$ph_x"],v["$ph_y"]=(90,0) if direction=="across" else (0,90)
            results.append(solve_row(h,out,f"ARTEFACT_{label}_{state}_{direction}",v,cores))
    export_modes(out/f"ARTEFACT_{label}_modes.csv",results)
    return results


LIVE = {"h": None}   # the currently live app handle (replaced on reconnect) so main() always releases the right one


def sweep(h,out,path,cores,reconnect=None):
    """Solve every row of a CSV. A failed row is retried ONCE after a full AEDT reconnect; a second failure is logged to
    FAILED_ROWS.json and the row is SKIPPED - one bad row never stops the sweep (user rule, 2026-11-10)."""
    rows=list(read_rows(path)); results=[]; failed=[]
    def fresh():
        nonlocal h
        try:
            h.release_desktop()
        except Exception:
            pass
        h=reconnect(); LIVE["h"]=h
    for run_id,v in rows:
        cache=out/(run_id+".json")
        if cache.exists():
            r=json.loads(cache.read_text())
            if r["variables"]!=v:
                raise RuntimeError("Run ID collision with different variables")
            results.append(r); continue
        errors=[]
        for attempt in (1,2):
            try:
                r=solve_row(h,out,run_id,v,cores)
                results.append(r); break
            except Exception as exc:
                errors.append(f"attempt {attempt}: {type(exc).__name__}: {str(exc)[:300]}")
                print(f"ROW_FAILED {run_id} {errors[-1]}",flush=True)
                if reconnect is None:
                    break
                fresh()          # clean AEDT state before the retry (and before the next row if this was the last attempt)
        if not any(r_["run_id"]==run_id for r_ in results):
            failed.append({"run_id":run_id,"attempts":len(errors),"errors":errors,"variables":v})
            write_json(out/"FAILED_ROWS.json",failed)
            print(f"ROW_SKIPPED_AS_FAILED {run_id}",flush=True)
        export_modes(out/(path.stem+"_modes.csv"),results)
    return results


C_LIGHT = 299792458.0
SLOPE_FLAT, SLOPE_SW = 0.005, 0.02          # GHz/deg (user, 2026-11-10 criterion fix)
TRACK_DF_GHZ, TRACK_QRATIO = 0.6, 8.0


def _track(m90, others):
    """Nearest mode at a neighbouring phase (|df|<=0.6 GHz, Q within 8x). Returns (mode|None, ambiguous)."""
    cand = sorted((m for m in others if m["q"] > 0 and 1/TRACK_QRATIO <= m["q"]/m90["q"] <= TRACK_QRATIO
                   and abs(m["f_ghz"]-m90["f_ghz"]) <= TRACK_DF_GHZ), key=lambda m: abs(m["f_ghz"]-m90["f_ghz"]))
    if not cand:
        return None, False
    d0 = abs(cand[0]["f_ghz"]-m90["f_ghz"])
    amb = len(cand) > 1 and abs(cand[1]["f_ghz"]-m90["f_ghz"]) < 1.5*d0 + 1e-9
    return cand[0], amb


def beta_over_k0(phis, fs, period_mm, f_target=5.8):
    """beta/k0 at f_target from a 3-point branch f(phi) (piecewise linear). Returns (value, 'BRACKETED'|'EXTRAPOLATED')."""
    pts = sorted(zip(phis, fs))
    phi_star, how = None, "BRACKETED"
    for (p0, f0), (p1, f1) in zip(pts, pts[1:]):
        if (f0-f_target)*(f1-f_target) <= 0 and f1 != f0:
            phi_star = p0 + (f_target-f0)*(p1-p0)/(f1-f0); break
    if phi_star is None:
        (p0, f0), (p1, f1) = (pts[0], pts[1]) if abs(pts[0][1]-f_target) < abs(pts[-1][1]-f_target) else (pts[-2], pts[-1])
        if f1 == f0:
            return None, "NO_DISPERSION"
        phi_star = p0 + (f_target-f0)*(p1-p0)/(f1-f0); how = "EXTRAPOLATED"
    beta = math.radians(phi_star)/(period_mm*1e-3)
    return beta/(2*math.pi*f_target*1e9/C_LIGHT), how


def classify_group(rs, windows):
    """Per state/direction: track every low-Q mode across phase 80/90/100 and classify by |df/dphi|."""
    by = {}
    for r in rs:
        by.setdefault((r["state"], r["direction"]), {})[int(round(r["phase_deg"]))] = r
    out = {}; period = rs[0]["variables"]["period"]
    for (state, direction), ph in sorted(by.items()):
        key = f"{state}_{direction}"
        if not all(k in ph for k in (80, 90, 100)):
            out[key] = {"status": "AWAITING_80_100", "phases_have": sorted(ph)}; continue
        modes = []
        for m in ph[90]["modes"]:
            if not (m["q"] < 1000 and 4.0 <= m["f_ghz"] <= 7.0) or excluded_mode(m, windows):
                continue
            m80, a80 = _track(m, ph[80]["modes"]); m100, a100 = _track(m, ph[100]["modes"])
            e = {"mode90": m["mode"], "f90_ghz": round(m["f_ghz"], 4), "q90": round(m["q"], 1)}
            if m80 is None or m100 is None:
                e.update(klass="MIXED", note="UNTRACKED across 80/100 (treated as MIXED, conservative)")
            else:
                slope = (m100["f_ghz"]-m80["f_ghz"])/20.0
                curv = abs((m100["f_ghz"]-m["f_ghz"])-(m["f_ghz"]-m80["f_ghz"]))/10.0   # GHz/deg change of slope
                a = abs(slope)
                klass = "RESONANCE" if a < SLOPE_FLAT else "SURFACE_WAVE" if a > SLOPE_SW else "MIXED"
                e.update(f80_ghz=round(m80["f_ghz"], 4), f100_ghz=round(m100["f_ghz"], 4), slope_ghz_per_deg=round(slope, 5),
                         klass=klass, ambiguous_track=bool(a80 or a100), curvature_ghz_per_deg=round(curv, 5))
                if a80 or a100:
                    e["note"] = "AMBIGUOUS_TRACK (another mode within 1.5x the matching distance)"
                if klass == "SURFACE_WAVE":
                    b, how = beta_over_k0([80, 90, 100], [m80["f_ghz"], m["f_ghz"], m100["f_ghz"]], period)
                    e.update(beta_over_k0_5p8=None if b is None else round(b, 4), beta_status=how)
            e["dist_from_5p8_ghz"] = round(abs(m["f_ghz"]-5.8), 3)
            e["counts_for_1ghz_criterion"] = e["klass"] in ("RESONANCE", "MIXED")
            e["violates_1ghz"] = e["counts_for_1ghz_criterion"] and abs(m["f_ghz"]-5.8) < 1.0
            modes.append(e)
        out[key] = {"status": "CLASSIFIED", "modes": modes}
    return out


def anisotropy_report(cls):
    """OFF and ON |beta_across - beta_along|/k0 from BRACKETED surface-wave branches (nearest to 5.8 GHz per direction)."""
    rep = {}
    for state in sorted({k.split("_")[0] for k in cls}):
        beta = {}
        for d in ("across", "along"):
            c = cls.get(f"{state}_{d}", {})
            sw = [m for m in c.get("modes", []) if m["klass"] == "SURFACE_WAVE" and m.get("beta_status") == "BRACKETED"
                  and m.get("beta_over_k0_5p8") is not None]
            if sw:
                beta[d] = min(sw, key=lambda m: m["dist_from_5p8_ghz"])["beta_over_k0_5p8"]
        rep[state] = {"beta_over_k0_by_direction": beta,
                      "anisotropy": round(abs(beta["across"]-beta["along"]), 4) if len(beta) == 2 else None}
    if rep.get("OFF", {}).get("anisotropy") and rep.get("ON", {}).get("anisotropy") is not None:
        rep["ON_over_OFF_anisotropy_ratio"] = round(rep["ON"]["anisotropy"]/rep["OFF"]["anisotropy"], 3)
    return rep


ABSORBER_ONLY = {"f_ghz": (6.69, 6.74), "q_max": 60.0}   # Q~21 mode at 6.712 GHz seen only with the 377-ohm floor


def exclusion_windows():
    """Modes excluded from SCREEN: the absorber-only mode, plus any frequency window the artefact test flagged ARTEFACT."""
    wins = [{"why": "absorber-only Q~21 mode", "f": ABSORBER_ONLY["f_ghz"], "q_max": ABSORBER_ONLY["q_max"]}]
    path = HERE / "ARTEFACT_TEST.csv"
    if path.exists():
        with path.open(newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["verdict"] == "ARTEFACT":
                    # one +-1% window per reported frequency (not one span covering all three)
                    for k in ("pec30_f_ghz", "pec60_f_ghz", "z377_f_ghz"):
                        fq = float(r[k])
                        wins.append({"why": f"ARTEFACT per artefact test ({r['state']} {r['direction']} m{r['pec30_mode']})",
                                     "f": (fq*0.99, fq*1.01), "q_max": None})
    return wins


def excluded_mode(m, windows):
    return any(w["f"][0] <= m["f_ghz"] <= w["f"][1] and (w["q_max"] is None or m["q"] <= w["q_max"]) for w in windows)


def analyse(out):
    windows=exclusion_windows()
    records=[]; dropped_waive=[]
    for path in sorted(out.glob("*.json")):
        r=json.loads(path.read_text())
        if "modes" in r and "FABWAIVE" in r["run_id"]:
            dropped_waive.append(r["run_id"]); continue      # user dropped all FABWAIVE rows (anti-pad 0.4 mm not fabricable)
        if "modes" in r and not r["run_id"].startswith("VERIFY_"):
            ok, why = readback_ok(r)
            if not ok:
                raise RuntimeError(f"ANALYSIS REFUSED row {r['run_id']}: {why}")
            records.append(r)
    # Screening is phase 90 only: beta at 5.8 and flatness are NOT identifiable.
    # Never match branches solely by frequency order through crossings.
    ranked=[]
    groups={}
    for r in records:
        geom={k:v for k,v in r["variables"].items() if k not in {"$ph_x","$ph_y","Rs","Xs"}}
        geom["_project_variant"]=r.get("project_variant","full")
        key=json.dumps(geom,sort_keys=True)
        groups.setdefault(key,[]).append(r)
    for key,rs in groups.items():
        cls=classify_group(rs, windows)
        awaiting=[k for k,v in cls.items() if v["status"]!="CLASSIFIED"]
        violating=[(k,m) for k,v in cls.items() if v["status"]=="CLASSIFIED" for m in v["modes"] if m["violates_1ghz"]]
        low=violating
        nearest=min((m["dist_from_5p8_ghz"] for k,m in violating),default=None)
        # Bias-path (via + stub/trace) modes: tracked as every low-Q (<1000) mode in 4.0-7.0 GHz per state/direction, lowest first.
        # Identity was established on the BAR50 cell by no-stub / no-via solves (see handoff); here it is a frequency-ordered list,
        # not a field-verified branch label. Criterion applies to them too: none within 1 GHz of 5.8 GHz in OFF and ON.
        bias={}
        for r in rs:
            lows=sorted(((m["f_ghz"],m["q"],m["mode"]) for m in r["modes"]
                         if m["q"]<1000 and 4.0<=m["f_ghz"]<=7.0 and not excluded_mode(m, windows)))
            bias[f"{r['state']}_{r['direction']}"]=[{"mode":i,"f_ghz":round(f,4),"q":round(q,1),
                                                    "dist_from_5p8_ghz":round(abs(f-5.8),3),"inside_1ghz_band":abs(f-5.8)<1.0}
                                                   for f,q,i in lows]
        ranked.append({"geometry":json.loads(key),"runs":[r["run_id"] for r in rs],"low_q_modes_by_state":bias,
            "floor_uncertainty_note":"cell-mode frequencies carry +-1 % floor uncertainty (artefact test)",
            "classification":cls,"anisotropy":anisotropy_report(cls),
            "awaiting_80_100":awaiting,"violations_1ghz":[{"case":k,**{x:m[x] for x in ("f90_ghz","q90","klass","dist_from_5p8_ghz")}} for k,m in violating],
            "criterion_1ghz":("PENDING_PHASES" if awaiting else "FAIL" if violating else "PASS"),
            "low_q_nearby_mode_count":len(low),"nearest_low_q_distance_ghz":nearest,
            "flatness":"NEEDS_FULL_PHASE_SWEEP_AND_EZ_BRANCH_LABELS",
            "beta_off_on":"NOT_AVAILABLE_FROM_CURRENT_DATA",
            "surface_wave_q":"NEEDS_EZ_BRANCH_IDENTIFICATION",
            "status":"PROVISIONAL_SCREEN_ONLY_NO_DESIGN_PASS"})
    ranked.sort(key=lambda r:(r["criterion_1ghz"]=="PENDING_PHASES",r["low_q_nearby_mode_count"],-(r["nearest_low_q_distance_ghz"] or 0)))
    # Same geometry across project variants (full stub / no stub / periodic trace): lowest low-Q mode per state, for the bias-path comparison.
    cmp={}
    for e in ranked:
        g=dict(e["geometry"]); var=g.pop("_project_variant")
        gk=json.dumps(g,sort_keys=True)
        cmp.setdefault(gk,{})[var]={k:(v[0]["f_ghz"] if v else None) for k,v in e["low_q_modes_by_state"].items()}
    variant_comparison=[{"geometry":json.loads(k),"lowest_low_q_mode_ghz_by_variant":v} for k,v in cmp.items() if len(v)>1]
    # KEY RESULT (user): spread of each bias-path (low-Q) mode across stub 0 / 0.3 / 1 / 2 mm and the periodic trace, same cell geometry.
    # Stub axis label: nostub project -> "stub0", trace project -> "trace", full project -> "stub<stub_l>". Modes are compared by RANK
    # (1st, 2nd lowest low-Q mode per state/direction) - a frequency ordering, not a field-verified branch label.
    axis={}
    for e in ranked:
        g=dict(e["geometry"]); var=g.pop("_project_variant")
        lab="stub0" if var=="nostub" else "trace" if var=="trace" else f"stub{g['stub_l']:g}"
        g.pop("stub_l",None)
        axis.setdefault(json.dumps(g,sort_keys=True),{})[lab]=e["low_q_modes_by_state"]
    spread=[]
    for gk,labs in axis.items():
        if len(labs)<2: continue
        row={"geometry":json.loads(gk),"labels":sorted(labs),"by_state":{}}
        for sd in sorted({k for v in labs.values() for k in v}):
            per={}
            for rank in (1,2):
                vals={lab:v[sd][rank-1]["f_ghz"] for lab,v in labs.items() if sd in v and len(v[sd])>=rank}
                if len(vals)>=2:
                    sp=round(max(vals.values())-min(vals.values()),4)
                    per[f"mode_rank{rank}"]={"f_ghz_by_label":vals,"spread_ghz":sp,"FLAG_spread_gt_0p3_ghz":sp>0.3}
            row["by_state"][sd]=per
        spread.append(row)
    report={"status":"NO_SOLUTIONS" if not records else "PROVISIONAL_SCREEN_RANKING",
       "criteria":{"frequency_exclusion_ghz":[4.8,6.8],"low_q_threshold":1000,
           "flat_span_85_to_95_ghz":.20,"surface_wave_q_min":1000,
           "off_anisotropy_min":.09,"on_to_off_anisotropy_max":.2,"fab_clearance_mm_min":.15},
       "excluded_mode_windows":windows,"variant_comparison":variant_comparison,"bias_path_stub_spread":spread,
       "dropped_fabwaive_rows":len(dropped_waive),"ranking":ranked,"warning":"8 modes from 4.5 GHz cannot prove absence of modes below 4.5 or above the highest exported mode. FULL requires top-surface Ez branch labels and 5.8 GHz bracketing; extrapolation forbidden."}
    write_json(out/"analysis.json",report)
    print(f"ANALYSIS {report['status']}: {len(records)} rows, {len(ranked)} geometry groups, "
          f"{len(spread)} stub-spread sets, flagged spreads: "
          f"{sum(v2.get('FLAG_spread_gt_0p3_ghz',False) for r_ in spread for st in r_['by_state'].values() for v2 in st.values())}, "
          f"dropped FABWAIVE rows: {len(dropped_waive)} -> {out/'analysis.json'}",flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    mode=p.add_mutually_exclusive_group(required=True)
    for flag in ["build","verify","analyse","all","artefact"]:
        mode.add_argument("--"+flag,action="store_true")
    mode.add_argument("--sweep",type=Path)
    mode.add_argument("--diag",type=Path,help="Diagnostic rows CSV (mode identification); no SCREEN gates, never SCREEN data")
    p.add_argument("--project",type=Path,default=HERE/"work"/"CellParam.aedt")
    p.add_argument("--output",type=Path,default=HERE/"runs")
    p.add_argument("--floor",choices=["pec","impedance377"],default="pec")
    p.add_argument("--air-below",type=float,choices=[30.0,60.0],default=30.0)
    p.add_argument("--floor-z",type=float,default=None,
                   help="Exact global floor z in mm; calibration VERIFY uses -33.21")
    p.add_argument("--variant",choices=["full","novia","nostub","trace"],default="full",help="Diagnostic geometry variants (BUILD only)")
    p.add_argument("--no-mesh-ops",action="store_true",help="BUILD without any manual mesh operations (matches a plain user setup)")
    p.add_argument("--cores",type=int,default=None,help="HPC cores requested per solve (default: all logical CPUs)")
    p.add_argument("--row-timeout-min",type=float,default=0,help="Hard wall-clock limit per solved row; 0 = none")
    p.add_argument("--diode",choices=sorted(PROFILES),default=None,
                   help="Required: bar50 = user calibration cell, bar90 = BAR90-02LRH")
    p.add_argument("--verify-gate",type=Path,default=None,
                   help="verify.json of the BAR50 calibration run required before SCREEN")
    p.add_argument("--artefact-gate",type=Path,default=HERE/"ARTEFACT_TEST.json",
                   help="Artefact-test summary that must exist before SCREEN")
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if a.diode is None:
        p.error("--diode bar50|bar90 is required")
    if a.analyse:
        set_profile(a.diode); analyse(a.output);return
    if (a.verify or a.all) and a.diode!="bar50":
        p.error("VERIFY targets are the user's BAR50 calibration; use --diode bar50 for --verify")
    if a.all and a.diode!="bar50":
        p.error("--all is the BAR50 calibration flow only; run SCREEN separately with --sweep")
    set_profile(a.diode)
    global ROW_TIMEOUT_S
    ROW_TIMEOUT_S=a.row_timeout_min*60
    cores=a.cores or os.cpu_count()
    global PROJECT_VARIANT
    for cand in (a.output/"build_summary.json",):
        if cand.exists():
            PROJECT_VARIANT=json.loads(cand.read_text()).get("variant","full")
    if a.build and not a.diag:
        PROJECT_VARIANT=a.variant
    if a.sweep:
        # All SCREEN gates are checked BEFORE AEDT is launched.
        gate=a.verify_gate or a.output/"verify.json"
        if not gate.exists() or json.loads(gate.read_text())["status"]!="PASS":
            p.error(f"VERIFY must pass before SCREEN (no PASS verify.json at {gate})")
        if not a.artefact_gate.exists():
            p.error(f"Artefact test must be completed and reported before SCREEN (missing {a.artefact_gate})")
        built=a.output/"build_summary.json"
        info=json.loads(built.read_text()) if built.exists() else {}
        if info.get("floor_boundary")!="IMPEDANCE_377_OHM_PER_SQUARE" or info.get("diode_profile")!=a.diode:
            p.error("SCREEN needs a project built with the absorbing floor and this --diode profile "
                    f"(build_summary says {info.get('floor_boundary')}/{info.get('diode_profile')})")
        if not a.project.exists():
            p.error(f"SCREEN project not found: {a.project}")
    a.project.parent.mkdir(parents=True,exist_ok=True)
    if (a.build or a.all) and a.project.exists():
        raise RuntimeError("Refusing to overwrite existing project; use a fresh --project")
    # Decide BEFORE connecting: connect() itself creates the project file, so testing existence afterwards
    # wrongly skipped the build and sent an empty design to the solver.
    need_build = a.build or a.all or (a.artefact and not a.project.exists())
    h=None
    try:
        h=connect(a.project); LIVE["h"]=h
        write_json(a.output/"environment.json",{"python":sys.version,"executable":sys.executable,
            "pyaedt":importlib.metadata.version("pyaedt"),"aedt":h.desktop_class.aedt_version_id,
            "logical_cpus":os.cpu_count(),"platform":platform.platform(),
            "requested_cores":cores,"no_mesh_ops":a.no_mesh_ops,"row_timeout_min":a.row_timeout_min,"non_graphical":True,
            "license_features":"See LICENSE_DIAGNOSIS.md and lmstat output; not inferred from CPU count"})
        if need_build:
            build(h,a.output,a.floor,a.air_below,a.floor_z,not a.no_mesh_ops,a.variant)
        if a.verify or a.all:
            verify(h,a.output,cores)
        if a.artefact:
            label=("PEC" if a.floor=="pec" else "Z377")+f"_{int(a.air_below)}mm"
            artefact_cases(h,a.output,cores,label,a.air_below)
        if a.sweep:
            sweep(h,a.output,a.sweep,cores,reconnect=lambda: connect(a.project))
        if a.diag:
            sweep(h,a.output,a.diag,cores,reconnect=lambda: connect(a.project))
        if a.all:
            analyse(a.output)
    except Exception as exc:
        write_json(a.output/"STOP.json",{"status":"STOPPED", "error":str(exc),"time_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())})
        raise
    finally:
        h=LIVE["h"] or h
        if h is not None:
            # Preserve the original error even if AEDT has already terminated.
            try:
                h.save_project()
            except Exception:
                pass
            try:
                h.release_desktop()
            except Exception:
                pass


if __name__ == "__main__":
    main()
