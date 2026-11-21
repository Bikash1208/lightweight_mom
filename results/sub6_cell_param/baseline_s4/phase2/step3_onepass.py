"""Step 3: saved ONE-pass solve (min=max=1) of a COPY of the stored BOTH archive on 16 cores. Saves + exports (S11 mag/phase, Z11, |E| through the feed) BEFORE releasing; confirms results folder."""
import os, sys, json, math, re, time, traceback
from pathlib import Path
import numpy as np
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\both_step3"); temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"step3.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id
    assert [s.name for s in h.setups] == ["Setup1"]
    st = h.setups[0]; assert str(st.props["DrivenSolverType"]).lower().startswith("direct")
    st.props["MinimumPasses"] = 1; st.props["MaximumPasses"] = 1; assert st.update()
    res["setup"] = {k: str(h.setups[0].props.get(k)) for k in ("Frequency", "MinimumPasses", "MaximumPasses", "DrivenSolverType")}; save()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=16, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project(); res["saved"] = True; res["results_dir_exists"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["delta_s"] = re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", txt); res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt); res["solve_types"] = sorted(set(re.findall(r"Type:\s*(\w+), Cores", txt)))
    res["max_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt); save()
    sw = "Setup1 : LastAdaptive"
    def val(expr):
        try:
            d = h.post.get_solution_data(expressions=expr, setup_sweep_name=sw)
            return float(d.get_expression_data(expression=expr, formula="real")[1][0]) if d else None
        except Exception as e: return "ERR " + str(e)[:100]
    for e in ("dB(S(Port1,Port1))", "mag(S(Port1,Port1))", "ang_deg(S(Port1,Port1))", "re(S(Port1,Port1))", "im(S(Port1,Port1))", "re(Z(Port1,Port1))", "im(Z(Port1,Port1))", "re(Zin(Port1))", "im(Zin(Port1))", "re(Z0(Port1))", "im(Z0(Port1))"):
        res[e] = val(e)
    res["accepted_power"] = 1 - (res["mag(S(Port1,Port1))"] ** 2) if isinstance(res["mag(S(Port1,Port1))"], float) else None; save()
    def fld(name, pts):
        f = W/f"{name}.fld"
        assert h.post.export_field_file(quantity="E", solution=sw, output_file=str(f), sample_points=pts, export_in_si_system=False, intrinsics={"Freq": "5.68GHz", "Phase": "0deg"})
        with f.open() as s:
            s.readline(); s.readline(); raw = np.loadtxt(s)
        e = np.sqrt((raw[:, 3::2]**2 + raw[:, 4::2]**2).sum(axis=1)); return raw[:, :3], e
    # x-z plane through the feed (y=0), x -40..40, z -3.1..20
    xs = np.arange(-40, 40.1, 2.0); zs = [-3.0, -1.5, 1.0, 3.0, 6.0, 9.0, 12.0, 20.0]
    xyz, e = fld("xz_feed", [[x, 0.0, z] for x in xs for z in zs]); fin = np.isfinite(e)
    res["E_xz_feed_V_per_m"] = {f"z={z}": [round(float(np.nanmax(e[np.isclose(xyz[:, 2], z) & (np.abs(xyz[:, 0]) <= 4)])), 2) if (np.isclose(xyz[:, 2], z) & (np.abs(xyz[:, 0]) <= 4) & fin).any() else None,
                                          round(float(np.sqrt(np.nanmean(e[np.isclose(xyz[:, 2], z) & (np.abs(xyz[:, 0]) > 10)]**2))), 3)] for z in zs}
    res["E_xz_note"] = "per z: [max |E| within |x|<=4 mm of the feed axis, rms |E| for |x|>10 mm]"
    # coax annulus (port plane region) and board-wide at z=2.035 by radial band (rms)
    g = np.arange(-180, 181, 12.0); pts = [[x, y, 2.035] for x in g for y in g if math.hypot(x, y) <= 184]
    xyz, e = fld("z2p035", pts); r = np.hypot(xyz[:, 0], xyz[:, 1]); fin = np.isfinite(e)
    res["E_z2p035_rms_by_band"] = {f"r{a}-{b}": float(np.sqrt(np.mean(e[fin & (r >= a) & (r < b)]**2))) for a, b in ((0, 12), (12, 40), (40, 80), (80, 120), (120, 150), (150, 186))}
    xyz, e = fld("port_annulus", [[r_*math.cos(math.radians(a)), r_*math.sin(math.radians(a)), -3.175] for r_ in (0.9, 1.3, 1.8) for a in range(0, 360, 45)])
    res["E_port_annulus_mean_V_per_m"] = float(np.nanmean(e)); save()
    res["setups_solved"] = [(s.name, s.is_solved) for s in h.setups]
    res["aedt_warnings"] = [str(x) for x in h.odesktop.GetMessages("BOARD_R150_V31_BOTH", "GATE_BOTH", 1)]
    h.save_project(); res["results_dir_exists_final"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
