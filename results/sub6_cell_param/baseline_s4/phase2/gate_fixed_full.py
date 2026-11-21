"""FULL BOTH gate with the feed fix (c1 mesh ops + c2 redrawn port annulus; no other geometry change) on a COPY, Direct, 16 cores, min 6 / max 12 / dS 0.02 as stored.
Saves the project, exports S11, Z11 and the far field (FFD + stored-format metrics) BEFORE releasing the desktop."""
import os, sys, json, math, re, time, shutil, traceback
from pathlib import Path
import numpy as np
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
sys.path.insert(0, r"C:\work\lwmom_s4\scripts\research"); sys.dont_write_bytecode = True
from full_board_gate_v31_extract import ffd, solution_point, direction_separation
W = Path(r"C:\work\s4_work\gate_fixed"); W.mkdir(exist_ok=True, parents=True)
if not (W/"BOARD_R150_V31_BOTH.aedt").exists(): shutil.copy(r"C:\work\s4_work\fix_c2\BOARD_R150_V31_BOTH.aedt", W/"BOARD_R150_V31_BOTH.aedt")
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {"cores": 16}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; m = h.modeler
    res["mesh_ops"] = [o.name for o in h.mesh.meshoperations]; res["port_objects"] = [(b.name, str(b.props.get("Objects"))) for b in h.boundaries if "Port" in b.name]
    assert sorted(res["mesh_ops"]) == ["ML_HOLE", "ML_PIN"]  # ML_PORT was deleted with the old PORT sheet in fix c2 and "PORT_NEW" in m.object_names and "PORT" not in m.object_names, res
    st = h.setups[0]; assert st.name == "Setup1" and str(st.props["DrivenSolverType"]).lower().startswith("direct")
    st.props["MaximumPasses"] = 12; st.props["MinimumPasses"] = 6; st.props["MaxDeltaS"] = 0.02; assert st.update()  # max first: min 6 > max 1 is rejected
    res["setup"] = {k: str(h.setups[0].props.get(k)) for k in ("Frequency", "MinimumPasses", "MaximumPasses", "MaxDeltaS", "DrivenSolverType", "BasisOrder")}; save()
    h.save_project()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=16, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project(); res["saved"] = True; res["results_dir_exists"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["delta_s"] = re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", txt); res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt)
    res["max_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt); res["converged"] = "Adaptive Passes converged" in txt; save()
    def val(expr): return solution_point(h, expr)
    for e in ("dB(S(Port1,Port1))", "mag(S(Port1,Port1))", "ang_deg(S(Port1,Port1))", "re(Z(Port1,Port1))", "im(Z(Port1,Port1))"): res[e] = val(e)
    res["accepted_power_fraction"] = 1 - 10 ** (res["dB(S(Port1,Port1))"] / 10); save()
    # far field in the stored metrics format (same logic as full_board_gate_v31_extract.py)
    target = W/"farfield"; target.mkdir(exist_ok=True)
    assert h.export_antenna_metadata(frequencies=[5.68e9], setup="Setup1 : LastAdaptive", sphere="Global", variations={}, output_dir=str(target), export_element_pattern=True, export_objects=False, export_touchstone=True, export_power=True), "FFD_EXPORT_FAILED"
    files = list(target.rglob("*.ffd")); assert len(files) == 1; res["ffd_path"] = str(files[0])
    theta, phi, et, ep = ffd(files[0]); res["grid"] = {"theta": [theta[0], theta[-1], len(theta)], "phi": [phi[0], phi[-1], len(phi)]}
    plus = np.abs(et + 1j*ep)**2; minus = np.abs(et - 1j*ep)**2
    ix_th = int(np.argmin(abs(theta-30))); ix_ph = int(np.argmin(abs(phi-0)))
    hfss_diff = solution_point(h, "dB(GainRHCP)", 30, 0) - solution_point(h, "dB(GainLHCP)", 30, 0)
    ffd_diff = 10*math.log10(max(plus[ix_th, ix_ph], 1e-300)/max(minus[ix_th, ix_ph], 1e-300))
    rh = plus if abs(ffd_diff-hfss_diff) < abs(-ffd_diff-hfss_diff) else minus; lh = minus if rh is plus else plus
    res["ffd_rhcp_convention"] = "et_plus_j_ep" if rh is plus else "et_minus_j_ep"
    for hand, inten in (("RHCP", rh), ("LHCP", lh)):
        vis = inten[theta <= 90, :-1]; k = np.unravel_index(int(np.argmax(vis)), vis.shape)
        pth, pph = float(theta[k[0]]), float(phi[k[1]]); pg = solution_point(h, f"dB(Gain{hand})", pth, pph)
        sc = pg - 10*math.log10(max(float(inten[k]), 1e-300)); grid = 10*np.log10(np.maximum(inten, 1e-300)) + sc
        sep = direction_separation(theta, phi, pth, pph); sll = float(np.max(grid[(sep >= 10) & (theta[:, None] <= 90)])) - pg
        res[hand] = {"peak_gain_db": pg, "peak_theta_deg": pth, "peak_phi_deg": pph, "sll_db_10deg_exclusion": sll,
                     "at_30_phi0_db": solution_point(h, f"dB(Gain{hand})", 30, 0), "at_30_phi180_db": solution_point(h, f"dB(Gain{hand})", 30, 180),
                     "at_36_phi180_db": solution_point(h, f"dB(Gain{hand})", 36, 180)}
    save(); h.save_project(); res["results_dir_exists_final"] = (W/"BOARD_R150_V31_BOTH.aedtresults").exists(); res["complete"] = True; save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass


