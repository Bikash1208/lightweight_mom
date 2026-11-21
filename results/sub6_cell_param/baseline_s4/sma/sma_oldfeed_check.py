"""SMA machine check: lab-B BOTH_oldfeed_251_1pass.aedtz, .aedt taken byte-for-byte from the zip; ONLY change min=max passes=1; saved one-pass solve on 16 cores; save + export S11 (mag/phase) and Z11 before releasing."""
import os, json, re, time, zipfile, hashlib, traceback
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\sma_oldfeed"); temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    arch = W/"BOTH_oldfeed_251_1pass.aedtz"; res["archive_sha256"] = hashlib.sha256(arch.read_bytes()).hexdigest()
    assert res["archive_sha256"] == "adbecbb255112e73a70db8e484411a234ca01c737c260415fb7f2c19f73fe6f8"
    dst = W/"copy"; dst.mkdir(exist_ok=True)
    with zipfile.ZipFile(arch) as z: data = z.read("$PROJECTDIR/BOARD_R150_GATE_V31_BOTH.aedt")
    aedt = dst/"BOARD_R150_GATE_V31_BOTH.aedt"; aedt.write_bytes(data); res["aedt_sha256"] = hashlib.sha256(data).hexdigest(); save()
    h = Hfss(project=str(aedt), solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; res["designs"] = h.design_list; res["design"] = h.design_name
    res["setups_on_open"] = [(s.name, str(s.props.get("Frequency")), str(s.props.get("DrivenSolverType")), str(s.props.get("MinimumPasses")), str(s.props.get("MaximumPasses")), s.is_solved) for s in h.setups]
    res["excitations"] = list(h.excitation_names); res["boundaries"] = [(b.name, b.type) for b in h.boundaries]; save()
    assert len(h.setups) == 1; st = h.setups[0]; res["setup_name"] = st.name
    st.props["MinimumPasses"] = 1; st.props["MaximumPasses"] = 1; assert st.update()   # the only change
    res["setup_after_edit"] = {k: str(h.setups[0].props.get(k)) for k in ("Frequency", "MinimumPasses", "MaximumPasses", "MaxDeltaS", "DrivenSolverType", "BasisOrder")}; save()
    t0 = time.monotonic(); ok = h.analyze_setup(st.name, cores=16, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project(); res["saved"] = True; res["results_dir_exists"] = any(p.name.endswith(".aedtresults") for p in dst.iterdir()); save()
    prof = Path(h.export_profile(st.name, output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt); res["max_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt)
    res["solve_types"] = sorted(set(re.findall(r"Type:\s*(\w+), Cores", txt))); res["delta_s"] = re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", txt); save()
    sw = f"{st.name} : LastAdaptive"; res["sweeps"] = h.existing_analysis_sweeps
    port = res["excitations"][0]; res["port_used"] = port
    def val(expr):
        try:
            d = h.post.get_solution_data(expressions=expr, setup_sweep_name=sw)
            return float(d.get_expression_data(expression=expr, formula="real")[1][0]) if d else None
        except Exception as e: return "ERR " + str(e)[:100]
    for e in (f"dB(S({port},{port}))", f"mag(S({port},{port}))", f"ang_deg(S({port},{port}))", f"re(S({port},{port}))", f"im(S({port},{port}))", f"re(Z({port},{port}))", f"im(Z({port},{port}))"): res[e] = val(e)
    save(); h.save_project()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass

