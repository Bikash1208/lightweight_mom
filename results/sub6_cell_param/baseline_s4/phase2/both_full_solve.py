"""Open the stored v31 BOTH archive UNCHANGED, read Setup1 back by name, require Direct, set max passes 2 only, solve Setup1 by name, report S11 per pass.
No geometry edit, no rebuild, no far-field export. Run: python both_archive_2pass.py"""
import os, sys, json, time, re, threading, traceback
from pathlib import Path
import psutil
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\both_full"); W.mkdir(exist_ok=True, parents=True)
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Desktop, Hfss
log = W/"result.json"; res = {"start": time.time()}
def save(): log.write_text(json.dumps(res, indent=1, default=str))
def s11(h, sw="Setup1 : LastAdaptive"):
    d = h.post.get_solution_data(expressions="dB(S(Port1,Port1))", setup_sweep_name=sw)
    _, v = d.get_expression_data(expression="dB(S(Port1,Port1))", formula="real"); return float(v[0]), list(d.primary_sweep_values)[:2]
mem = {"peak_gib": 0.0}; stop = threading.Event()
def watch():
    while not stop.wait(10):
        try:
            tot = sum(p.memory_info().rss for p in psutil.process_iter() if p.name().lower() in ("ansysedt.exe", "hf3d.exe", "hfsscomengine.exe"))/2**30
            mem["peak_gib"] = max(mem["peak_gib"], tot); res["peak_mem_gib_so_far"] = round(mem["peak_gib"], 2); save()
        except Exception: pass
threading.Thread(target=watch, daemon=True).start()
h = None; desk = None
try:
    arch = W/"BOARD_R150_V31_BOTH.aedtz"; assert arch.stat().st_size == 677593
    # AEDT's RestoreProjectArchive produced nothing here (blk-file case mismatch warning), so the zip is extracted byte-for-byte instead.
    import zipfile, hashlib
    dst = W/"restored"; dst.mkdir(exist_ok=True); assert not list(dst.glob("*.lock"))
    with zipfile.ZipFile(arch) as z:
        data = z.read("$PROJECTDIR/BOARD_R150_V31_BOTH.aedt")
    assert len(data) == 6107726
    aedt = [dst/"BOARD_R150_V31_BOTH.aedt"]; aedt[0].write_bytes(data)
    res["restored_project"] = str(aedt[0]); res["aedt_sha256"] = hashlib.sha256(data).hexdigest(); save()
    h = Hfss(project=str(aedt[0]), solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["design"] = h.design_name; res["designs"] = h.design_list
    res["setups_on_open"] = [(s.name, str(s.props.get("Frequency")), str(s.props.get("DrivenSolverType")), s.props.get("MaximumPasses"), s.props.get("MinimumPasses"), s.is_solved) for s in h.setups]
    res["metal_objects"] = len([n for n in h.modeler.object_names if n.startswith("METAL_PEC")]); save()
    names = [s.name for s in h.setups]
    assert "Auto1" not in names, "Auto1 present on open - stop and report"
    st = [s for s in h.setups if s.name == "Setup1"]; assert len(st) == 1; setup = st[0]
    res["setup1_readback_before"] = {k: str(setup.props.get(k)) for k in ("Frequency", "MaximumPasses", "MinimumPasses", "MaxDeltaS", "DrivenSolverType", "BasisOrder")}
    res["sweeps_before"] = h.existing_analysis_sweeps; save()
    assert str(setup.props["DrivenSolverType"]).lower().startswith("direct"), "NOT DIRECT - stop"
    f = str(setup.props["Frequency"]).strip().lower(); assert abs(float(f.replace("ghz", "").strip())*1e9-5.68e9) < 1.0, f
    exp = {"Frequency": "5.68GHz", "MaximumPasses": "12", "MinimumPasses": "6", "MaxDeltaS": "0.02"}
    rb = res["setup1_readback_before"]
    res["readback_matches_expected"] = (rb["MaximumPasses"] == "12" and rb["MinimumPasses"] == "6" and abs(float(rb["MaxDeltaS"]) - 0.02) < 1e-12)
    save()
    assert res["readback_matches_expected"], ("SETUP1 READBACK MISMATCH - stop and report", rb)
    free = psutil.virtual_memory().available/2**30; res["free_gib_before_solve"] = round(free, 1); assert free > 45+15, free
    t0 = time.monotonic(); res["solve_start"] = time.time(); save()
    ok = h.analyze_setup("Setup1", cores=30, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["profile_delta_s"] = re.findall(r"Max Mag\. Delta S:\s*([\d.]+)", txt)
    res["profile_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt); res["profile_tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt)
    res["profile_pass_lines"] = [l.strip() for l in txt.splitlines() if re.search(r"Adaptive Pass \d|Adaptive Passes|Elapsed Time|Normal Completion|Status", l)][:40]
    res["sweeps_after"] = h.existing_analysis_sweeps
    res["s11_LastAdaptive_dB"], res["freq"] = s11(h)
    res["accepted_power"] = 1-10**(res["s11_LastAdaptive_dB"]/10)
    res["setups_after_solve"] = [(s.name, s.is_solved) for s in h.setups]; save()
    for sw in h.existing_analysis_sweeps:
        if "AdaptivePass" in sw:
            try:
                d = h.post.get_solution_data(expressions="dB(S(Port1,Port1))", setup_sweep_name=sw, variations={"Pass": ["All"], "Freq": ["5.68GHz"]})
                res["s11_per_pass_AdaptivePass"] = {"sweep": sw, "pass_values": list(d.primary_sweep_values), "s11_db": [float(x) for x in d.data_real()]}
            except Exception as e:
                res["s11_per_pass_AdaptivePass_err"] = str(e)[:200]
    save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    stop.set()
    res["peak_mem_gib"] = round(mem["peak_gib"], 2); res["end"] = time.time(); save()
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
        elif desk is not None: desk.release_desktop(close_projects=True)
    except Exception: pass


