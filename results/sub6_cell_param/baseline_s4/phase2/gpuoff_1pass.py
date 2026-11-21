"""GPU-off one-pass solve of a FRESH copy of the stored BOTH archive (.aedt sha256 6af57096...). min=max=1, 16 cores, gpus=0 explicit. Saves and exports S11 mag/phase + Z11 before release."""
import os, json, re, time, zipfile, hashlib, traceback
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"; os.environ["ANSYSLMD_LICENSE_FILE"] = "1055@10.100.105.13"
W = Path(r"C:\work\s4_work\gpuoff_1pass"); W.mkdir(exist_ok=True, parents=True)
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"result.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    arch = Path(r"C:\work\s4_work\both_full\BOARD_R150_V31_BOTH.aedtz"); assert arch.stat().st_size == 677593
    dst = W/"copy"; dst.mkdir(exist_ok=True)
    with zipfile.ZipFile(arch) as z: data = z.read("$PROJECTDIR/BOARD_R150_V31_BOTH.aedt")
    (dst/"BOARD_R150_V31_BOTH.aedt").write_bytes(data); res["aedt_sha256"] = hashlib.sha256(data).hexdigest(); assert res["aedt_sha256"].startswith("6af57096"); save()
    h = Hfss(project=str(dst/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id
    assert [s.name for s in h.setups] == ["Setup1"]; st = h.setups[0]; assert str(st.props["DrivenSolverType"]).lower().startswith("direct")
    st.props["MinimumPasses"] = 1; st.props["MaximumPasses"] = 1; assert st.update()
    res["setup"] = {k: str(h.setups[0].props.get(k)) for k in ("Frequency", "MinimumPasses", "MaximumPasses", "MaxDeltaS", "DrivenSolverType")}; save()
    t0 = time.monotonic(); ok = h.analyze_setup("Setup1", cores=16, gpus=0, use_auto_settings=False, blocking=True)
    res["wall_s"] = round(time.monotonic()-t0, 1); res["analyze_ok"] = bool(ok); save()
    h.save_project(); res["saved"] = True; res["results_dir_exists"] = (dst/"BOARD_R150_V31_BOTH.aedtresults").exists(); save()
    prof = Path(h.export_profile("Setup1", output_file=str(W/"solve.prof"))); txt = prof.read_text(errors="replace")
    res["tets"] = re.findall(r"Tetrahedra:\s*(\d+)", txt)[:8]; res["max_mem_gb"] = re.findall(r"Max memory/process:\s*([\d.]+)\s*GB", txt); res["solve_types"] = sorted(set(re.findall(r"Type:\s*(\w+), Cores", txt))); res["gpu_lines"] = [l.strip() for l in txt.splitlines() if "GPU" in l]; save()
    acf = list(dst.rglob("pyaedt_config.acf")); res["acf_NumGPUs"] = re.findall(r"NumGPUs=(\d+)", acf[0].read_text()) if acf else None
    def val(e):
        try:
            d = h.post.get_solution_data(expressions=e, setup_sweep_name="Setup1 : LastAdaptive")
            return float(d.get_expression_data(expression=e, formula="real")[1][0]) if d else None
        except Exception as ex: return "ERR " + str(ex)[:80]
    for e in ("dB(S(Port1,Port1))", "mag(S(Port1,Port1))", "ang_deg(S(Port1,Port1))", "re(S(Port1,Port1))", "im(S(Port1,Port1))", "re(Z(Port1,Port1))", "im(Z(Port1,Port1))"): res[e] = val(e)
    save(); h.save_project()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.save_project(); h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass

