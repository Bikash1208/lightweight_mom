import os, json, shutil, math
from pathlib import Path
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
W = Path(r"C:\work\s4_work\sma_t3_zo"); W.mkdir(exist_ok=True)
src = Path(r"C:\work\s4_work\sma_t3\run")
if not (W/"BOARD_R150_V31_BOTH_SMA.aedt").exists():
    shutil.copy(src/"BOARD_R150_V31_BOTH_SMA.aedt", W/"BOARD_R150_V31_BOTH_SMA.aedt"); shutil.copytree(src/"BOARD_R150_V31_BOTH_SMA.aedtresults", W/"BOARD_R150_V31_BOTH_SMA.aedtresults")
temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH_SMA.aedt"), design="SMA_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["setups"] = [(s.name, s.is_solved) for s in h.setups]; port = list(h.excitation_names)[0]; res["port"] = port; res["sweeps"] = h.existing_analysis_sweeps
    sw = "Setup1 : LastAdaptive"
    def val(e):
        try:
            d = h.post.get_solution_data(expressions=e, setup_sweep_name=sw)
            return float(d.get_expression_data(expression=e, formula="real")[1][0]) if d else None
        except Exception as ex: return "ERR " + str(ex)[:80]
    P = port
    for e in (f"dB(S({P},{P}))", f"mag(S({P},{P}))", f"ang_deg(S({P},{P}))", f"re(S({P},{P}))", f"im(S({P},{P}))", f"re(Z0({P}))", f"im(Z0({P}))", f"re(Zo({P}))", f"re(Z({P},{P}))", f"im(Z({P},{P}))", f"re(Zin({P}))", f"re(PortZ0({P}))", f"re(Gamma({P}))", f"im(Gamma({P}))"): res[e] = val(e)
    try: res["port_props"] = {k: str(v)[:200] for k, v in [b for b in h.boundaries if b.name == "Port1"][0].props.items()}
    except Exception as ex: res["port_props_err"] = str(ex)[:100]
    (W/"zo.json").write_text(json.dumps(res, indent=1, default=str))
finally:
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
