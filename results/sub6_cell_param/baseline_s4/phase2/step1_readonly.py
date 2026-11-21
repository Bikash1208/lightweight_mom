"""Step 1, read-only on a COPY of the converged GATE_BOTH solution: S11 phase, port Z0/gamma, |E| on planes."""
import os, json, math, traceback
from pathlib import Path
import numpy as np
os.environ["ANSYSEM_ROOT261"] = r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
W = Path(r"C:\work\s4_work\both_step1"); temp = W/"aedt_temp"; temp.mkdir(exist_ok=True); os.environ.update(TEMP=str(temp), TMP=str(temp))
from ansys.aedt.core import Hfss
res = {}; out = W/"step1.json"
def save(): out.write_text(json.dumps(res, indent=1, default=str))
h = None
try:
    h = Hfss(project=str(W/"BOARD_R150_V31_BOTH.aedt"), design="GATE_BOTH", solution_type="Modal", version="2026.1", non_graphical=True, new_desktop=True, close_on_exit=False)
    res["pid"] = h.desktop_class.aedt_process_id; res["setups"] = [(s.name, s.is_solved) for s in h.setups]
    sw = "Setup1 : LastAdaptive"
    def val(expr):
        try:
            d = h.post.get_solution_data(expressions=expr, setup_sweep_name=sw)
            if not d: return None
            return [float(x) for x in d.get_expression_data(expression=expr, formula="real")[1]][:3]
        except Exception as e: return "ERR " + str(e)[:120]
    for e in ("dB(S(Port1,Port1))", "ang_deg(S(Port1,Port1))", "mag(S(Port1,Port1))", "re(S(Port1,Port1))", "im(S(Port1,Port1))",
              "re(Z0(Port1))", "im(Z0(Port1))", "re(Zo(Port1))", "re(Gamma(Port1))", "im(Gamma(Port1))", "re(Zin(Port1))", "im(Zin(Port1))",
              "re(Z(Port1,Port1))", "im(Z(Port1,Port1))", "re(PortZO(Port1))", "mag(Y(Port1,Port1))"):
        res[e] = val(e)
    save()
    # |E| on planes: z=2.035 (above copper), z=-1.6 (mid substrate), z=14 (air above), x-z cut at y=0
    def fld(name, pts):
        f = W/f"{name}.fld"
        ok = h.post.export_field_file(quantity="E", solution=sw, output_file=str(f), sample_points=pts, export_in_si_system=False, intrinsics={"Freq": "5.68GHz", "Phase": "0deg"})
        with f.open() as s:
            s.readline(); s.readline(); raw = np.loadtxt(s)
        e2 = raw[:, 3::2]**2 + raw[:, 4::2]**2
        return raw[:, :3], np.sqrt(e2.sum(axis=1))
    g = np.arange(-180, 181, 12.0)
    for name, z in (("z2p035", 2.035), ("zmid_sub", -1.6), ("zair14", 14.0)):
        pts = [[x, y, z] for x in g for y in g if math.hypot(x, y) <= 184]
        xyz, e = fld(name, pts); r = np.hypot(xyz[:, 0], xyz[:, 1]); fin = np.isfinite(e)
        bands = {}
        for a, b in ((0, 12), (12, 40), (40, 80), (80, 120), (120, 150), (150, 156), (156, 186)):
            m = fin & (r >= a) & (r < b); bands[f"r{a}-{b}"] = [int(m.sum()), float(np.sqrt(np.mean(e[m]**2))) if m.any() else None]
        res["E_" + name] = {"rms_|E|_V/m_by_radial_band": bands, "max": float(np.nanmax(e)), "n": int(fin.sum())}
        save()
    # feed-region |E| at z=-1.6 on rings
    pts = [[r*math.cos(math.radians(a)), r*math.sin(math.radians(a)), -1.6] for r in (1, 1.5, 2.0, 3, 5, 8, 15) for a in range(0, 360, 45)]
    xyz, e = fld("feed_ring", pts); rr = np.hypot(xyz[:, 0], xyz[:, 1])
    res["E_feed_rings_mid_sub"] = {f"{r:g}mm": float(np.nanmean(e[np.isclose(rr, r)])) for r in (1, 1.5, 2.0, 3, 5, 8, 15)}
    save()
except BaseException:
    res["error"] = traceback.format_exc(); save()
finally:
    try:
        if h is not None: h.release_desktop(close_projects=True, close_desktop=True)
    except Exception: pass
