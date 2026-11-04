"""Post-process results/sidequest_hybrid/gate_cell_final_raw.json -> gate_cell_final.json (+ printed summary)."""
import json, sys, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
import sidequest_gate_cell as gc
import sidequest_gate_final as gf

raw = json.loads(Path("results/sidequest_hybrid/gate_cell_final_raw.json").read_text(encoding="utf-8"))["rows"]
Y = lambda r, i=0, j=0: complex(*r["y_ms"][i][j]) * 1e-3  # noqa: E731
sel = lambda st, mesh, f=10.54, sc=None: next((r for r in raw if r["state"] == st and r["mesh"] == mesh and abs(r["f_ghz"] - f) < 1e-6 and r["scale"] == sc), None)  # noqa: E731
out = {"cand_mm": gf.CAND, "mesh": {}, "tolerance": {}, "freq32": [], "freq64": [], "rows": {}, "angle": {}}


def modal(y, f=10.54e9, rows=3):
    r = gf.line(y, f, rows)
    b = None
    old = (gc.F0, gc.K0, gc.W0)
    gc.F0, gc.K0, gc.W0 = f, gc.free_space_wavenumber(f), gc.angular_frequency(f)
    try:
        bb = gc.beta_complex(y)
        b = None if bb is None else [float(bb.real / gc.K0), float(-bb.imag)]
    finally:
        gc.F0, gc.K0, gc.W0 = old
    return {"y_ms": [y.real * 1e3, y.imag * 1e3], "beta_over_k0_alpha_np_m": b, "strip": r}


for mesh in (32, 64):
    for st in ("OFF", "ON"):
        r = sel(st, mesh)
        if r:
            out["mesh"].setdefault(str(mesh), {})[st] = {"xx": modal(Y(r)), "yy_ms": r["y_ms"][1][1]}
for mesh in (32, 64):
    for st, key in (("ON", "R"), ("ON", "R"), ("OFF", "C"), ("OFF", "C")):
        pass
    for st, key in (("ON", "R"), ("OFF", "C")):
        for fac in (0.8, 1.2):
            r = sel(st, mesh, sc={key: fac})
            if r:
                out["tolerance"].setdefault(str(mesh), {})[f"{st}_{key}x{fac}"] = modal(Y(r))
for r in raw:
    if r["mesh"] == 32 and r["scale"] is None:
        out["freq32"].append({"state": r["state"], "f_ghz": r["f_ghz"], **modal(Y(r), r["f_ghz"] * 1e9)})
    if r["mesh"] == 64 and r["scale"] is None and abs(r["f_ghz"] - 10.54) > 1e-6:
        out["freq64"].append({"state": r["state"], "f_ghz": r["f_ghz"], **modal(Y(r), r["f_ghz"] * 1e9)})
# rows needed (mesh-64 ON if present, else mesh 32)
on = sel("ON", 64) or sel("ON", 32)
off = sel("OFF", 64) or sel("OFF", 32)
for n in (2, 3, 4, 5, 6, 8):
    out["rows"][str(n)] = {"ON": gf.line(Y(on), 10.54e9, n), "OFF": gf.line(Y(off), 10.54e9, n)}
# obliquity: Y_eff(psi) = yxx cos^2 + yyy sin^2 (diagonal tensor), psi = angle between propagation and the diode axis
for psi in (0, 15, 30, 45, 60, 75, 90):
    c, s = np.cos(np.radians(psi)), np.sin(np.radians(psi))
    yo, yf = Y(on, 0, 0) * c * c + Y(on, 1, 1) * s * s, Y(off, 0, 0) * c * c + Y(off, 1, 1) * s * s
    out["angle"][str(psi)] = {"ON_S21_3rows_db": (gf.line(yo, 10.54e9, 3) or {}).get("S21_db"), "OFF_S21_3rows_db": (gf.line(yf, 10.54e9, 3) or {}).get("S21_db")}
Path("results/sidequest_hybrid/gate_cell_final.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
print(json.dumps({k: out[k] for k in ("mesh", "tolerance", "rows", "angle")}, indent=1)[:6000])
for k in ("freq32", "freq64"):
    for r in out[k]:
        s = r["strip"]
        print(k, r["state"], r["f_ghz"], [round(v, 1) for v in r["y_ms"]], None if s is None else (round(s["S21_db"], 1), round(s["S11_db"], 1)))
