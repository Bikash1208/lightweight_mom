"""Finite-radius launch pad: FD (`sidequest_launcher_fd`) for loaded annuli rho >= rho_in; launch estimated RELATIVE to the uniformly loaded FD case times the spectral uniform value (the FD energy balance is not closed, absolute FD launch is not used).
Baseline check: RT5880 with rho_in = 6.5 mm (cells from 6.5 mm) against the HFSS-based 0.627.  MODEL_ONLY.  -> results/sidequest_hybrid/k2_pad_study.json"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sidequest_launcher_fd as FD

F = 10.5e9
CASES = {"RT5880_1.575": (2.2, 0.0009, 1.575e-3, 7.38e-3, 0.846 * 0.627 / 0.611), "RO3006_1.27": (6.15, 0.002, 1.27e-3, 6.9e-3, 0.96)}   # (er, tand, h, Ybar, spectral uniform launch, scaled to the HFSS baseline)
RIN = (0.000635, 0.0026, 0.0035, 0.0045, 0.0065)


def job(a):
    name, rin = a
    e, t, h, ys, _l = CASES[name]
    r = FD.solve(F, e, t, h, 6e-3, ys, rin)
    return {"stack": name, "rho_in_mm": rin * 1e3, **{k: (complex(v) if isinstance(v, complex) else float(v)) for k, v in r.items()}}


if __name__ == "__main__":
    jobs = [(n, r) for n in CASES for r in RIN]
    with ProcessPoolExecutor(len(jobs)) as ex:
        res = list(ex.map(job, jobs))
    out = []
    for n in CASES:
        rows = [x for x in res if x["stack"] == n]
        uni = rows[0]
        for x in rows:
            x["launch_rel"] = CASES[n][4] * (x["P_sw"] / uni["P_sw"]) / (x["P_in"] / uni["P_in"])
            x["Z_in"] = [x["Z_in"].real, x["Z_in"].imag]
            out.append(x)
            print(n, "rho_in %.2f mm: launch_rel %.3f  R_in %.1f X_in %.1f  P_sw %.1f P_in %.1f  resid %.2f" % (x["rho_in_mm"], x["launch_rel"], x["R_in"], x["Z_in"][1], x["P_sw"], x["P_in"], x["fit_resid"]), flush=True)
    Path("results/sidequest_hybrid/k2_pad_study.json").write_text(json.dumps(out, indent=1), encoding="utf-8", newline="\n")
