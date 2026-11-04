"""Run K = 2 and K = 4 grouped coded holograms: baseline, phase-offset + greedy optimisation (<= 1 dB), Pareto, full-pattern evaluation. MODEL_ONLY."""
import itertools
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np

import sidequest_coded_model as M


def eps_search(board, budget_db, grid):
    best = None
    for e in itertools.product(grid, repeat=board.K - 1):
        eps = np.array((0.0,) + e)
        bits, s0, r, a0 = M.greedy(board, eps, budget_db, 1.0, merge=False)
        sm = M.summarize(board, bits)
        key = (sm["static_frac"], -sm["lines_with_inverters"])
        if best is None or key > best[0]:
            best = (key, eps)
    return best[1]


def run(name, phis, mode, grid):
    b = M.Board(phis, mode)
    out = {"name": name, "cells": b.n, "alpha_beam0": M.A_BEAM0, "alpha_abs0": M.A_ABS0, "beta_k0": M.BETA_K0, "pareto": [], "final": {}}
    eps0 = np.zeros(b.K)
    bits0 = b.ideal_bits(eps0)
    out["baseline"] = M.summarize(b, bits0)
    print(name, "baseline", out["baseline"], flush=True)
    eps_opt = eps_search(b, 1.0, grid)
    out["eps_opt_deg"] = list(np.degrees(eps_opt))
    print(name, "eps_opt", np.round(np.degrees(eps_opt)), flush=True)
    for budget in (0.25, 0.5, 1.0):
        for fs in (1.0, 0.5, 0.25):
            for eps_lab, eps in (("eps0", eps0), ("eps_opt", eps_opt)):
                bits, s0, r, a0 = M.greedy(b, eps, budget, fs, merge=True)
                sm = M.summarize(b, bits)
                row = {"budget_db": budget, "frac_static_phase": fs, "eps": eps_lab, **sm, "loss_db_est": [float(-20 * np.log10(1 - x)) for x in r]}
                out["pareto"].append(row)
                print(name, row, flush=True)
    # full evaluation: baseline and the two 1 dB solutions with eps_opt
    res0, idl = M.evaluate(b, bits0, eps0)
    out["final"]["baseline"] = {"summary": out["baseline"], "states": res0}
    print(name, "baseline gains", [round(x["gain_dbi"], 2) for x in res0], [round(x["SLL_db"], 1) for x in res0], flush=True)
    cands = []
    for fs in (0.25, 0.5, 0.75, 1.0):
        bits, s0, r, a0 = M.greedy(b, eps_opt, 1.0, fs, merge=True)
        res, _ = M.evaluate(b, bits, eps_opt, idl)
        lab = f"fs{fs}_1dB"
        out["final"][lab] = {"summary": M.summarize(b, bits), "states": res, "loss_est_db": [float(-20 * np.log10(1 - x)) for x in r]}
        print(name, lab, out["final"][lab]["summary"], [round(x["gain_dbi"], 2) for x in res], [round(x["SLL_db"], 1) for x in res], flush=True)
        cands.append((lab, out["final"][lab]))
    base_g = min(x["gain_dbi"] for x in res0)
    base_sll = max(x["SLL_db"] for x in res0)
    ok = [(lab, c) for lab, c in cands if base_g - min(x["gain_dbi"] for x in c["states"]) <= 1.05 and max(x["SLL_db"] for x in c["states"]) <= base_sll + 3.0]
    ok.sort(key=lambda t: (-t[1]["summary"]["static_frac"], t[1]["summary"]["lines_with_inverters"]))
    out["recommended"] = ok[0][0] if ok else None
    print(name, "recommended (gain loss <= 1 dB, SLL <= baseline + 3 dB):", out["recommended"], flush=True)
    return out


if __name__ == "__main__":
    res = {}
    res["K2"] = run("K2", [0, 180], "sector", np.linspace(0, 2 * np.pi, 16, endpoint=False))
    res["K4"] = run("K4", [0, 90, 180, 270], "sector", np.linspace(0, 2 * np.pi, 6, endpoint=False))
    Path("results/sidequest_hybrid/coded_hologram_results.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
