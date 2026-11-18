"""Audit D3p6 rotated-supercell mode candidates from saved HFSS rows."""

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"
CELL = "B2p2_D3p6"
PHASES = (15, 30, 45)
P_SUPER = 0.0928
P_CELL = 0.0116
C0 = 299792458.0
TRACKS = {
    "x_low_mean_candidate": (-1, 1, (3, 3, 2)),
    "x_low_across_alternative": (-1, 1, (4, 4, 4)),
    "x_high_candidate": (-1, 2, (8, 8, 7)),
    "x_high_45_alternative": (-1, 2, (8, 8, 8)),
    "positive_fold_candidate": (+1, 1, (7, 7, 8)),
    "flat_candidate": (+1, 1, (5, 5, 5)),
}


def row(phase):
    path = R / f"SUPER_{CELL}_ph{phase}_PML40_resume4_modes.json"
    return json.loads(path.read_text(encoding="utf-8"))


def mag_profile(phase, mode):
    path = R / f"SUPER_{CELL}_ph{phase}_PML40_resume4_mode{mode}_MagE.fld"
    if not path.exists():
        return None
    samples = np.loadtxt(path, skiprows=2)[:, 3].reshape(8, 2)
    profile = np.sqrt(np.mean(samples * samples, axis=1))
    return profile / np.linalg.norm(profile)


def main():
    modes = {p: row(p)["modes"] for p in PHASES}
    baseline = json.loads((R / "folded_qbase.json").read_text(encoding="utf-8"))["rows"]
    out = []
    for name, (sign, branch, numbers) in TRACKS.items():
        selected = [modes[p][n - 1] for p, n in zip(PHASES, numbers)]
        freq = [m["f_ghz"] for m in selected]
        vg = [abs((freq[1] - freq[0]) * 1e9 / 15) * 360 * P_SUPER,
              abs((freq[2] - freq[0]) * 1e9 / 30) * 360 * P_SUPER,
              abs((freq[2] - freq[1]) * 1e9 / 15) * 360 * P_SUPER]
        field = [mag_profile(p, n) for p, n in zip(PHASES, numbers)]
        correlation = [None if field[i] is None or field[i+1] is None else
                       float(field[i] @ field[i+1]) for i in (0, 1)]
        for i, (p, n, mode) in enumerate(zip(PHASES, numbers, selected)):
            base = next(x for x in baseline if x["cell"] == CELL
                        and x["super_phase_deg"] == p and x["fold_sign"] == sign
                        and x["branch_index"] == branch)
            f = mode["f_ghz"]
            beta = base["per_cell_phase_deg"] * C0 / (360 * P_CELL * f * 1e9)
            harmonic = beta - C0 / (P_SUPER * f * 1e9)
            angle = math.degrees(math.asin(harmonic)) if abs(harmonic) <= 1 else None
            qinv = 1 / mode["Q"] - 1 / base["Q_base"]
            alpha = math.pi * f * 1e9 * qinv / vg[i] if vg[i] > 0 else None
            out.append({"cell": CELL, "track": name, "fold_sign": sign,
                        "branch_index": branch, "phase_deg": p, "mode": n,
                        "f_ghz": f, "Q": mode["Q"], "Q_base": base["Q_base"],
                        "per_cell_phase_deg": base["per_cell_phase_deg"],
                        "beta_over_k0": beta, "vg_m_s": vg[i],
                        "alpha_np_m": alpha, "radiated_harmonic_angle_deg": angle,
                        "field_cosine_prev": correlation[i-1] if i else None,
                        "field_cosine_next": correlation[i] if i < 2 else None,
                        "max_solved_tets": {15:102994, 30:100780, 45:118236}[p],
                        "mode_span_reaches_6p5": row(p)["mode_span_reaches_6p5"]})
    (R / "d3p6_resume4_tracks.json").write_text(json.dumps(out, indent=2) + "\n",
                                                encoding="utf-8")
    with (R / "d3p6_resume4_tracks.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out[0].keys())
        writer.writeheader()
        writer.writerows(out)
    for line in out:
        print(line["track"], line["phase_deg"], line["mode"],
              f"f={line['f_ghz']:.4f}", f"Q={line['Q']:.1f}",
              f"vg={line['vg_m_s']/1e6:.2f}", f"alpha={line['alpha_np_m']:.3f}",
              f"angle={line['radiated_harmonic_angle_deg']:.1f}")


if __name__ == "__main__":
    main()
