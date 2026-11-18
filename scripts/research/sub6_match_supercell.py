"""List HFSS supercell mode candidates against folded unit-cell predictions."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results" / "sub6_meander_model" / "supercell_hfss"


def load_modes(cell, phase):
    stem = f"SUPER_{cell}_ph{phase}_PML40"
    for suffix in ("", "_retry1", "_resume1", "_resume4"):
        path = R / f"{stem}{suffix}_modes.json"
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))["modes"]
    raise FileNotFoundError(f"No completed mode row for {stem}")


def main(cell):
    baseline = json.loads((R / "folded_qbase.json").read_text(encoding="utf-8"))["rows"]
    phases = (15, 30, 45)
    spectra = {p: load_modes(cell, p) for p in phases}
    result = {"cell": cell, "method": "folded unrotated frequency midpoint; candidates retained",
              "modes_by_phase": spectra, "predictions": []}
    for branch in (1, 2):
        for sign in (-1, 1):
            prediction = {"branch_index": branch, "fold_sign": sign, "points": []}
            for phase in phases:
                base = next(x for x in baseline if x["cell"] == cell and
                            x["super_phase_deg"] == phase and x["fold_sign"] == sign and
                            x["branch_index"] == branch)
                f_pred = (base["f_across_ghz"] + base["f_along_ghz"]) / 2
                candidates = sorted((x for x in spectra[phase] if x["Q"] < 100000),
                                    key=lambda x: abs(x["f_ghz"] - f_pred))[:5]
                prediction["points"].append({"phase_deg": phase,
                                             "per_cell_phase_deg": base["per_cell_phase_deg"],
                                             "predicted_f_ghz": f_pred,
                                             "Q_base": base["Q_base"],
                                             "candidates": [{"mode": x["mode"], "f_ghz": x["f_ghz"],
                                                             "Q": x["Q"],
                                                             "delta_f_ghz": x["f_ghz"] - f_pred}
                                                            for x in candidates]})
            result["predictions"].append(prediction)
    out = R / f"SUPER_{cell}_mode_candidates.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for p in result["predictions"]:
        print("branch", p["branch_index"], "sign", p["fold_sign"])
        for x in p["points"]:
            print(x["phase_deg"], round(x["predicted_f_ghz"], 4),
                  [(c["mode"], round(c["f_ghz"], 4), round(c["Q"]))
                   for c in x["candidates"][:3]])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cell", choices=("ANCHOR", "B2p2_D3p0", "B2p2_D3p6"))
    main(ap.parse_args().cell)
