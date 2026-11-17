"""Summarise the completed MODEL_ONLY 5.8 GHz meander screen."""

import json
import math
from pathlib import Path


def fmt(value, digits=4):
    return "—" if value is None else f"{value:.{digits}f}"


def main(folder):
    folder = Path(folder)
    rows = json.loads((folder / "screen.json").read_text())["rows"]
    if any("ok" not in row for row in rows):
        raise ValueError("incomplete screen")
    good = [row for row in rows if row["ok"]]
    straight = next((row for row in good if row["geometry"]["kind"] == "straight"), None)
    if straight is None:
        raise ValueError("straight-slot comparator has no valid model result")
    ranked = sorted((row for row in good if row["geometry"]["kind"] == "meander"
                     and abs(row["geometry"]["orientation"]) < 1e-12),
                    key=lambda r: (-r["M_t"], -r["radiating_potential_to_dielectric_ratio"]))
    finalists, seen = [], set()
    for row in ranked:
        g = row["geometry"]
        key = (g["tooth_count"], g["tooth_width"], g["slot_width"])
        if key in seen:
            continue  # rotation swaps x/y; do not call it another geometry finalist
        seen.add(key)
        finalists.append({"name": row["name"], "label": "MODEL_ONLY / NOT_HFSS_VALIDATED",
                         "period_mm": row["pitch_mm"], "disc_radius_mm": g["diameter"] * 500,
                         "tooth_count": g["tooth_count"], "tooth_width_mm": g["tooth_width"] * 1000,
                         "slot_width_mm": g["slot_width"] * 1000, "tooth_depth_mm": g["slot_span"] * 1000,
                         "orientation_deg": math.degrees(g["orientation"]),
                         "beta_over_k0_along": row["beta_over_k0_along"],
                         "beta_over_k0_across": row["beta_over_k0_across"],
                         "M_t": row["M_t"], "alpha_radiating_potential_np_per_m": row["alpha_radiating_potential_np_per_m"],
                         "radiating_potential_to_dielectric_ratio": row["radiating_potential_to_dielectric_ratio"]})
        if len(finalists) == 3:
            break
    (folder / "finalists.json").write_text(json.dumps({"label": "MODEL_ONLY / NOT_HFSS_VALIDATED",
        "selection": "M_t then radiating-potential/dielectric ratio; rotational duplicates removed; HFSS must verify",
        "cell_parameters": finalists}, indent=2) + "\n", encoding="utf-8")
    lines = ["# 5.8 GHz meander cell screen (MODEL_ONLY)", "",
             "RT5880 εr 2.2, tanδ 0.0009, 3.175 mm over PEC; period 11.6 mm, disc radius 5.6 mm. The straight slot has a 3.0 × 0.4 mm cut. The proposed baseline has 3 turns, a 1.2 mm copper finger, 0.4 mm slot, and 3.0 mm full horizontal run along x. This full run is how `tooth_depth` was interpreted. MoM uses a uniform rooftop mesh with at least three cells across the narrowest gap (88 cells/axis for 0.4 mm; 0.2 mm would require at least 174 and was not solved). The β values use the normal-incidence tensor-sheet approximation, not a Bloch eigenmode. No HFSS cell solve has checked these rankings.", "",
             "| Cell | β/k₀ x | β/k₀ y | Δβ vs straight x/y | M_t | α rad potential (Np/m) | α diel mean (Np/m) | rad/diel | Mesh |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        if not row["ok"]:
            lines.append(f"| {row['name']} | FAILED: {row['error']} | | | | | | | |")
            continue
        dx = row["beta_over_k0_along"] - straight["beta_over_k0_along"]
        dy = row["beta_over_k0_across"] - straight["beta_over_k0_across"]
        diel = 0.5 * (row["alpha_dielectric_along_np_per_m"] + row["alpha_dielectric_across_np_per_m"])
        lines.append(f"| {row['name']} | {fmt(row['beta_over_k0_along'])} | {fmt(row['beta_over_k0_across'])} | {dx:+.4f}/{dy:+.4f} | {fmt(row['M_t'])} | {fmt(row['alpha_radiating_potential_np_per_m'])} | {fmt(diel)} | {fmt(row['radiating_potential_to_dielectric_ratio'], 2)} | {row['mesh']['nx']}² |")
    lines += ["", "M_t = |X_x - X_y|/|X_x + X_y|, where X is the effective surface reactance. The α rad value is a first-order **30° rotating-axis modulation potential** inferred from the local tensor. For the infinite uniform bound-mode cell, the modeled radiating leakage is zero; this nonzero potential requires spatial hologram modulation. The dielectric α comes from the lossy-slab transverse-resonance root; copper is PEC, so conductor loss is omitted. The ratio is a potential, not a measured radiation efficiency. The large four-tooth α is particularly sensitive to the first-order approximation. At p=11.6 mm, spatial dispersion and three-cell-wide gaps need HFSS and mesh-convergence checks.", "",
              "The removed-tooth island row is geometry only: no diode, bias or switching claim. Relative to the static baseline, its β/k₀ shifts by -0.0031 along x and +0.0001 across y in this model; this does not establish a radiating diode state. [Finalists](finalists.json) list three unique static parameter sets, selected by M_t then model loss ratio. The [ideal aperture estimate](../aperture/REPORT.md) is separate and does not use these cell responses.", "",
              "## HANDOFF", "", "- STATUS: Cell screen complete, MODEL_ONLY.",
              "- Cell JSON: `finalists.json` (2–3 static cells); removed-tooth row is not a finalist.",
              "- β/k₀ along/across, α, M_t, loss ratio: see table and per-cell JSON.",
              "- HFSS or MODEL_ONLY: MODEL_ONLY; HFSS cell and array checks pending.",
              "- Next request: HFSS eigenmode/periodic-cell check of the finalist β and loss, then a physical CP/leakage board study if the cell passes.", ""]
    (folder / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
