"""Optional HFSS spot-check set: 6 single nylon cells across the phase range (Phase 19 Part B).

    .venv/bin/python examples/export_nylon_cell_checkset_28ghz.py

Reads results/nylon_cells_28ghz/library.csv; writes hfss_exports/nylon_cells_28ghz/ (STEP per cell,
rcwa_reference.csv, README.md).  The STEP files are gitignored; everything else is committed.
"""

import csv
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dielectric_cell.fourier import octagon_vertices  # noqa: E402
from src.dielectric_cell.library import best_window, order_by_fill  # noqa: E402

LIB = ROOT / "results" / "nylon_cells_28ghz" / "library.csv"
OUT = ROOT / "hfss_exports" / "nylon_cells_28ghz"
P_MM, H_MM, PLATE_MM = 5.0, 25.0, 1.0


def cell_solid(kind, radius_mm):
    import cadquery as cq

    base = cq.Workplane("XY").box(P_MM, P_MM, H_MM, centered=(True, True, False))
    if radius_mm == 0:
        body = base if kind == "hole" else cq.Workplane("XY").box(P_MM, P_MM, PLATE_MM, centered=(True, True, False))
        return body
    pts = [tuple(round(float(c), 6) for c in v) for v in octagon_vertices(radius_mm, 22.5)]
    prism = cq.Workplane("XY").workplane(offset=PLATE_MM).polyline(pts).close().extrude(H_MM - PLATE_MM)
    if kind == "pillar":
        plate = cq.Workplane("XY").box(P_MM, P_MM, PLATE_MM, centered=(True, True, False))
        return plate.union(prism)
    return base.cut(prism)


def main():
    with open(LIB) as fh:
        rows = list(csv.DictReader(fh))
    tab, ref = [], {}
    for r in rows:
        key = (r["kind"], round(float(r["radius_mm"]), 6))
        ref.setdefault(key, []).append(r)
        if float(r["theta_deg"]) == 0:
            t = complex(float(r["t_pp_re"]), float(r["t_pp_im"]))
            tab.append(
                (r["kind"], float(r["radius_mm"]), float(r["fill"]), math.degrees(np.angle(t)), 20 * math.log10(abs(t)))
            )
    ordered = order_by_fill(tab)
    i0, i1, *_ = best_window(ordered)
    win = ordered[i0 : i1 + 1]
    targets = np.linspace(win[0][3], win[-1][3], 6)
    picks = [min(win, key=lambda x, tg=tg: abs(x[3] - tg)) for tg in targets]
    OUT.mkdir(parents=True, exist_ok=True)
    import cadquery as cq

    with open(OUT / "rcwa_reference.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cell", "kind", "radius_mm", "theta_deg", "phi_deg", "t_pp_re", "t_pp_im", "t_ss_re", "t_ss_im"])
        for n, (kind, rad, *_x) in enumerate(picks, 1):
            name = f"cell{n}_{kind}_R{rad:.3f}"
            cq.exporters.export(cell_solid(kind, rad), str(OUT / f"{name}.step"))
            for r in ref[(kind, round(rad, 6))]:
                w.writerow(
                    [
                        name,
                        kind,
                        r["radius_mm"],
                        r["theta_deg"],
                        r["phi_deg"],
                        r["t_pp_re"],
                        r["t_pp_im"],
                        r["t_ss_re"],
                        r["t_ss_im"],
                    ]
                )
    (OUT / "README.md").write_text(
        f"""# Nylon cell spot-check set, 28 GHz (optional)

Six single cells spanning the library's 360-degree window (pitch {P_MM} mm, total thickness {H_MM} mm, {PLATE_MM} mm base
plate, flats on the x/y axes).  One STEP per cell, z = 0 is the bottom (exit) face, z = {H_MM} the top (feed-side) face.
`rcwa_reference.csv` has the model's complex transmission (project convention exp(+jwt)) at the sampled angles.

HFSS setup: nylon eps_r 2.45, tan_d 0.0185 (single material); unit-cell boundaries (master/slave, 5 mm square) on the
four side faces, Floquet ports above and below, scan angle theta/phi from the CSV.  S21 phase must be referenced to
the **top face of the cell for the incident wave and the bottom face for the transmitted wave** (port de-embedding of
the air above/below, no extra offset).  `t_pp` is p (in-plane E) polarisation, `t_ss` is s; at theta = 0, p = x.

Targets: phase within +-3.6 deg and |S21| within +-0.1 dB of the CSV.  HFSS comparison status:
NOT_AVAILABLE_FROM_CURRENT_DATA (model-only until you run it).

| cell | kind | radius [mm] | RCWA phase [deg] | RCWA |S21| [dB] (theta=0, p) |
|---|---|---|---|---|
"""
        + "\n".join(f"| {n} | {k} | {rd:.3f} | {ph:.1f} | {ls:.2f} |" for n, (k, rd, _f, ph, ls) in enumerate(picks, 1))
        + "\n"
    )
    print(*picks, sep="\n")


if __name__ == "__main__":
    main()
