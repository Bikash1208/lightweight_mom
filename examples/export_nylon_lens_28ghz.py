"""Part D: export the designed 28 GHz nylon dual-beam lens (STEP + STL + layout + README).

    .venv/bin/python examples/export_nylon_lens_28ghz.py

Reads results/transmitarray_dual30/{layout,summary,local_periodicity}.json (scripts/research/phase19_design.py,
phase19_local_periodicity.py).  Writes hfss_exports/nylon_lens_28ghz_dual30/.  STEP/STL are gitignored (regenerable).
"""

import csv
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.transmitarray import export as ex  # noqa: E402

RES = ROOT / "results" / "transmitarray_dual30"
OUT = ROOT / "hfss_exports" / "nylon_lens_28ghz_dual30"


def main():
    layout = json.loads((RES / "layout.json").read_text())
    summ = json.loads((RES / "summary.json").read_text())
    inp = summ["inputs"]
    cells = layout["cells"]
    d, p, plate, total = inp["d_flats_mm"], inp["pitch_mm"], inp["plate_mm"], inp["total_mm"]
    OUT.mkdir(parents=True, exist_ok=True)
    body = ex.build_lens(cells, d, p, plate, total)
    ex.save_step_named(body, OUT / "lens.step")
    ex.save_stl(body, OUT / "lens.stl")
    ver = ex.verify_step(OUT / "lens.step", cells, d, p, plate, total)
    ver.pop("outline_vertices")
    (OUT / "verification.json").write_text(json.dumps(ver, indent=1))
    print(json.dumps(ver, indent=1))

    keys = [
        "x_mm",
        "y_mm",
        "kind",
        "radius_mm",
        "theta_deg",
        "phi_deg",
        "target_deg",
        "residual_deg",
        "power_transmission",
    ]
    with open(OUT / "layout.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for c in cells:
            w.writerow([f"{c[k]:.5f}" if isinstance(c[k], float) else c[k] for k in keys])
    (OUT / "layout.json").write_text(
        json.dumps(dict(inputs=inp, cells=[{k: c[k] for k in keys} for c in cells]), indent=1)
    )

    pat = summ[summ["chosen"]]["pattern"]
    res = summ[summ["chosen"]]["residual_deg"]
    lp = RES / "local_periodicity.json"
    lp_txt = ""
    if lp.exists():
        s = json.loads(lp.read_text())
        a = s["adjacency"]
        lp_txt = (
            f"- Pillar|hole adjacencies in this layout: {a['counts']['hole-pillar']} of {a['total_internal_edges']} internal cell edges "
            f"({100 * a['pillar_hole_fraction']:.1f} %); median target-phase step across them {a['pillar_hole_target_phase_step_deg']['median']:.0f} deg.\n"
            f"- **Local-periodicity caveat** (1x4 strips, order 0, normal incidence; `results/transmitarray_dual30/local_periodicity.json`): the model "
            f"treats each cell as if it sat in an infinite array of itself. Against exact supercell RCWA the error vector, relative to the cell amplitude, "
            f"is ~0.07 at a 60 deg neighbour step, ~0.2 at 110 deg and 0.4-0.5 at 150-175 deg, at pillar|hole and same-family boundaries alike; "
            f"two dense-hole|thin-pillar pairs were worse (up to {s['summary']['mixed_order0_error_over_cell_amp_max']:.2f}). "
            f"Expect the HFSS pattern to depart from the numbers above (lower gain, higher side lobes) -- they are an optimistic, model-only bound.\n"
        )
    (OUT / "README.md").write_text(
        f"""# 28 GHz nylon dual-beam (+-30 deg) transmitarray lens -- model-only, awaiting your HFSS

Files: `lens.step` (ONE fused solid named **Lens**, mm), `lens.stl` (printing), `layout.csv` / `layout.json` (per-cell x, y, kind, octagon
circumradius, incidence angle, target phase, residual, power transmission), `verification.json` (STEP re-import checks).
Regenerate: `examples/export_nylon_lens_28ghz.py` (STEP/STL are gitignored).

## Geometry
Regular-octagon outline, {d:g} mm across flats (flats on the x and y axes), {p:g} mm cell pitch, {total:g} mm total thickness,
{plate:g} mm base plate; z = 0 is the **exit** face, z = {total:g} the **feed-side** face. {sum(c["kind"] == "pillar" for c in cells)} pillar
cells (octagonal nylon pillar on the plate) and {sum(c["kind"] == "hole" for c in cells)} hole cells (solid nylon with an octagonal
through-hole cut through plate and lens). Min feature {inp.get("d_seam_deg") and 0.6} mm. One material: nylon, eps_r {inp["eps_r"]}, tan_d {inp["tan_delta"]}.

## Horn
Phase centre on the axis at **z = +{inp["focal_mm"]:g} mm** (F = 1.2 D) above the feed-side face, pointing at the lens, **x-polarised**
(E along x). Model horn: cos^q pattern, q = {inp["q"]:.2f}, {inp["edge_taper_db"]:g} dB edge taper at the lens edge. The beams point to
theta = +-{inp["beam_deg"]:g} deg in the x-z plane (phi = 0 and 180 deg) on the exit side (z < 0).

## Model predictions to compare with HFSS (all model-only: HFSS comparison NOT_AVAILABLE_FROM_CURRENT_DATA)
- Peak realised gain per beam (this model, horn power = 1, includes spillover, dissipation, model radiation efficiency): **{pat["beam_phi0"]["peak_gain_dbi"]:.1f} dBi**
  at theta = {pat["beam_phi0"]["theta_deg"]:.2f} deg (phi 0) and {pat["beam_phi180"]["theta_deg"]:.2f} deg (phi 180); level difference **{pat["level_difference_db"]:.2f} dB**
  (zero by the layout's mirror symmetry, not an independent measurement).
- Co-pol side-lobe level relative to the beams: **{pat["sll_db"]:.1f} dB**; max cross-pol {pat["max_cross_dbi"]:.1f} dBi (off-plane).
- Aperture-averaged power transmission of the cells: {pat["aperture_transmission_mean"]:.3f} ({10 * __import__("math").log10(pat["aperture_transmission_mean"]):.2f} dB);
  fraction of horn power reaching the exit aperture: {pat["aperture_power_fraction_of_horn"]:.3f}; radiated into propagating waves: {pat["radiated_power_fraction_of_horn"]:.3f}.
- Per-cell phase residual vs requested: max {res["max_abs"]:.2f} deg, rms {res["rms"]:.3f} deg (cells above 1 deg: {res["n_over_1deg"]}).
{lp_txt}
## Suggested HFSS setup
Material nylon (eps_r {inp["eps_r"]}, tan_d {inp["tan_delta"]}); horn at the phase centre above; radiation (or PML) boundary >= lambda/4 (2.7 mm)
beyond the structure on all sides; far-field sphere, Ludwig-3 co/cross for x polarisation; fine mesh in the nylon (lambda/(10 n) ~ 0.7 mm).
Because the cell responses are local-periodicity predictions, expect the pillar|hole boundaries and the lens edge to be the first
places a full-wave run departs from the model. The lens has {summ[summ["chosen"]]["n_cells"]} cells; a symmetry plane (x and y) can be used.
"""
    )
    print("wrote", OUT)


if __name__ == "__main__":
    main()
