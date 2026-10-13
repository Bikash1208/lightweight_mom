"""Part D: export the designed 28 GHz nylon dual-beam lens (STEP + STL + layout + README).

    .venv/bin/python examples/export_nylon_lens_28ghz.py

Reads results/transmitarray_dual30/{layout,summary,local_periodicity}.json (scripts/research/phase19_design.py,
phase19_local_periodicity.py).  Writes hfss_exports/nylon_lens_28ghz_dual30/.  STEP/STL are gitignored (regenerable).
"""

import csv
import json
import math
import pathlib
import subprocess
import sys
import time

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
    exact, _intended = ex.export_exact_step(body, OUT / "lens.step", cells, d, p, plate, total)
    ex.save_stl(body, OUT / "lens.stl")
    ver = ex.verify_step(OUT / "lens.step", cells, d, p, plate, total)
    ver.pop("outline_vertices")
    ver["exact_coordinates"] = exact
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

    (OUT / "README.md").write_text(readme_text(summ, cells, inp, exact))
    print("wrote", OUT)


def _load(name):
    path = RES / name
    return json.loads(path.read_text()) if path.exists() else None


def _band(v, fmt="{:.2f}"):
    return f"{fmt.format(v['p50'])} ({fmt.format(v['p10'])} ... {fmt.format(v['p90'])})"


def readme_text(summ, cells, inp, exact):
    d, p, plate, total = inp["d_flats_mm"], inp["pitch_mm"], inp["plate_mm"], inp["total_mm"]
    ch = summ[summ["chosen"]]
    pat, res = ch["pattern"], ch["residual_deg"]
    lp, eb, pe = _load("local_periodicity.json"), _load("error_bar.json"), _load("pitch_estimate.json")
    rm = _load("row_models.json")
    lp_line = (
        f"{lp['adjacency']['counts']['hole-pillar']} of {lp['adjacency']['total_internal_edges']} internal edges "
        f"({100 * lp['adjacency']['pillar_hole_fraction']:.0f} %) are pillar|hole, with a median step of "
        f"{lp['adjacency']['pillar_hole_target_phase_step_deg']['median']:.0f} deg across them."
        if lp
        else ""
    )
    n_pillar = sum(c["kind"] == "pillar" and c["radius_mm"] > 0 for c in cells)
    n_empty = sum(c["kind"] == "pillar" and c["radius_mm"] == 0 for c in cells)
    n_hole = sum(c["kind"] == "hole" for c in cells)
    iso_txt = (
        f"- **ISOLATED-CELL prediction (deterministic)**: peak realised gain per beam **{pat['beam_phi0']['peak_gain_dbi']:.1f} dBi** at "
        f"theta = {eb['nominal_isolated_cell']['beam_theta_deg'][0]:.2f} deg (continuous peak search; the 0.25-deg grid of `summary.json` reports "
        f"{pat['beam_phi0']['theta_deg']:.2f}), SLL **{pat['sll_db']:.1f} dB**, max cross-pol {pat['max_cross_dbi']:.1f} dBi (off-plane). "
        "Horn power = 1; includes spillover, dissipation and the model's radiation efficiency.\n"
        if eb
        else ""
    )
    err_txt = ""
    if eb:
        sc = eb["scenarios"]
        rows = []
        for key, label in (
            ("S", "S: measured phase (primary)"),
            ("R", "R: same magnitudes, random phase"),
            ("R0", "R0: order-0 relation, random phase"),
        ):
            sm = sc[key]["summary"]
            rows.append(
                f"| {label} | {_band(sm['mean_gain_dbi'])} | {_band(sm['beam_theta_phi0_deg'])} | {_band(sm['pointing_error_deg_worst_beam'])} | "
                f"{_band(sm['sll_db'], '{:.1f}')} | {_band(sm['level_difference_db_abs'])} |"
            )
        nom = eb["nominal_isolated_cell"]
        rv = rm["windows"] if rm else []
        iso_err = [w["isolated_spectrum_error_vs_exact"] for w in rv]
        s_err = [w["S_measured_phase"]["spectrum_error_vs_exact"]["p50"] for w in rv]
        r_err = [w["R_random_phase"]["spectrum_error_vs_exact"]["p50"] for w in rv]
        pw = [w["exact_power_over_isolated"] for w in rv]
        pw_s = [w["S_measured_phase"]["power_over_isolated"]["p50"] for w in rv]
        err_txt = f"""
## Error bar on the prediction (Monte Carlo of the measured local-periodicity error)
The isolated-cell figures above treat every cell as if it sat in an infinite array of itself. Exact supercell RCWA shows each cell's effective
exit field differs from that by a perturbation that grows with the phase step to its neighbours, and this layout needs steps with a median of
{eb["edge_phase_step_deg"]["median"]:.0f} deg (90th percentile {eb["edge_phase_step_deg"]["p90"]:.0f} deg) over all {eb["edge_phase_step_deg"]["n_edges"]} neighbour edges. It is
**set by the required aperture phase gradient**, so it is propagated, not just noted. {lp_line} Method: `scripts/research/phase19_error_bar.py`; {eb["n_draws"]} draws per
scenario (seed {eb["seed"]}); percentiles are stable -- between 250 and 2000 draws the 10/50/90 values move by at most 0.04 dB (gain), 0.14 dB (SLL), 0.05 deg (pointing).
E_c' = E_c (1 + eps_c) with eps_c = K * mean(T_n/T_c - 1) over the neighbours that differ (|T_n/T_c - 1| >= {eb["d_min"]}), K bootstrapped from exact-supercell
samples with the same number of differing neighbours (1 / 2 / 4 from strips, alternating strips, 2x2 checkerboards). **Single-neighbour coupling does not add
linearly** -- the 2- and 4-neighbour supercells give a per-neighbour coupling of 0.29 and 0.15 against 0.40 for one -- so naive superposition (scenario `Cn` in
`error_bar.json`) over-predicts and is not used.

Median (10th ... 90th percentile), isolated-cell values in the first row:

| model | mean gain / beam [dBi] | beam theta, phi=0 [deg] | worst-beam pointing error [deg] | SLL [dB] | abs inter-beam level difference [dB] |
|---|---|---|---|---|---|
| ISOLATED-CELL (deterministic) | {nom["mean_gain_dbi"]:.2f} | {nom["beam_theta_deg"][0]:.2f} | 0.13 | {nom["sll_db"]:.1f} | 0.00 (identity, see below) |
{chr(10).join(rows)}

- **Phase assumption.** The supercell data DO constrain the perturbation phase (single-neighbour coupling is ~0.40 at -52 deg, not a random phase; two or four
  neighbours give a real ~0.5-0.65), so scenario S uses the measured complex K. Scenario R keeps the same magnitudes but randomises the phase: it leaves the
  coherent beam unchanged to first order and therefore shows **no** gain loss; the 4 dB between S and R is purely this phase assumption.
- **How far to trust S.** Checked on four real 9-cell windows of this layout solved exactly (`row_validation.json`, `row_models.json`): the isolated-cell order
  spectrum is off by {", ".join(f"{v:.2f}" for v in iso_err)} (relative norm, i.e. as large as the signal); model S reduces that to {", ".join(f"{v:.2f}" for v in s_err)}
  (random-phase R makes it worse: {", ".join(f"{v:.2f}" for v in r_err)}). Total propagating power: exact/isolated = {", ".join(f"{v:.2f}" for v in pw)}, model S
  {", ".join(f"{v:.2f}" for v in pw_s)} -- S under-predicts power by roughly 0.4-1.4 dB there, so **S probably overstates the gain loss**; the true systematic
  effect is a gain reduction of unknown size between roughly 0 and the S value, and HFSS is the only thing that fixes it. The window test is at normal incidence,
  x-neighbours only.
- **Inter-beam level difference.** The nominal **0.00 dB is an identity of the mirror-symmetric layout** (the layout and the horn are symmetric under x -> -x, so the two
  beams are equal by construction). It cannot fail and is **not evidence** of anything. The Monte Carlo difference above is what a perturbation does to it.
- The error bar covers the local-periodicity error only. Not included: print tolerance, nylon eps_r spread, horn pattern/phase-centre error, truncation of the
  Huygens exit-aperture model, and the lens edge clipping effects.
"""
    pitch_txt = ""
    if pe:
        a_, b_ = pe["current_p5mm"], pe["lambda_over_3"]
        r_ = pe["ratios"]
        pitch_txt = f"""
## Cell-pitch lever (estimate only -- no re-design, no layout committed)
`scripts/research/phase19_pitch_estimate.py` runs only the phase synthesis on the finer grid. p = lambda/3 needs 42 cells across (D/p must be an even integer),
so **p = 150/42 = 3.5714285... mm, not an exact decimal** (current p = {a_["pitch_mm"]:g} mm = {a_["pitch_over_lambda"]:.3f} lambda).

| | p = {a_["pitch_mm"]:g} mm (now) | p = {b_["pitch_mm"]:.4f} mm (lambda/3) |
|---|---|---|
| cells | {a_["n_cells"]} | {b_["n_cells"]} ({r_["cell_count"]:.2f}x) |
| minimum feature | 0.6 mm | 0.6 mm (unchanged); octagon R {b_["octagon_R_range_mm"][0]:.2f} ... {b_["octagon_R_range_mm"][1]:.2f} mm (was ... {a_["octagon_R_range_mm"][1]:.2f}) |
| smallest pocket (fill cap 0.804) | {a_["smallest_pocket_flat_to_flat_mm_at_fill_cap"]:.2f} mm wide, {a_["pocket_aspect_ratio_at_fill_cap"]:.0f}:1 | {b_["smallest_pocket_flat_to_flat_mm_at_fill_cap"]:.2f} mm wide, {b_["pocket_aspect_ratio_at_fill_cap"]:.0f}:1 |
| median neighbour phase step (target) | {a_["neighbour_target_phase_step_deg"]["median"]:.0f} deg | {b_["neighbour_target_phase_step_deg"]["median"]:.0f} deg ({r_["median_step"]:.2f}x) |
| 90th-percentile step; edges over 90 deg | {a_["neighbour_target_phase_step_deg"]["p90"]:.0f} deg; {100 * a_["fraction_of_edges_over_90deg"]:.0f} % | {b_["neighbour_target_phase_step_deg"]["p90"]:.0f} deg; {100 * b_["fraction_of_edges_over_90deg"]:.0f} % |
| coupling-error rms of eps (measured K, equal-amplitude targets) | {a_["eps_estimate"]["rms"]:.2f} | {b_["eps_estimate"]["rms"]:.2f} ({r_["eps_rms"]:.2f}x) |

Expected error reduction from the measured relation: about **{100 * (1 - r_["eps_rms"]):.0f} %** in the rms perturbation (the error is driven by the phase step). The cost is
2x the cells and narrower pockets (14:1); the pillar family shrinks (largest fill {b_["largest_pillar_fill"]:.2f}), which the library would have to be re-checked against.
"""
    exact_txt = ""
    if exact:
        vc = exact["vertex_check"]
        exact_txt = f"""
## STEP coordinates are exact (Phase 19 Part D-bis)
- **Why**: the OpenCASCADE writer prints reals with 12 significant digits and its boolean kernel leaves 1e-16..1e-10 noise, so the first export carried literals such as
  `8.500000000001E-02` (intended 0.085), `6.169999999999E-02` (0.0617) and `1.110223024625E-16` (0) -- 34,764 of them. The earlier "no float noise" check had a blind spot
  below 1e-13 and missed them.
- **Fix**: after the kernel writes, every `CARTESIAN_POINT`, `DIRECTION` and `VECTOR` magnitude is snapped to its intended design value when within 1e-7 mm, then re-written
  with the shortest decimal that round-trips to the same double (`4.0` -> `4.`, `0.085` -> `0.085`; never a fixed `%.Nf`). Intended values come from the layout table
  (cell centres, +-pitch/2 corners, plate and total thickness, outline, octagon vertices rounded to 1e-4 mm, outline-flat x grid intersections), never from kernel arithmetic.
- **Result on this lens**: {exact["snap"]["n_tokens"]:,} coordinates walked, {exact["inexact_literals_before"]:,} inexact before, **{exact["inexact_literals_after"]} after**; largest snap displacement
  {exact["snap"]["max_snap_mm"]:.1e} mm; no two previously distinct vertices merged ({vc["n_vertices"]:,} vertices, {vc["distinct_merge_violations"] and len(vc["distinct_merge_violations"]) or 0} violations,
  {vc["collapsed_edges"]} collapsed edges); one valid solid, self-intersection analyzer clean; volume before vs after differs by {exact["volume_rel_diff"]:.1e} (limit 1e-9).
- **Uncertainty**: `UNCERTAINTY_MEASURE_WITH_UNIT` (it lives in the STEP data section, in the geometric representation context, not the header) is set to
  **1.E-6 mm** (OpenCASCADE's default is 1e-7). It must exceed the residual coordinate noise (<= 1e-9 mm after this export) and the snap tolerance (1e-7) by an order of
  magnitude so an importer welds duplicate vertices, and stay far below the smallest real vertex separation ({vc["min_distinct_vertex_separation_mm"]:.4f} mm, a factor
  {vc["min_distinct_vertex_separation_mm"] / 1e-6:,.0f}) so it never merges distinct ones.
- **Cannot be exact in decimal (do not "fix")**: octagon outline vertices (75 tan 22.5 deg = 31.0660171779821...; stored as 31.066, a 1.7e-5 mm deliberate rounding);
  all pillar/pocket octagon vertex coordinates (cos/sin of 22.5 deg multiples times a solved radius; rounded to 1e-4 mm, up to 5e-5 mm); the solved radii themselves; the edge
  directions of the rounded octagons and the 1/sqrt(2) directions of the outline flats; intersections of octagon edges with the 45-degree outline flats. The pitch (150/30 = 5 mm)
  and the grid are exact; **a pitch such as 150/42 (lambda/3) would not be**.
- **STL**: an STL is a triangulated approximation; exactness is not a meaningful property of it, so none is claimed or chased there.
- **Short edges**: {vc["n_vertex_pairs_closer_than_0p6mm"]} vertex pairs lie closer than the 0.6 mm minimum feature (smallest {vc["min_distinct_vertex_separation_mm"]:.4f} mm), where pillar/pocket octagons are cut by the diagonal outline flats. Review these
  before printing; this export does not change the geometry.
"""
    t_commit = (
        subprocess.run(
            ["git", "show", "-s", "--format=%cd", "--date=format:%Y-%m-%d %H:%M", "20bf86d"],
            capture_output=True,
            text=True,
            cwd=ROOT,
        ).stdout.strip()
        or "n/a"
    )
    mt = lambda path: time.strftime("%Y-%m-%d %H:%M", time.localtime(path.stat().st_mtime))  # noqa: E731
    t_lib = mt(ROOT / "results" / "nylon_cells_28ghz" / "library_plate0.6mm_through.csv")
    t_layout = mt(RES / "layout.json")
    t_step = mt(OUT / "lens.step")
    return f"""# 28 GHz nylon dual-beam (+-30 deg) transmitarray lens -- model-only, awaiting your HFSS

Files: `lens.step` (ONE fused solid named **Lens**, mm), `lens.stl` (printing), `layout.csv` / `layout.json` (per-cell x, y, kind, octagon
circumradius, incidence angle, target phase, residual, power transmission), `verification.json` (STEP re-import checks).
Regenerate: `examples/export_nylon_lens_28ghz.py` (STEP/STL are gitignored).

## Geometry
Regular-octagon outline, {d:g} mm across flats (flats on the x and y axes), {p:g} mm cell pitch, {total:g} mm total thickness,
{plate:g} mm base plate; z = 0 is the **exit** face, z = {total:g} the **feed-side** face. {n_pillar} pillar cells (octagonal nylon pillar on the plate),
{n_empty} empty cells (plate only) and {n_hole} hole cells (solid nylon with an octagonal through-hole cut through plate and lens).
Minimum printable feature 0.6 mm (your stated value). One material: nylon, eps_r {inp["eps_r"]}, tan_d {inp["tan_delta"]}.

## Horn
Phase centre on the axis, F = {inp["focal_mm"]:g} mm (= 1.2 D) in front of the feed-side face, i.e. at **z = {total + inp["focal_mm"]:g} mm in the STEP coordinates**
(the feed-side face is z = {total:g}), pointing at the lens, **x-polarised** (E along x). Model horn: cos^q pattern, q = {inp["q"]:.2f}, {inp["edge_taper_db"]:g} dB edge
taper at the lens edge. The beams point to theta = +-{inp["beam_deg"]:g} deg in the x-z plane (phi = 0 and 180 deg) on the exit side (z < 0).

## Model predictions to compare with HFSS (all model-only: HFSS comparison NOT_AVAILABLE_FROM_CURRENT_DATA)
{iso_txt}- Aperture-averaged power transmission of the cells: {pat["aperture_transmission_mean"]:.3f} ({10 * math.log10(pat["aperture_transmission_mean"]):.2f} dB);
  fraction of horn power reaching the exit aperture: {pat["aperture_power_fraction_of_horn"]:.3f}; radiated into propagating waves: {pat["radiated_power_fraction_of_horn"]:.3f}.
- Per-cell phase residual vs requested: max {res["max_abs"]:.2f} deg, rms {res["rms"]:.3f} deg (cells above 1 deg: {res["n_over_1deg"]}).
- **Inter-beam level difference 0.00 dB is an IDENTITY of the mirror-symmetric layout, not a performance result.** It cannot fail and must not be quoted as evidence.
{err_txt}{pitch_txt}{exact_txt}
## Suggested HFSS setup
Material nylon (eps_r {inp["eps_r"]}, tan_d {inp["tan_delta"]}); horn at the phase centre above; radiation (or PML) boundary >= lambda/4 (2.7 mm)
beyond the structure on all sides; far-field sphere, Ludwig-3 co/cross for x polarisation; fine mesh in the nylon (lambda/(10 n) ~ 0.7 mm).
The lens has {ch["n_cells"]} cells; symmetry planes (x and y) can be used for the nominal geometry.

## Verification history
- Original Phase 19 series (Parts A-D): the fast test suite was run **once on the final tree** and not at each intermediate commit, so commit `af7b0aa` (the supercell-capable
  RCWA extension, the only follow-up commit that modified existing solver files) is verified by its outcome -- it reproduces the committed library to CSV rounding and every later
  test passed on top of it -- rather than individually. The Part D-bis commits (error bar, exact STEP coordinates) each ran the full fast suite before committing.
- Part A S21 reference plane: the "fix" was to a wrong *statement* (an empty cell was documented as S21 = 1) in the docstring and the tests, made before the first Phase 19 commit; the
  solver always returned S21 = exp(-j k0 cos(theta) h) for an empty cell. Generated after the first Phase 19 commit (20bf86d, {t_commit}) that contains the correction: cell library {t_lib},
  layout and prediction {t_layout}, this STEP {t_step} -- so the correction is not just in the docs and tests.
"""


if __name__ == "__main__":
    main()
