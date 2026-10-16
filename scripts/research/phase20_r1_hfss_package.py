"""Phase 20 R1 item 4: the one-cell HFSS check package for the loaded, oblique MoM (`hfss_exports/phase20_r1_loaded_cell/`).

    python scripts/research/phase20_r1_hfss_package.py [out_dir]

STEP (Substrate, Metal = the two arms, Ground, Load_1 sheet bridging the gap), a MoM reference CSV (predictions written before any
HFSS run: HFSS lumped RLC boundaries are PARALLEL R||L||C, so only pure-R and pure-C loads, plus short and open, are listed),
the bare-slab analytic references for calibrating the TE/TM phase convention, and README with the Floquet / RLC setup.
All geometry numbers are the demo's parameters (`phase20_r1_demo.py`), not a design. Labels: MODEL_ONLY, NOT_HFSS_VALIDATED.
"""

import csv
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataclasses import replace  # noqa: E402

import phase20_r1_demo as demo  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.cad_export import export_unit_cell_step, import_named_parts  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.loads import LumpedLoad  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads, oblique_stack_s_parameters  # noqa: E402

OUT = ROOT / "hfss_exports/phase20_r1_loaded_cell"
ANGLES = (0.0, 30.0, 60.0)
PARAMS = dict(period_mm=10.0, thickness_mm=1.575, eps_r=2.2, tan_d=0.0009, freq_ghz=10.0, patch_length_mm=8.0, patch_width_mm=8.0, gap_mm=1.0)
# (label, HFSS assignment, LumpedLoad kwargs extra)
LOADS = [
    ("short", "PerfE sheet on Load_1 (or no gap: one patch)", dict(R=0.0)),
    ("R=10", "Lumped RLC on Load_1: R = 10 ohm (L, C unchecked)", dict(R=10.0)),
    ("R=50", "Lumped RLC on Load_1: R = 50 ohm", dict(R=50.0)),
    ("R=200", "Lumped RLC on Load_1: R = 200 ohm", dict(R=200.0)),
    ("C=0.5pF", "Lumped RLC on Load_1: C = 0.5 pF (R, L unchecked)", dict(C=0.5e-12)),
    ("C=0.1pF", "Lumped RLC on Load_1: C = 0.1 pF", dict(C=0.1e-12)),
    ("open", "no boundary on Load_1 (delete the sheet): empty gap", dict(open_circuit=True)),
]


MESHES = ((60, 32), (80, 40), (100, 50))  # all aligned with the 1 mm gap and the 8 mm patch edges


def solve_rows(mesh_cells, floquet):
    """List of dicts with the complex co-pol S11 of every (angle, load, mode) at one mesh."""
    cell, gap, sub = demo.build(*(PARAMS[k] for k in ("period_mm", "thickness_mm", "eps_r", "tan_d", "patch_length_mm", "patch_width_mm", "gap_mm")))
    f = PARAMS["freq_ghz"] * 1e9
    planes = ReferencePlanes(0.0, None, sub.thickness)
    rows = []
    for th in ANGLES:
        sets = [(LumpedLoad(**dict(gap, transverse_conducts=False, **kw)),) if not kw.get("open_circuit") else (LumpedLoad(**dict(gap, **kw)),) for _, _, kw in LOADS]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = compute_s_parameters_oblique_over_loads(cell, f, MeshSettings(mesh_cells, mesh_cells, floquet), planes, math.radians(th), 0.0, sets)
        for (label, hfss, _), s in zip(LOADS, res, strict=True):
            for i, m in enumerate(demo.MODES):
                rows.append({"load": label, "hfss_assignment": hfss, "theta_deg": th, "floquet_mode": f"{i + 1} ({m})", "s": complex(s.s11[i, i])})
    return rows, cell, sub


def richardson(s_a, n_a, s_b, n_b, p=2.0):
    """Extrapolate a complex quantity assumed to behave as S(n) = S_inf + c / n^p from two meshes (n_b > n_a)."""
    wa, wb = n_a**p, n_b**p
    return (wb * s_b - wa * s_a) / (wb - wa)


def _order(s60, s80, s100):
    """Observed convergence order p of S(n) = S_inf + c / n^p from three meshes (None if the sequence is not monotone enough)."""
    d1, d2 = abs(s60 - s80), abs(s80 - s100)
    if d1 < 1e-9 or d2 < 1e-9:
        return None
    from scipy.optimize import brentq

    ratio = d1 / d2
    f = lambda p: (60.0**-p - 80.0**-p) / (80.0**-p - 100.0**-p) - ratio  # noqa: E731
    try:
        return float(brentq(f, 0.2, 8.0))
    except ValueError:
        return None


def _db(z):
    return 20 * math.log10(abs(z))


def _deg(z):
    return math.degrees(np.angle(z))


def _dphase(a, b):
    return abs(((_deg(a) - _deg(b) + 180.0) % 360.0) - 180.0)


def main(out_dir=OUT):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    per_mesh = {}
    for n, fl in MESHES:
        per_mesh[n], cell, sub = solve_rows(n, fl)
    coarse, _, _ = solve_rows(20, 16)
    fine = []
    for k, r in enumerate(per_mesh[100]):
        e_hi = richardson(per_mesh[80][k]["s"], 80, r["s"], 100)
        e_lo = richardson(per_mesh[60][k]["s"], 60, per_mesh[80][k]["s"], 80)
        fine.append({
            "load": r["load"], "hfss_assignment": r["hfss_assignment"], "theta_deg": r["theta_deg"], "phi_deg": 0.0, "floquet_mode": r["floquet_mode"],
            "mom_S11_mag_dB": _db(e_hi), "mom_S11_phase_deg": _deg(e_hi),
            "mom_uncertainty_phase_deg": _dphase(e_hi, e_lo), "mom_uncertainty_mag_dB": abs(_db(e_hi) - _db(e_lo)),
            "mom_mesh100_mag_dB": _db(r["s"]), "mom_mesh100_phase_deg": _deg(r["s"]),
            "mom_mesh80_phase_deg": _deg(per_mesh[80][k]["s"]), "mom_mesh60_phase_deg": _deg(per_mesh[60][k]["s"]),
            "observed_order_p": _order(per_mesh[60][k]["s"], per_mesh[80][k]["s"], r["s"]),
            "mom_mesh20_mag_dB": _db(coarse[k]["s"]), "mom_mesh20_phase_deg": _deg(coarse[k]["s"]),
        })
    with open(out / "mom_reference.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fine[0]))
        w.writeheader()
        w.writerows(fine)
    green = LayeredGreenFunction(sub, PARAMS["freq_ghz"] * 1e9, backing="pec")
    bare = []
    for th in ANGLES:
        s0 = oblique_stack_s_parameters(green, ReferencePlanes(0.0, None, sub.thickness), math.radians(th), 0.0)
        for i, m in enumerate(demo.MODES):
            bare.append({"theta_deg": th, "floquet_mode": f"{i + 1} ({m})", "S11_mag_dB": 20 * math.log10(abs(s0.s11[i, i])), "S11_phase_deg": math.degrees(np.angle(s0.s11[i, i]))})
    with open(out / "bare_grounded_slab_reference.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bare[0]))
        w.writeheader()
        w.writerows(bare)
    gap = demo.build(*(PARAMS[k] for k in ("period_mm", "thickness_mm", "eps_r", "tan_d", "patch_length_mm", "patch_width_mm", "gap_mm")))[1]
    res = export_unit_cell_step(replace(cell, loads=(LumpedLoad(**dict(gap, R=10.0, transverse_conducts=False)),)), out / "cell.step")
    parts = import_named_parts(out / "cell.step")
    (out / "step_parts.json").write_text(json.dumps({"parts": parts, "load_bboxes_mm": res.load_bboxes_mm}, indent=1))
    p_all = [r['observed_order_p'] for r in fine if r['observed_order_p'] is not None]
    p_med, p_lo, p_hi = float(np.median(p_all)), min(p_all), max(p_all)
    ut = {m: max(r["mom_uncertainty_phase_deg"] for r in fine if r["floquet_mode"].endswith(f"({m})")) for m in demo.MODES}
    d20 = {m: max(abs(((r["mom_S11_phase_deg"] - r["mom_mesh20_phase_deg"] + 180) % 360) - 180) for r in fine if r["floquet_mode"].endswith(f"({m})")) for m in demo.MODES}
    (out / "README.md").write_text(README.format(**PARAMS, unc_te=ut["TE"], unc_tm=ut["TM"], d20_te=d20["TE"], d20_tm=d20["TM"],
                                                 unc_mag=max(r["mom_uncertainty_mag_dB"] for r in fine), p_med=p_med, p_lo=p_lo, p_hi=p_hi))
    return fine


README = """# phase20_r1_loaded_cell -- one-cell HFSS check of the loaded, oblique MoM (Phase 20 R1, item 4)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY).** The loaded-cell and the oblique-incidence MoM have no HFSS number yet; this package
lets you produce one. It is a **demonstration cell, not a design**: all dimensions are parameters of `scripts/research/phase20_r1_demo.py`.

## The cell
Period {period_mm} x {period_mm} mm; substrate eps_r {eps_r}, tan d {tan_d}, h {thickness_mm} mm; PEC ground under the substrate; one metal patch
{patch_length_mm} (x) x {patch_width_mm} (y) mm centred in the cell, split by ONE gap of {gap_mm} mm across x at the centre (two arms), bridged by one lumped
load; f = {freq_ghz} GHz. `cell.step` parts: `Substrate` (z -1.575..0), `Metal` (the two arms, 35 um slab z 0..0.035), `Ground`, `Load_1` (a single
sheet, z = 0.0175 mm, x -0.5..+0.5, y -4..+4 mm, bridging the arms' side faces).

## HFSS setup
1. Import `cell.step`. Assign **Perfect E** to `Metal` and `Ground` (the MoM metal has zero thickness: if you prefer, replace the 35 um slab by sheets at z = 0; the slab
   shifts phases by a fraction of a degree to a few degrees -- state which you used).
2. Boundaries: master/slave on the four side walls of the 10 x 10 mm cell with the **scan angle** theta = 0, 30, 60 deg, phi = 0 (plane of incidence = xz); air box on top
   (>= lambda/4 above the metal) with a **Floquet port** on its top face, modes 1 and 2 (these are TE and TM with respect to the plane of incidence: at phi = 0 TE has E along y and
   TM has E along x, i.e. across the gap). **De-embed** the port to the pattern plane z = 0 (the S11 phases in `mom_reference.csv` are referenced to z = 0).
3. **Load_1**: assign a **Lumped RLC** boundary with the current line along +x across the sheet (from x = -0.5 to +0.5 mm, y = 0). **HFSS RLC boundaries combine R, L and C in
   PARALLEL, so the reference lists only pure-R and pure-C loads**; the MoM models the series impedance R + jwL + 1/(jwC) (identical for a single element). For `short` assign Perfect E to
   `Load_1` (or fill the gap); for `open` delete the sheet (empty gap).
4. Solve at 10 GHz (adaptive, delta S <= 0.01); export S11 for Floquet modes 1 and 2 (magnitude dB and phase), at each angle and load.

## What to compare (`mom_reference.csv`)
`mom_S11_mag_dB`, `mom_S11_phase_deg`: **Richardson-extrapolated** from the 80 x 80 and 100 x 100 meshes (Floquet order 40 / 50; assumed second-order convergence, S(n) = S_inf + c/n^2; the observed order of the
60/80/100 sequence is listed per row in `observed_order_p`; median {p_med:.1f}, range {p_lo:.1f}-{p_hi:.1f}: where p < 2 the p = 2 extrapolation under-corrects and the stated uncertainty is understated); `mom_uncertainty_phase_deg` / `mom_uncertainty_mag_dB` = the difference to the same extrapolation from the 60 and 80 meshes; the raw 100 x 100 and the 20 x 20
values are listed too. MoM S are power-wave normalised, like HFSS Floquet-port S-parameters.
1. **Convention first:** run the bare grounded slab (no metal, no load; `bare_grounded_slab_reference.csv`, analytic) at each angle and mode. If HFSS's TM phase differs from the
   listed value by a constant 180 deg (a sign convention of the TM mode), apply the same shift to every TM entry before comparing; TE should agree directly. Report what you saw.
2. Then the loads. Expected qualitative result (see `results/phase20_r1_demo/README.md`): the TM (E across the gap) phase swings by ~160 deg at theta = 0 (less at larger angles) between short and open,
   and a resistive load absorbs most of the power near R ~ 10-20 ohm at theta = 0. TE is NOT load-independent: the load is a full-width (8 mm) strip, so the odd (TE) mode drives circulating currents
   through it, which absorb up to ~0.7 of the TE power at 60 deg. A narrow real diode would not carry those currents; this check compares the full-width RLC sheet with the full-width MoM strip.
3. Honest limits of the model for this comparison: zero-thickness metal; load region modelled as a conducting strip carrying Z_L along x only; no vias; rooftop MoM. This cell sits near a patch
   resonance and the rooftop discretisation error in PHASE is large: the 20 x 20 mesh is off by up to {d20_te:.0f} deg (TE) and {d20_tm:.0f} deg (TM) from the extrapolated value, the extrapolation itself is
   uncertain by up to {unc_te:.1f} deg (TE) and {unc_tm:.1f} deg (TM) (listed per row), magnitudes by up to {unc_mag:.2f} dB. Treat an HFSS-vs-MoM phase difference of that order as discretisation, not as a load-model error.
"""

if __name__ == "__main__":
    main(*(sys.argv[1:2] or []))
