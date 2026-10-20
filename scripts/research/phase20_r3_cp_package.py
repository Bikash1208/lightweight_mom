"""Phase 20 R3 item 5: the HFSS check package of the CP-2 cell (`hfss_exports/phase20_r3_cp_cell/`) and the LP-interim numbers (LP as 1-bit codes of the
same cell). Predictions first; nothing here has an HFSS number yet. MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

    python scripts/research/phase20_r3_cp_package.py [out_dir]            # all (STEP, predictions, README, LP interim)
    python scripts/research/phase20_r3_cp_package.py [out_dir] --readme-only
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
sys.path.insert(0, str(ROOT / "scripts" / "research"))

from phase20_r3_finalist import build  # noqa: E402
from phase20_r3_opt import corner_kits, wrap180  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.cad_export import export_unit_cell_step, import_named_parts  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.ris.cells import madp_kit  # noqa: E402
from src.ris.polarization import cartesian_reflection, circular_reflection  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads, oblique_stack_s_parameters  # noqa: E402

OUT = ROOT / "hfss_exports" / "phase20_r3_cp_cell"
FINAL = ROOT / "results" / "phase20_r3" / "opt" / "cp2_madp_final.json"
PL = ReferencePlanes(0.0, None, 1.575e-3)
FREQS = np.arange(8.0, 12.01, 0.5) * 1e9
ANGLES = (0.0, 30.0)
MESHES = (60, 90)  # pitch 0.25 and 0.167 mm for the 15 mm cell
QUANTITIES = ("S_TE_TE", "S_TM_TM", "S_TE_TM", "S_TM_TE", "r_RR", "r_LR", "r_RL", "r_LL")


def _final():
    d = json.loads(FINAL.read_text())
    return build("cp2", d["params_mm"], d["period_mm"]), d


def _db(z):
    return 20 * math.log10(max(abs(z), 1e-12))


def _deg(z):
    return math.degrees(np.angle(z))


def _entries(s11):
    c = circular_reflection(s11)
    return {
        "S_TE_TE": s11[0, 0],
        "S_TM_TM": s11[1, 1],
        "S_TE_TM": s11[0, 1],
        "S_TM_TE": s11[1, 0],
        "r_RR": c[0, 0],
        "r_RL": c[0, 1],
        "r_LR": c[1, 0],
        "r_LL": c[1, 1],
    }


def solve_sets(cell, f, th, n, load_sets):
    base = cell.unit_cell(cell.states()[0], madp_kit())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique_over_loads(
            base, f, MeshSettings(n, n, n // 2 + 1), PL, math.radians(th), 0.0, load_sets
        )


def predictions(cell, kit):
    rows = []
    loads = [cell.unit_cell(st, kit).loads for st in cell.states()]
    for th in ANGLES:
        for f in FREQS:
            per = {n: solve_sets(cell, f, th, n, loads) for n in MESHES}
            for k in range(len(loads)):
                e90, e60 = _entries(per[90][k].s11), _entries(per[60][k].s11)
                for q in QUANTITIES:
                    rows.append(
                        dict(
                            state=k + 1,
                            f_ghz=f / 1e9,
                            theta_deg=th,
                            phi_deg=0.0,
                            quantity=q,
                            mom_mag_dB=_db(e90[q]),
                            mom_phase_deg=_deg(e90[q]),
                            mesh60_mag_dB=_db(e60[q]),
                            mesh60_phase_deg=_deg(e60[q]),
                            mesh_delta_phase_deg=abs(wrap180(_deg(e90[q]) - _deg(e60[q]))),
                            mesh_delta_mag_dB=abs(_db(e90[q]) - _db(e60[q])),
                        )
                    )
    return rows


def tolerance_rows(cell, kit):
    """10 GHz, theta = 0 and 30, mesh 90: the four states at the diode corners (C_T +-20 %, MADP L_s +0.1 nH ASSUMED)."""
    rows = []
    for th in ANGLES:
        kits = corner_kits(kit)
        names = list(kits)
        sets = [cell.unit_cell(st, kits[nm]).loads for nm in names for st in cell.states()]
        res = solve_sets(cell, 10e9, th, 90, sets)
        i = 0
        for nm in names:
            for k in range(4):
                e = _entries(res[i].s11)
                rows.append(
                    dict(
                        corner=nm,
                        state=k + 1,
                        theta_deg=th,
                        r_RR_mag_dB=_db(e["r_RR"]),
                        r_RR_phase_deg=_deg(e["r_RR"]),
                        r_LR_mag_dB=_db(e["r_LR"]),
                    )
                )
                i += 1
    return rows


def lp_interim(cell, kit):
    """LP as 1-bit codes of the CP-2 cell (10 GHz, mesh 90). theta = 0: Cartesian Jones r_xx, r_yx (x-polarised incidence) and r_yy, r_xy
    (y-polarised); theta = 30: TM (x-like) and TE (y-like) co/cross from the 2x2 power-wave matrix."""
    loads = [cell.unit_cell(st, kit).loads for st in cell.states()]
    out = {}
    for th in ANGLES:
        res = solve_sets(cell, 10e9, th, 90, loads)
        rows = []
        for k, r in enumerate(res):
            s = r.s11  # (TE, TM)
            if th == 0.0:
                j = cartesian_reflection(s, 0.0)  # rows reflected (x, y), columns incident (x, y)
                rows.append(
                    {
                        "state": k + 1,
                        "x_inc": {"co": [abs(j[0, 0]), _deg(j[0, 0])], "cross": [abs(j[1, 0]), _deg(j[1, 0])]},
                        "y_inc": {"co": [abs(j[1, 1]), _deg(j[1, 1])], "cross": [abs(j[0, 1]), _deg(j[0, 1])]},
                    }
                )
            else:
                rows.append(
                    {
                        "state": k + 1,
                        "TM_inc": {"co": [abs(s[1, 1]), _deg(s[1, 1])], "cross": [abs(s[0, 1]), _deg(s[0, 1])]},
                        "TE_inc": {"co": [abs(s[0, 0]), _deg(s[0, 0])], "cross": [abs(s[1, 0]), _deg(s[1, 0])]},
                    }
                )
        out[f"theta{th:g}"] = rows
    return out


def bare_slab():
    out = []
    for th in ANGLES:
        for f in FREQS:
            g = LayeredGreenFunction(cell_substrate(), f, backing="pec")
            s0 = oblique_stack_s_parameters(g, PL, math.radians(th), 0.0)
            for i, m in enumerate(("TE", "TM")):
                out.append(
                    dict(
                        f_ghz=f / 1e9,
                        theta_deg=th,
                        floquet_mode=f"{i + 1} ({m})",
                        S11_mag_dB=_db(s0.s11[i, i]),
                        S11_phase_deg=_deg(s0.s11[i, i]),
                    )
                )
    return out


def cell_substrate():
    from src.ris.cells import RT5880

    return RT5880


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(out_dir=OUT, readme_only=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cell, final = _final()
    kit = madp_kit()
    if not readme_only:
        rows = predictions(cell, kit)
        write_csv(out / "mom_reference.csv", rows)
        write_csv(out / "mom_tolerance_10GHz.csv", tolerance_rows(cell, kit))
        write_csv(out / "bare_grounded_slab_reference.csv", bare_slab())
        uc = cell.unit_cell(cell.states()[0], kit)
        res = export_unit_cell_step(uc, out / "cell.step")
        parts = import_named_parts(out / "cell.step")
        (out / "step_parts.json").write_text(
            json.dumps({"parts": parts, "load_bboxes_mm": res.load_bboxes_mm}, indent=1)
        )
        lp = lp_interim(cell, kit)
        (ROOT / "results" / "phase20_r3" / "lp_interim.json").write_text(json.dumps(lp, indent=1))
    else:
        with open(out / "mom_reference.csv", newline="") as fh:
            rows = [{k: (float(v) if k != "quantity" else v) for k, v in r.items()} for r in csv.DictReader(fh)]
    (out / "README.md").write_text(readme(cell, final, rows))


def readme(cell, final, rows):
    sites = cell.diode_sites()
    table = "\n".join(
        f"| Load_{k + 1} | stub at {a:g} deg | ({x * 1e3:+.3f}, {y * 1e3:+.3f}) | {ln * 1e3:.2f} x {w * 1e3:.2f} | {ang:g} |"
        for k, ((x, y, ln, w, ang), a) in enumerate(zip(sites, cell.stub_angles_deg, strict=True))
    )

    def at(state, th, q, f=10.0):
        r = next(
            r
            for r in rows
            if r["state"] == state and r["theta_deg"] == th and r["quantity"] == q and abs(r["f_ghz"] - f) < 1e-9
        )
        return r

    lad = [at(k, 0.0, "r_RR") for k in (1, 2, 3, 4)]
    ladder = " / ".join(f"{r['mom_phase_deg']:.0f}" for r in lad)
    amps = " / ".join(f"{r['mom_mag_dB']:.2f}" for r in lad)
    opp = " / ".join(
        f"{max(at(k, 0.0, 'r_LR')['mom_mag_dB'], at(k, 0.0, 'r_RL')['mom_mag_dB']):.1f}" for k in (1, 2, 3, 4)
    )
    dmax = max(r["mesh_delta_phase_deg"] for r in rows if r["quantity"] == "r_RR" and abs(r["f_ghz"] - 10.0) < 1e-9)
    fx = final["params_mm"]
    return README.format(
        period=cell.period * 1e3,
        circumradius=fx[0],
        line=cell.line * 1e3,
        stubw=fx[2],
        span=cell.diode_span * 1e3,
        pad=cell.pad * 1e3,
        overlap=cell.overlap * 1e3,
        loads=table,
        ladder=ladder,
        amps=amps,
        opp=opp,
        dmax=dmax,
    )


README = """# phase20_r3_cp_cell -- HFSS check of the CP-2 octagon cell with MACOM MADP-000907-14020 diodes (Phase 20 R3)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY), diode values PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** Nothing here has an HFSS number yet; this package
is meant to produce one. It is a **candidate cell, not a design**: the geometry is the optimiser result `results/phase20_r3/OPTIMISATION.md` (CP-2, finalist A), which meets the
phase / amplitude / tolerance targets only within a model uncertainty of order 10 deg and **misses the opposite-handed target** (-13.6 dB at normal incidence against -15 dB).
This same package is also the check for the **LP interim** (LP as 1-bit codes of the CP-2 cell): the Floquet port gives TE and TM, so the linear co/cross numbers come out of the
same eight solves; there is no separate LP geometry (2-bit LP is deferred, see `docs/STATUS.md`).

Diode (MADP-000907-14020, MACOM datasheet Rev V8, DC-0004063, p. 2): ON R_s = 5.2 ohm (typ., 10 mA, 1 GHz), OFF C_T = 0.025 pF (typ., -10 V, 1 MHz; max 0.030 pF = +20 %).
**L_s is NOT_AVAILABLE_FROM_CURRENT_DATA**: the nominal model has L_s = 0 (an assumption); the tolerance file adds L_s = +0.1 nH (ASSUMED_SENSITIVITY_NOT_DEVICE_DATA). The part is a
**flip-chip** device: pads, solder and the real assembly geometry are NOT modelled (the diode is one lumped boundary), and flip-chip assembly is not yet confirmed with fabrication.
Because L_s = 0 each diode is **one single-type boundary** (ON: R only; OFF: C only): the parallel/series question of HFSS Lumped RLC boundaries (R2 package, item 0) does not
enter this package.

## The cell
Period {period:.1f} x {period:.1f} mm; RT5880 eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground under the substrate. Metal: a regular **octagon** (circumradius {circumradius:.3f} mm, flat edges
facing 0/45/90/135 deg ...) plus four radial stubs at **0, 45, 90, 135 deg** (width {stubw:.3f} mm), each = {overlap:.1f} mm overlap into the patch | {pad:.2f} mm pad | **{span:.2f} mm diode span** |
{pad:.2f} mm pad | {line:.3f} mm line. The centre via is omitted (VIA_OMITTED_APPROXIMATION); no bias network. `cell.step` parts (see `step_parts.json`): `Substrate`, `Metal`
(35 um slab: the octagon and the four stubs, five solids), `Ground`, and **`Load_1` ... `Load_4`**: single sheets bridging the diode span of each stub (a rectangle `span x stubw`
rotated with the stub; current line ALONG the stub axis):

| sheet | stub | centre (x, y) [mm] | span x width [mm] | current-line angle [deg from +x] |
|---|---|---|---|---|
{loads}

## HFSS setup
1. Import `cell.step`. **Perfect E** on `Metal` and `Ground` (the MoM metal has zero thickness; state whether you replaced the 35 um slab by sheets at z = 0).
2. Master/slave on the four side walls, scan angle **theta = 0 and 30 deg**, phi = 0 (plane of incidence xz); ONE **Floquet port** on top of an air box (>= lambda/4 at 8 GHz), modes 1 and 2
   (TE: E along y, TM: E along x at phi = 0). **De-embed** to the pattern plane z = 0.
3. Diode boundaries: **Lumped RLC** on each load sheet, current line along the stub (table above), ONE element type each:
   * diode ON : R = 5.2 ohm (L, C unchecked); diode OFF : C = 0.025 pF (R, L unchecked).
   * **State k** (k = 1..4) = `Load_k` ON and the other three OFF (one-hot).  The four states are the same cell turned by 0 / 45 / 90 / 135 deg in the Pancharatnam-Berry sense.
4. Sweep 8-12 GHz (0.5 GHz steps or interpolating), adaptive at 10 GHz, delta S <= 0.01. Export Floquet-mode complex S11 (modes 1, 2: the four entries TE/TM co and cross), 4 states x 2 angles
   = 8 solves. Optional (tolerance, 10 GHz only, 4 states): C_T = 0.030 pF (datasheet max, +20 %), C_T = 0.020 pF (-20 %), and ON with an extra L_s = 0.1 nH (second element, series).

## What to compare (`mom_reference.csv`, `mom_tolerance_10GHz.csv`)
* Columns: `quantity` is `S_TE_TE`, `S_TM_TM`, `S_TE_TM` (TE received, TM incident), `S_TM_TE`, and the circular terms `r_RR`, `r_RL`, `r_LR`, `r_LL`
  (reflected row, incident column; **propagation-relative handedness**, `docs/theory/CONVENTIONS.md` "Circular basis for the RIS cells": incident R = (TM + jTE)/sqrt2, reflected R = (TM - jTE)/sqrt2
  on the power-wave amplitudes). A metal plate has r_LR = r_RL = -1 and r_RR = r_LL = 0; the PB cell wants r_RR = exp(+j 2 alpha) with alpha = state turn.
* `mom_*` = MoM on the 90 x 90 mesh (pitch 0.167 mm); `mesh60_*` the 60 x 60 mesh; `mesh_delta_*` their difference. **The mesh uncertainty is of order 10 deg**: the octagon diagonals are staircased
  and the convergence order was not usable (see `OPTIMISATION.md`); at 10 GHz the 60 -> 90 phase change of r_RR is up to {dmax:.1f} deg here. Do not read a smaller HFSS-MoM difference as agreement to better than that.
* **Convention first:** run the bare grounded slab (no metal, no loads; `bare_grounded_slab_reference.csv`, analytic) at theta = 0, 30 deg and apply any 180 deg TM sign shift HFSS shows to every TM entry
  before comparing; then build r_RR etc. from the corrected TE/TM entries with the formula above. Report what you saw.
* Predictions at 10 GHz, normal incidence, states 1-4: **r_RR phase** {ladder} deg (a rotating ladder, nominally 90 deg steps), **|r_RR|** {amps} dB, **opposite-handed max(|r_LR|, |r_RL|)** {opp} dB.
  The opposite-handed term carries no PB phase (it is the same for every state), so in an array it goes to the specular lobe, not the steered beam (`results/phase20_r3/LOCALITY.md`).
* What would be news: an r_RR ladder step differing from ~90 deg by more than ~20 deg (the model uncertainty is ~10 deg), |r_RR| below -2 dB (diode loss: R_s = 5.2 ohm is the dominant loss), or an
  opposite-handed level much worse than -13 dB. The expected bandwidth is limited by the diode loss: |r_RR| >= -1.5 dB over about 10.0-10.75 GHz (MoM grid), the ladder steps staying within 90 +- 20 deg across 8.5-11.5 GHz.
* LP interim (`results/phase20_r3/lp_interim.json`): with x-polarised incidence the co-LP and cross-LP terms of the four states follow from the same S11 (TM = x at phi = 0).
* Honest limits: zero-thickness metal; staircased octagon in the MoM vs exact in HFSS; omitted via and bias network; lumped diode without flip-chip pad geometry; L_s assumed 0.
"""

if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(*(a[:1] or [OUT]), readme_only="--readme-only" in sys.argv)
