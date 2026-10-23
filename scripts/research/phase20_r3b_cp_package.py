"""Phase 20 R3b: the HFSS check package of the re-optimised CP-2 cell (`hfss_exports/phase20_r3b_cp_cell/`). Predictions first; nothing here has an HFSS number yet.
MODEL_ONLY (raw MoM at F = 2n, mesh 60 / 90) with a clearly separated HFSS_CALIBRATED_TRANSFERRED prediction (baseline-cell HFSS delta transferred), PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.

    python scripts/research/phase20_r3b_cp_package.py [out_dir] [--final results/phase20_r3b/cp2_r3b_final.json] [--smoke]

`--smoke` solves one frequency only (test). Frequencies: 8-12 GHz in 0.25 GHz steps; theta = 0 and 30 deg, phi = 0; 4 one-hot states; mesh 60 and 90 at F = 2n.
"""

import json
import math
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase20_r3_cp_package as r3  # noqa: E402
import phase20_r3b_opt as opt  # noqa: E402

from src.basis_functions import MeshSettings  # noqa: E402
from src.cad_export import export_unit_cell_step, import_named_parts  # noqa: E402
from src.ris.cells import DiodeKit, madp_kit  # noqa: E402
from src.ris.hfss_calibration import DeltaTable  # noqa: E402
from src.ris.polarization import circular_reflection, to_structure_frame  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads  # noqa: E402

OUT = ROOT / "hfss_exports" / "phase20_r3b_cp_cell"
FINAL = ROOT / "results" / "phase20_r3b" / "cp2_r3b_final.json"
FREQS = np.arange(8.0, 12.01, 0.25)
ANGLES = (0.0, 30.0)
MESHES = (60, 90)


def corner_set():
    nom = madp_kit()
    return {
        "nominal": nom,
        "C_T+20% (0.030 pF, datasheet max)": DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=1.2)),
        "C_T-20% (0.020 pF, ASSUMED)": DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=0.8)),
        "C_T-10% (0.0225 pF, ASSUMED)": DiodeKit(nom.name, nom.on, nom.off.scaled(C_T=0.9)),
        "L_s+0.1nH (ASSUMED)": madp_kit(0.1e-9),
    }


def solve(task):
    x, f_ghz, th, n, names = task
    cell = opt.make_cell(x)
    ks = corner_set()
    sets = [cell.unit_cell(st, ks[k]).loads for k in names for st in cell.states()]
    base = cell.unit_cell(cell.states()[0], ks["nominal"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_s_parameters_oblique_over_loads(
            base, f_ghz * 1e9, MeshSettings(n, n, 2 * n), opt.PL, math.radians(th), 0.0, sets, warn_grating=False
        )
    return np.array([r.s11 for r in res])


def rows_for(x, freqs, workers):
    tasks = [(tuple(x), float(f), th, n, ["nominal"]) for th in ANGLES for f in freqs for n in MESHES]
    with Pool(workers) as pool:
        res = pool.map(solve, tasks, chunksize=1)
    by = {(t[1], t[2], t[3]): r for t, r in zip(tasks, res, strict=True)}
    delta = DeltaTable.load(ROOT / "results/phase20_r3b/hfss_delta_table.npz")
    rows, cal_rows = [], []
    for th in ANGLES:
        for f in freqs:
            s90, s60 = by[(float(f), th, 90)], by[(float(f), th, 60)]
            for k in range(4):
                e90, e60 = r3._entries(s90[k]), r3._entries(s60[k])
                # transferred HFSS calibration (baseline-cell delta added to the mesh-90 circular response, phi = 0) -- NOT an HFSS number
                c = to_structure_frame(circular_reflection(s90[k]), 0.0)
                try:
                    c_cal, label = delta.apply(c, float(f), th, k + 1)
                except Exception:
                    c_cal, label = c, "NOT_AVAILABLE_FROM_CURRENT_DATA"
                for q in r3.QUANTITIES:
                    rows.append(
                        dict(
                            state=k + 1,
                            f_ghz=float(f),
                            theta_deg=th,
                            phi_deg=0.0,
                            quantity=q,
                            mom_mag_dB=r3._db(e90[q]),
                            mom_phase_deg=r3._deg(e90[q]),
                            mesh60_mag_dB=r3._db(e60[q]),
                            mesh60_phase_deg=r3._deg(e60[q]),
                            mesh_delta_phase_deg=abs(r3.wrap180(r3._deg(e90[q]) - r3._deg(e60[q]))),
                            mesh_delta_mag_dB=abs(r3._db(e90[q]) - r3._db(e60[q])),
                        )
                    )
                for q, (i, j) in {"r_RR": (0, 0), "r_RL": (0, 1), "r_LR": (1, 0), "r_LL": (1, 1)}.items():
                    cal_rows.append(
                        dict(
                            state=k + 1,
                            f_ghz=float(f),
                            theta_deg=th,
                            phi_deg=0.0,
                            quantity=q,
                            transferred_mag_dB=r3._db(c_cal[i, j]),
                            transferred_phase_deg=r3._deg(c_cal[i, j]),
                            label="HFSS_CALIBRATED_TRANSFERRED" if label == "HFSS_CALIBRATED" else label,
                        )
                    )
    return rows, cal_rows


def tolerance_rows(x, workers):
    ks = corner_set()
    names = list(ks)
    tasks = [(tuple(x), 10.5, th, 90, names) for th in ANGLES]
    with Pool(min(workers, 2)) as pool:
        res = pool.map(solve, tasks, chunksize=1)
    rows = []
    for th, s in zip(ANGLES, res, strict=True):
        for ci, nm in enumerate(names):
            for k in range(4):
                e = r3._entries(s[4 * ci + k])
                rows.append(
                    dict(
                        corner=nm,
                        state=k + 1,
                        theta_deg=th,
                        f_ghz=10.5,
                        r_RR_mag_dB=r3._db(e["r_RR"]),
                        r_RR_phase_deg=r3._deg(e["r_RR"]),
                        r_LR_mag_dB=r3._db(e["r_LR"]),
                        r_RL_mag_dB=r3._db(e["r_RL"]),
                    )
                )
    return rows


def main(out_dir=OUT, final=FINAL, smoke=False, workers=11):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    d = json.loads(Path(final).read_text())
    x = d["x_mm"]
    cell = opt.make_cell(x)
    freqs = [10.5] if smoke else list(FREQS)
    rows, cal = rows_for(x, freqs, workers)
    r3.write_csv(out / "mom_reference.csv", rows)
    r3.write_csv(out / "mom_reference_hfss_calibrated_transferred.csv", cal)
    r3.write_csv(out / "mom_tolerance_10p5GHz.csv", tolerance_rows(x, workers))
    r3.write_csv(out / "bare_grounded_slab_reference.csv", bare_slab(freqs))
    uc = cell.unit_cell(cell.states()[0], madp_kit())
    res = export_unit_cell_step(uc, out / "cell.step")
    parts = import_named_parts(out / "cell.step")
    (out / "step_parts.json").write_text(json.dumps({"parts": parts, "load_bboxes_mm": res.load_bboxes_mm}, indent=1))
    (out / "README.md").write_text(readme(cell, d, rows))


def bare_slab(freqs):
    from src.green_function import LayeredGreenFunction
    from src.sparameters import oblique_stack_s_parameters

    out = []
    for th in ANGLES:
        for f in freqs:
            g = LayeredGreenFunction(r3.cell_substrate(), f * 1e9, backing="pec")
            s0 = oblique_stack_s_parameters(g, opt.PL, math.radians(th), 0.0)
            for i, m in enumerate(("TE", "TM")):
                out.append(
                    dict(
                        f_ghz=f,
                        theta_deg=th,
                        floquet_mode=f"{i + 1} ({m})",
                        S11_mag_dB=r3._db(s0.s11[i, i]),
                        S11_phase_deg=r3._deg(s0.s11[i, i]),
                    )
                )
    return out


def readme(cell, d, rows):
    sites = cell.diode_sites()
    table = "\n".join(
        f"| Load_{k + 1} | stub at {a:g} deg | ({x * 1e3:+.3f}, {y * 1e3:+.3f}) | {ln * 1e3:.2f} x {w * 1e3:.2f} | {ang:g} |"
        for k, ((x, y, ln, w, ang), a) in enumerate(zip(sites, cell.stub_angles_deg, strict=True))
    )

    def at(state, th, q, f=10.5):
        return next(
            r
            for r in rows
            if r["state"] == state and r["theta_deg"] == th and r["quantity"] == q and abs(r["f_ghz"] - f) < 1e-9
        )

    ladder = " / ".join(f"{at(k, 0.0, 'r_RR')['mom_phase_deg']:.0f}" for k in (1, 2, 3, 4))
    amps = " / ".join(f"{at(k, 0.0, 'r_RR')['mom_mag_dB']:.2f}" for k in (1, 2, 3, 4))
    opp = " / ".join(
        f"{max(at(k, 0.0, 'r_LR')['mom_mag_dB'], at(k, 0.0, 'r_RL')['mom_mag_dB']):.1f}" for k in (1, 2, 3, 4)
    )
    return README.format(
        period=cell.period * 1e3,
        circumradius=cell.circumradius * 1e3,
        line_axis=cell.line * 1e3,
        w_axis=cell.stub_width * 1e3,
        line_diag=cell.line_at(45.0) * 1e3,
        w_diag=cell.width_at(45.0) * 1e3,
        span=cell.diode_span * 1e3,
        pad=cell.pad * 1e3,
        overlap=cell.overlap * 1e3,
        loads=table,
        ladder=ladder,
        amps=amps,
        opp=opp,
        source=d.get("source", ""),
    )


README = """# phase20_r3b_cp_cell -- HFSS check of the re-optimised CP-2 octagon cell, 10.5 GHz (Phase 20 R3b)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY), diode values PROVISIONAL_DATASHEET, VIA_OMITTED_APPROXIMATION.** Generated at Floquet cutoff F = 2n (converged; see `results/phase20_r4/FLOQUET_CONVERGENCE.md`).
This is the second HFSS check of the RIS cell. The first (`hfss_results/phase20_r3_cp_cell/`) showed the PB ladder working, best point 10.5 GHz, but unequal steps (98/92/102/68 deg) and theta-30 opposite-handed
-11.9 dB; R3b changed the cell (`docs/briefs/PHASE20_R3B_CP2_REOPT.md`): same 15 mm period and diode, **separate axis and diagonal stub length and width**. {source}

Diode (MACOM MADP-000907-14020): ON R_s = 5.2 ohm, OFF C_T = 0.025 pF (max 0.030 pF = +20 %; the datasheet gives **no minimum**, so the C_T -10 % / -20 % cases are our ASSUMPTION); L_s NOT_AVAILABLE_FROM_CURRENT_DATA
(nominal 0, +0.1 nH ASSUMED sensitivity). Flip-chip assembly not confirmed with fabrication; one lumped boundary per diode. Each diode is one single-type boundary (ON: R only; OFF: C only).

## The cell
Period {period:.1f} x {period:.1f} mm; RT5880 eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground. Metal: regular **octagon** (circumradius {circumradius:.3f} mm) plus four radial stubs at 0, 45, 90, 135 deg, each
{overlap:.1f} mm overlap | {pad:.2f} mm pad | **{span:.2f} mm diode span** | {pad:.2f} mm pad | line. **Axis stubs (0, 90 deg): line {line_axis:.3f} mm, width {w_axis:.3f} mm. Diagonal stubs (45, 135 deg): line {line_diag:.3f} mm,
width {w_diag:.3f} mm.** Centre via omitted (VIA_OMITTED_APPROXIMATION); no bias network. `cell.step` parts (`step_parts.json`): `Substrate`, `Metal`, `Ground`, `Load_1` ... `Load_4` (sheets bridging the diode span; current
line ALONG the stub axis):

| sheet | stub | centre (x, y) [mm] | span x width [mm] | current-line angle [deg from +x] |
|---|---|---|---|---|
{loads}

## HFSS setup
1. Import `cell.step`; **Perfect E** on `Metal` and `Ground`.
2. Master/slave on the side walls, scan angle **theta = 0 and 30 deg**, phi = 0; ONE Floquet port, modes 1 (TE) and 2 (TM). **Fix the de-embedding distance** (the first run stopped ~1.24 mm short of the patch plane); state the
   distance you used so that no offset has to be fitted (we still fit it from the bare slab).
3. Lumped RLC on each load sheet along the stub; ON: R = 5.2 ohm; OFF: C = 0.025 pF (one element type each). State k = `Load_k` ON, the other three OFF.
4. Sweep **8-12 GHz in 0.25 GHz steps** (17 points; the first run's 0.5 GHz sampling undersampled the sharp features), adaptive at 10.5 GHz. Export the Floquet-mode complex S (modes 1, 2: all four entries), 4 states x 2 angles,
   and the bare slab (no metal, no loads) for the convention check.
5. **Tolerance cases at 10.5 GHz, theta = 0 and 30, 4 states each** (`mom_tolerance_10p5GHz.csv`): C_T = 0.030 pF (+20 %, datasheet max), C_T = 0.020 pF (-20 %, assumed), C_T = 0.0225 pF (-10 %, assumed), and ON with
   an extra L_s = 0.1 nH (series element, assumed).

## What to compare
* `mom_reference.csv`: raw MoM at F = 2n, mesh 90 (`mom_*`) and mesh 60 (`mesh60_*`), long format as the first package (`quantity` in S_TE_TE, S_TM_TM, S_TE_TM (TE received, TM incident), S_TM_TE, r_RR, r_RL, r_LR, r_LL; propagation-relative handedness).
* `mom_reference_hfss_calibrated_transferred.csv`: the same circular entries with the **baseline cell's HFSS delta added** (label `HFSS_CALIBRATED_TRANSFERRED`: an assumption that the HFSS-vs-MoM offsets of the first cell carry over; NOT an HFSS
  number). The second HFSS run tells whether the transfer holds (`src/ris/hfss_calibration.py`).
* Predictions at 10.5 GHz, normal incidence, states 1-4: **r_RR phase** {ladder} deg, **|r_RR|** {amps} dB, **opposite-handed max(|r_LR|, |r_RL|)** {opp} dB (raw MoM, mesh 90).
* Bare slab convention first (`bare_grounded_slab_reference.csv`), as before. What would be news: ladder steps differing from ~90 deg by more than ~20 deg, |r_RR| below -2 dB, opposite-handed worse than -11 dB at theta = 30.
* Honest limits: zero-thickness metal; staircased octagon / 45 deg stubs in the MoM; omitted via and bias network; lumped diode without flip-chip pad geometry; L_s assumed 0; model uncertainty about 10 deg.
"""

if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    fin = sys.argv[sys.argv.index("--final") + 1] if "--final" in sys.argv else FINAL
    a = [x for x in a if x != str(fin)]
    wk = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 11
    a = [x for x in a if x != str(wk)]
    main(*(a[:1] or [OUT]), final=fin, smoke="--smoke" in sys.argv, workers=wk)
