"""Phase 20 R2 item 5 (user-approved after item 4): the one-cell HFSS check package for the SMP1345-079LF ON and OFF cell
(`hfss_exports/phase20_r2_diode_cell/`). Predictions first, nothing here has an HFSS number yet.

    python scripts/research/phase20_r2_hfss_package.py [out_dir]

Cell: R1 demo cell (10 mm, RT5880 eps_r 2.2 tan d 0.0009 h 1.575 mm, PEC ground), 8 x 8 mm patch with a 2 mm slot across x; the slot is bridged
by a 0.5 mm wide strip = 0.5 mm pad | diode chain (span 1.0 mm) | 0.5 mm pad. The chain is `src.ris.diodes.hfss_chain`: Load_1 (0.25 mm) |
island (0.5 mm) | Load_2 (0.25 mm); ON: Load_1 = R 1.5 ohm, Load_2 = L 0.7 nH; OFF: Load_1 = C 0.19 pF, Load_2 = L 0.7 nH. Each HFSS boundary has ONE
element type, so the parallel/series question of HFSS RLC boundaries does not enter. Labels: MODEL_ONLY, NOT_HFSS_VALIDATED, PROVISIONAL_DATASHEET.
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

from src.basis_functions import MeshSettings  # noqa: E402
from src.cad_export import export_unit_cell_step, import_named_parts  # noqa: E402
from src.composite_geometry import AxisRect, CompositeGeometry  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.geometry import UnitCell  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.periodic_bc import RectangularLattice  # noqa: E402
from src.ris.diodes import diode_load, hfss_chain, smp1345_off, smp1345_on  # noqa: E402
from src.sparameters import compute_s_parameters_oblique_over_loads, oblique_stack_s_parameters  # noqa: E402

OUT = ROOT / "hfss_exports/phase20_r2_diode_cell"
SUB = Substrate("rt5880", 2.2, 0.0009, 1.575e-3)
SPAN, PAD, STRIP_W, ELEMENT, ISLAND = 1.0e-3, 0.5e-3, 0.5e-3, 0.25e-3, 0.5e-3
FREQS = np.arange(8.0, 12.01, 0.5) * 1e9
ANGLES = (0.0, 30.0)
MODES = ("TE", "TM")
MESHES = (
    (40, 24),
    (80, 41),
)  # edges at 0.25 mm multiples: only n a multiple of 40 is aligned with the strip (see README)


def cell(loads=()):
    half = SPAN / 2 + PAD
    arms = (
        AxisRect(-(half + 1.5e-3), 0.0, 3.0e-3, 8e-3),
        AxisRect(half + 1.5e-3, 0.0, 3.0e-3, 8e-3),
        AxisRect(0.0, 0.0, 2 * half, STRIP_W),
    )
    return UnitCell(
        RectangularLattice(10e-3, 10e-3), CompositeGeometry(arms), SUB, Metal(), backing="pec", loads=tuple(loads)
    )


def states():
    on, off = smp1345_on(), smp1345_off()
    return {"ON": on, "OFF": off}


def chain(model):
    return tuple(hfss_chain(model, 0.0, 0.0, ELEMENT, ISLAND, STRIP_W)[0])


def solve_all(f, th, mesh_n, floquet):
    """Co-pol S11 (2x2 TE/TM) for [ON chain, OFF chain, ON ideal, OFF ideal] at one frequency and angle."""
    ms = states()
    sets = [
        chain(ms["ON"]),
        chain(ms["OFF"]),
        (diode_load(ms["ON"], 0, 0, SPAN, STRIP_W),),
        (diode_load(ms["OFF"], 0, 0, SPAN, STRIP_W),),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return compute_s_parameters_oblique_over_loads(
            cell(),
            f,
            MeshSettings(mesh_n, mesh_n, floquet),
            ReferencePlanes(0.0, None, SUB.thickness),
            math.radians(th),
            0.0,
            sets,
        )


def _db(z):
    return 20 * math.log10(abs(z))


def _deg(z):
    return math.degrees(np.angle(z))


def _extrap(a, b, p):
    """S_inf from meshes 40 and 80 assuming S(n) = S_inf + c / n^p."""
    return (80.0**p * b - 40.0**p * a) / (80.0**p - 40.0**p)


P_ASSUMED = (
    1.5  # median observed order at 10 GHz (meshes 40/80/120): 1.2 (OFF TM) ... 1.65, see results/phase20_r2/DEMO.md
)


def predictions():
    rows = []
    for th in ANGLES:
        for f in FREQS:
            per = {n: solve_all(f, th, n, fl) for n, fl in MESHES}
            for si, (state, kind) in enumerate(
                (("ON", "hfss_chain"), ("OFF", "hfss_chain"), ("ON", "ideal_series"), ("OFF", "ideal_series"))
            ):
                for i, m in enumerate(MODES):
                    a, b = per[40][si].s11[i, i], per[80][si].s11[i, i]
                    e = _extrap(a, b, P_ASSUMED)
                    lo, hi = _extrap(a, b, 1.0), _extrap(a, b, 2.0)
                    rows.append(
                        dict(
                            state=state,
                            model=kind,
                            f_ghz=f / 1e9,
                            theta_deg=th,
                            phi_deg=0.0,
                            floquet_mode=f"{i + 1} ({m})",
                            mom_S11_mag_dB=_db(e),
                            mom_S11_phase_deg=_deg(e),
                            uncertainty_phase_deg=abs(_deg(lo / hi)),
                            uncertainty_mag_dB=abs(_db(lo) - _db(hi)),
                            mesh80_mag_dB=_db(b),
                            mesh80_phase_deg=_deg(b),
                            mesh40_mag_dB=_db(a),
                            mesh40_phase_deg=_deg(a),
                        )
                    )
    return rows


def bare_slab():
    out = []
    for th in ANGLES:
        for f in FREQS:
            g = LayeredGreenFunction(SUB, f, backing="pec")
            s0 = oblique_stack_s_parameters(g, ReferencePlanes(0.0, None, SUB.thickness), math.radians(th), 0.0)
            for i, m in enumerate(MODES):
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


def rlc_check_numbers(f=10e9, r=50.0, c=0.1e-12):
    w = 2 * math.pi * f
    zp = 1.0 / (1.0 / r + 1j * w * c)
    zs = r + 1.0 / (1j * w * c)
    s = lambda z: (z - 50.0) / (z + 50.0)  # noqa: E731
    return dict(zp=zp, zs=zs, sp=s(zp), ss=s(zs))


def contrast_sentence(chain_rows):
    def wrap(x):
        return (x + 180.0) % 360.0 - 180.0

    tm = {
        (r["state"], r["f_ghz"], r["theta_deg"]): r["mom_S11_phase_deg"]
        for r in chain_rows
        if r["floquet_mode"].endswith("(TM)")
    }
    parts = []
    for th in ANGLES:
        diff = {f: wrap(tm[("ON", f, th)] - tm[("OFF", f, th)]) for f in sorted({k[1] for k in tm})}
        fmax = max(diff, key=lambda f: abs(diff[f]))
        parts.append(
            f"theta {th:g} deg: {diff[10.0]:+.0f} deg at 10 GHz, largest {diff[fmax]:+.0f} deg at {fmax:g} GHz"
        )
    return "TM ON-OFF phase difference of the extrapolated MoM prediction (hfss_chain): " + "; ".join(parts) + "."


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(out_dir=OUT, readme_only=False):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if readme_only:
        with open(out / "mom_reference.csv", newline="") as fh:
            rows = [
                {k: (float(v) if k not in ("state", "model", "floquet_mode") else v) for k, v in r.items()}
                for r in csv.DictReader(fh)
            ]
    else:
        rows = predictions()
        write_csv(out / "mom_reference.csv", rows)
    if not readme_only:
        write_csv(out / "bare_grounded_slab_reference.csv", bare_slab())
        c = cell(chain(states()["ON"]))
        res = export_unit_cell_step(c, out / "cell.step")
        parts = import_named_parts(out / "cell.step")
        (out / "step_parts.json").write_text(
            json.dumps({"parts": parts, "load_bboxes_mm": res.load_bboxes_mm}, indent=1)
        )
    ch = rlc_check_numbers()
    chain_rows = [r for r in rows if r["model"] == "hfss_chain"]
    ideal_rows = {
        (r["state"], r["f_ghz"], r["theta_deg"], r["floquet_mode"]): r for r in rows if r["model"] == "ideal_series"
    }
    d = max(
        abs(
            (
                (
                    r["mom_S11_phase_deg"]
                    - ideal_rows[(r["state"], r["f_ghz"], r["theta_deg"], r["floquet_mode"])]["mom_S11_phase_deg"]
                    + 180
                )
                % 360
            )
            - 180
        )
        for r in chain_rows
        if r["state"] == "OFF" and r["floquet_mode"].endswith("(TM)")
    )
    (out / "README.md").write_text(
        README.format(
            zp=abs(ch["zp"]),
            zpa=math.degrees(np.angle(ch["zp"])),
            zs=abs(ch["zs"]),
            zsa=math.degrees(np.angle(ch["zs"])),
            sp=abs(ch["sp"]),
            ss=abs(ch["ss"]),
            unc_te=max(r["uncertainty_phase_deg"] for r in chain_rows if r["floquet_mode"].endswith("(TE)")),
            unc_tm=max(r["uncertainty_phase_deg"] for r in chain_rows if r["floquet_mode"].endswith("(TM)")),
            chain_vs_ideal=d,
            contrast=contrast_sentence(chain_rows),
        )
    )
    return rows


README = """# phase20_r2_diode_cell -- one-cell HFSS check of the SMP1345-079LF ON and OFF cell (Phase 20 R2)

**State: PREDICTION_NOT_HFSS_VALIDATED (MODEL_ONLY), diode values PROVISIONAL_DATASHEET.** Nothing here has an HFSS number yet; this package is meant
to produce one. It is a **demonstration cell, not a design**. Diode values: Skyworks SMP1345 Series datasheet 200046U (2018-06-19): R_s 1.5 ohm
(10 mA, 100 MHz, Table 3), C_T 0.19 pF (1 V, 1 MHz, Table 3; C_T at 0 V is not given), L_s 0.7 nH (SC-79, Table 1); R_p not available (OFF is lossless).
The datasheet covers 10 MHz-6 GHz: these values are **extrapolated** to 8-12 GHz. Details and sources: `results/phase20_r2/README.md`.

## 0. First: the one-line RLC-boundary check (do this before anything else)
HFSS Lumped RLC boundaries are documented here as R, L, C in **parallel**; I could not confirm that for your HFSS version (the one source found does not say).
The package avoids depending on it (every boundary carries ONE element type), but verify the semantic anyway, because it decides whether an R_p could ever be added:
terminate a 50 ohm lumped port with a 1 mm x 0.5 mm sheet carrying one Lumped RLC boundary with **R = 50 ohm and C = 0.1 pF both checked**, and read S11 at 10 GHz:
* parallel R||C: |Z| = {zp:.1f} ohm at {zpa:.1f} deg  ->  |S11| = {sp:.3f}
* series R+C:   |Z| = {zs:.1f} ohm at {zsa:.1f} deg  ->  |S11| = {ss:.3f}
Report which one you get. Expected for the package itself: nothing changes either way (single-type boundaries).

## The cell
Period 10 x 10 mm; RT5880 eps_r 2.2, tan d 0.0009, h 1.575 mm; PEC ground under the substrate; a metal patch 8 x 8 mm (x, y) with a **2.0 mm slot across x** at the centre,
bridged by a **0.5 mm wide strip** (the two 0.5 mm **pads** at the slot edges, the 1.0 mm diode span in the middle). Pads are required: with the boundaries touching the wide arms
the 2-D current crowding changes the OFF phase by tens of degrees against the ideal series model (`results/phase20_r2/README.md`). `cell.step` parts (see `step_parts.json`):
`Substrate` (z -1.575..0), `Metal` (35 um slab z 0..0.035: the two arms plus the pads and the 0.5 mm **island** between the boundaries, three solids), `Ground`,
`Load_1` and `Load_2` (single sheets at z = 0.0175 mm, 0.25 mm long in x, 0.5 mm wide in y, centred at x = -0.375 and +0.375 mm, bridging the metal faces).

## HFSS setup
1. Import `cell.step`. **Perfect E** on `Metal` and `Ground` (the MoM metal has zero thickness; state whether you replaced the 35 um slab by sheets at z = 0).
2. Master/slave on the four side walls, **scan angle theta = 0 and 30 deg**, phi = 0 (plane of incidence xz); air box above (>= lambda/4 at 8 GHz) with ONE **Floquet port**
   on its top face, modes 1 and 2 (TE and TM w.r.t. the plane of incidence: at phi = 0 TE has E along y, TM has E along x, i.e. across the slot). **De-embed** to the pattern plane z = 0.
3. Diode boundaries: **Lumped RLC** on each load sheet, current line along +x across the sheet. One element type each:

| state | Load_1 (x = -0.375 mm) | Load_2 (x = +0.375 mm) |
|---|---|---|
| ON  | R = 1.5 ohm (L, C unchecked) | L = 0.7 nH (R, C unchecked) |
| OFF | C = 0.19 pF (R, L unchecked) | L = 0.7 nH (R, C unchecked) |

4. Frequency sweep 8-12 GHz (0.5 GHz steps or interpolating); adaptive solve at 10 GHz, delta S <= 0.01. Export Floquet-mode S11 (magnitude dB and phase), modes 1 and 2, theta = 0, 30 deg, ON and OFF.
   That is 8 solves-worth of data: 2 states x 2 angles, TE and TM from the same port.

## What to compare (`mom_reference.csv`)
Rows `model = hfss_chain` are the MoM of exactly this geometry (two loads and the island); `ideal_series` is the single-load series model of the same diode (the model used in `results/phase20_r2/DEMO.md`).
`mom_S11_*` are **extrapolated** from the 40 x 40 and 80 x 80 meshes (both aligned with the 0.25 mm strip edges; Floquet order 24 / 41) assuming S(n) = S_inf + c/n^p with p = 1.5, the median order
observed at 10 GHz on meshes 40/80/120 (1.2-1.65 per row, `results/phase20_r2/DEMO.md`); `uncertainty_*` is the spread between the p = 1 and p = 2 extrapolations (up to {unc_te:.1f} deg TE,
{unc_tm:.1f} deg TM); the raw 40 and 80 mesh values are listed too. TE is nearly state-independent for this narrow bridge (an HFSS TE ON/OFF difference of more than a few degrees would be news),
but the TE phase itself is the least converged quantity (several degrees between meshes).
1. **Convention first:** run the bare grounded slab (no metal, no loads; `bare_grounded_slab_reference.csv`, analytic) at theta = 0, 30 deg. If HFSS's TM phase differs from the listed value by 180 deg (sign
   convention of the TM mode), apply the same shift to every TM entry before comparing; TE should agree directly. Report what you saw.
2. Then ON and OFF. {contrast} (The OFF series resonance of L_s and C_T is 13.8 GHz nominal.) This cell's contrast is NOT that of the item-4 demo cell (`results/phase20_r2/DEMO.md`, direct 0.5 mm bridge, no pads): the 2 mm slot and the pads change it.
3. **OFF is a sensitive state**: +-20 % on C_T or L_s moves its phase by tens of degrees (up to >100 deg). A difference of that size between HFSS and the MoM is most likely the diode values (C_T at the OFF bias, the
   lead inductance of the pads and island), not a solver error: do the HFSS OFF run with C_T 0.19 pF and also 0.23 pF (+20 %) so the sensitivity can be read directly.
4. Honest limits: zero-thickness metal in the MoM; the island and pads are conductor strips, the loads are series sheets (`transverse_conducts=False`); the chain and the ideal series model differ by up to {chain_vs_ideal:.1f} deg
   (OFF, TM) in the MoM itself; the MoM mesh is aligned only at n = 40, 80, 120 (mesh 160 does not fit in 36 GiB)
"""

if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    main(*(a[:1] or [OUT]), readme_only="--readme-only" in sys.argv)
