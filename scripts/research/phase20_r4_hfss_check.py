"""First HFSS check of the CP-2 RIS cell: the user's HFSS run (`hfss_results/phase20_r3_cp_cell/`) against the F = 2n state library and the regenerated F = 2n `mom_reference.csv`.
COMPARISON ONLY (no calibration is applied here). Writes results/phase20_r4/hfss_check.json and HFSS_CHECK.md.

    python scripts/research/phase20_r4_hfss_check.py [hfss_dir]

Conventions fitted from the HFSS bare grounded slab (state 1 rows of the bare-slab file): Floquet mode 1 = TE, mode 2 = TM; reference-plane offset d with HFSS phase = analytic - 2 k0 d cos(theta)
(removed by multiplying every S by exp(+j 2 k0 d cos theta)); cross-term sign tested (+1 / -1) on r_RR against the library.
"""

import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.ris.hfss_calibration import compare, fit_reference_plane_mm, proposed_offsets, read_hfss_reim  # noqa: E402
from src.ris.polarization import circular_reflection  # noqa: E402
from src.ris.state_library import StateLibrary  # noqa: E402

PKG = ROOT / "hfss_exports" / "phase20_r3_cp_cell"
LIB = ROOT / "results" / "phase20_r4" / "cp2_library_v1.npz"
FREQS = [8.0 + 0.5 * i for i in range(9)]
STATES = (1, 2, 3, 4)


def db(x):
    return 20 * np.log10(max(abs(x), 1e-9))


def ph(x):
    return float(np.degrees(np.angle(x)))


def steps(phases):
    return [float((phases[(i + 1) % 4] - phases[i]) % 360.0) for i in range(4)]


def metrics(cs):
    """cs: four 2x2 circular matrices [reflected, incident] of the four states."""
    st = steps([ph(c[0, 0]) for c in cs])
    return {
        "min_rr_db": min(db(c[0, 0]) for c in cs),
        "opp_db": max(max(db(c[0, 1]), db(c[1, 0])) for c in cs),
        "steps_deg": st,
        "step_error_deg": max(abs(s - 90.0) for s in st),
    }


def mom_reference(path):
    """{(state, f, theta, phi): 2x2 TE/TM} from the package's long-format `mom_reference.csv` (mesh-90 values at F = 2n, no Richardson correction)."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["quantity"] not in ("S_TE_TE", "S_TM_TM", "S_TE_TM", "S_TM_TE"):
                continue
            key = (int(float(r["state"])), round(float(r["f_ghz"]), 6), float(r["theta_deg"]), float(r["phi_deg"]))
            m = out.setdefault(key, np.zeros((2, 2), dtype=complex))
            i, j = {"S_TE_TE": (0, 0), "S_TM_TM": (1, 1), "S_TE_TM": (0, 1), "S_TM_TE": (1, 0)}[r["quantity"]]
            m[i, j] = 10 ** (float(r["mom_mag_dB"]) / 20) * np.exp(1j * np.radians(float(r["mom_phase_deg"])))
    return out


def shift_fit(h, lib, what, theta=None, window=(9.5, 11.0)):
    """Frequency shift df (GHz) minimising the rms dB difference between HFSS(f) and library(f - df) for |r_RR| ('amp') or the opposite-handed level ('opp')."""

    def rms(df):
        e = []
        for (st, f, th, phi), m in h.items():
            if not window[0] <= f <= window[1] or (theta is not None and th != theta):
                continue
            cl = lib.structure_circular(f - df, th, phi, st)[0]
            c = circular_reflection(m)
            if what == "amp":
                e.append(db(c[0, 0]) - db(cl[0, 0]))
            else:
                e.append(max(db(c[0, 1]), db(c[1, 0])) - max(db(cl[0, 1]), db(cl[1, 0])))
        return float(np.sqrt(np.mean(np.square(e))))

    xs = np.arange(-0.5, 0.501, 0.05)
    es = [rms(x) for x in xs]
    i = int(np.argmin(es))
    return {"best_df_ghz": float(xs[i]), "rms_db_best": es[i], "rms_db_at_zero": rms(0.0)}


def main(hdir):
    hdir = Path(hdir)
    d_mm, resid = fit_reference_plane_mm(
        hdir / "hfss_bare_slab_theta0_30.csv", PKG / "bare_grounded_slab_reference.csv"
    )
    h = read_hfss_reim(hdir / "hfss_cp2_states1-4_theta0_30.csv", d_mm)
    lib = StateLibrary.load(LIB)
    mom = mom_reference(PKG / "mom_reference.csv")

    sign = {}
    for sx in (1, -1):
        e = []
        for (st, f, th, phi), m in h.items():
            m2 = m.copy()
            m2[0, 1] *= sx
            m2[1, 0] *= sx
            e.append(abs(circular_reflection(m2)[0, 0] - lib.structure_circular(f, th, phi, st)[0][0, 0]))
        sign[str(sx)] = float(np.mean(e))

    comps, summary = compare(lib, h)
    table = []
    for th in (0.0, 30.0):
        for f in FREQS:
            hc = [circular_reflection(h[(s, round(f, 6), th, 0.0)]) for s in STATES]
            lc = [lib.structure_circular(f, th, 0.0, s)[0] for s in STATES]
            mc = [circular_reflection(mom[(s, round(f, 6), th, 0.0)]) for s in STATES]
            table.append(
                {"f_ghz": f, "theta_deg": th, "hfss": metrics(hc), "library": metrics(lc), "mom_mesh90": metrics(mc)}
            )
    at105 = []
    for th in (0.0, 30.0):
        for s in STATES:
            c = circular_reflection(h[(s, 10.5, th, 0.0)])
            lc = lib.structure_circular(10.5, th, 0.0, s)[0]
            at105.append(
                {
                    "state": s,
                    "theta_deg": th,
                    "rr_level_hfss_db": db(c[0, 0]),
                    "rr_level_library_db": db(lc[0, 0]),
                    "rr_phase_diff_deg": ph(c[0, 0] / lc[0, 0]),
                    "opp_hfss_db": max(db(c[0, 1]), db(c[1, 0])),
                    "opp_library_db": max(db(lc[0, 1]), db(lc[1, 0])),
                }
            )
    shifts = {f"{w}_{t}": shift_fit(h, lib, w, t) for w in ("amp", "opp") for t in (0.0, 30.0, None)}
    report = {
        "labels": ["MODEL_ONLY", "HFSS_COMPARISON_ONLY_NOT_CALIBRATED", "FIRST_HFSS_CHECK_OF_THE_RIS_CELL"],
        "plane_offset_mm": d_mm,
        "plane_fit_max_residual_deg": resid,
        "cross_sign_mean_abs_drr": sign,
        "summary_vs_library": summary,
        "candidate_offsets_not_applied": proposed_offsets(comps),
        "table": table,
        "at_10p5_ghz": at105,
        "frequency_shift_fits": shifts,
        "points_vs_library": [asdict(c) for c in comps],
    }
    (ROOT / "results/phase20_r4/hfss_check.json").write_text(
        json.dumps(report, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    )
    lines = [
        "| theta | f [GHz] | HFSS min/opp [dB] | HFSS steps | library min/opp [dB] | library steps | mom mesh 90 steps |",
        "|---|---|---|---|---|---|---|",
    ]
    for t in table:
        a, b, c = t["hfss"], t["library"], t["mom_mesh90"]
        lines.append(
            f"| {t['theta_deg']:g} | {t['f_ghz']:g} | {a['min_rr_db']:.2f} / {a['opp_db']:.1f} | {' '.join(f'{s:.0f}' for s in a['steps_deg'])} "
            f"| {b['min_rr_db']:.2f} / {b['opp_db']:.1f} | {' '.join(f'{s:.0f}' for s in b['steps_deg'])} | {' '.join(f'{s:.0f}' for s in c['steps_deg'])} |"
        )
    (ROOT / "results/phase20_r4/HFSS_CHECK.md").write_text("\n".join(lines) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in (
                    "plane_offset_mm",
                    "plane_fit_max_residual_deg",
                    "cross_sign_mean_abs_drr",
                    "summary_vs_library",
                    "frequency_shift_fits",
                )
            },
            indent=1,
        )
    )
    print("\n".join(lines))
    for r in at105:
        print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "hfss_results" / "phase20_r3_cp_cell")
