"""Phase 19 Part C: design the 28 GHz nylon dual-beam (+-30 deg) lens.

    .venv/bin/python scripts/research/phase19_design.py

D = 150 mm octagon (across flats), pitch 5 mm, F = 1.2 D, 0.6 mm plate with through-hole pockets, total 25 mm.
Pass 1 uses the feed-tapered aperture amplitude; pass 2 re-synthesises with the actual |T| of the pass-1 cells.
Writes results/transmitarray_dual30/{layout.json, summary.json}.  Model-only (HFSS: NOT_AVAILABLE_FROM_CURRENT_DATA).
"""

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.transmitarray.cells import CellFamilyModel  # noqa: E402
from src.transmitarray.design_inputs import (  # noqa: E402
    BEAM_DEG,
    EDGE_TAPER_DB,
    FREQ_HZ,
    make_layout,
)
from src.transmitarray.pipeline import DesignConfig, predict_pattern, run_design  # noqa: E402

OUT = ROOT / "results" / "transmitarray_dual30"
LIB = ROOT / "results" / "nylon_cells_28ghz" / "library_plate0.6mm_through.csv"
PLATE_M, TOTAL_M, EPS, TAND, MIN_FEATURE_M = 0.6e-3, 25.0e-3, 2.45, 0.0185, 0.6e-3
D_SEAM_DEG, CAP_FILL, N_HARM = 230.0, 0.8044, 8


def config():
    lay = make_layout()
    model = CellFamilyModel(lay.pitch_m, TOTAL_M, PLATE_M, EPS, TAND, True, FREQ_HZ, N_HARM)
    return DesignConfig(FREQ_HZ, lay, BEAM_DEG, EDGE_TAPER_DB, model, MIN_FEATURE_M, D_SEAM_DEG, CAP_FILL, str(LIB))


def cells_record(cfg, res, tau):
    syn = res["syn"]
    out = []
    for i, c in enumerate(res["cells"]):
        j = c["J"]
        out.append(
            dict(
                x_mm=float(syn["x"][i] * 1e3),
                y_mm=float(syn["y"][i] * 1e3),
                theta_deg=float(syn["theta"][i]),
                phi_deg=float(syn["phi"][i]),
                kind=c["kind"],
                radius_mm=float(c["R"] * 1e3),
                residual_deg=float(c["residual_deg"]),
                n_solves=int(c["n_solves"]),
                target_deg=float(((res["target_deg"][i] + 180) % 360) - 180),
                power_transmission=float(tau[i]),
                J=[[float(j[a, b].real), float(j[a, b].imag)] for a in range(2) for b in range(2)],
            )
        )
    return out


def summarize(cfg, res, metrics, tag):
    r = np.array([c["residual_deg"] for c in res["cells"]])
    kinds = [c["kind"] for c in res["cells"]]
    return dict(
        pass_name=tag,
        psi0_deg=res["psi0_deg"],
        n_cells=len(kinds),
        n_pillar=kinds.count("pillar"),
        n_hole=kinds.count("hole"),
        residual_deg=dict(
            max_abs=float(np.abs(r).max()), rms=float(np.sqrt(np.mean(r**2))), n_over_1deg=int(np.sum(np.abs(r) > 1.0))
        ),
        symmetry_dev_deg=res["syn"]["symmetry_dev_deg"],
        synthesis=res["syn"]["gs"],
        pattern=metrics,
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = config()
    res1 = run_design(cfg)
    ff1, m1, ex1, ey1, tau1 = predict_pattern(cfg, res1)
    print("pass1", json.dumps(summarize(cfg, res1, m1, "pass1"), indent=1), flush=True)
    amp = np.sqrt(tau1)
    res2 = run_design(cfg, cell_amp=amp, psi0_deg=res1["psi0_deg"])
    ff2, m2, ex2, ey2, tau2 = predict_pattern(cfg, res2)
    s2 = summarize(cfg, res2, m2, "pass2")
    print("pass2", json.dumps(s2, indent=1), flush=True)
    final = (
        (res2, m2, tau2) if m2["beam_phi0"]["peak_gain_dbi"] >= m1["beam_phi0"]["peak_gain_dbi"] else (res1, m1, tau1)
    )
    chosen = "pass2" if final[0] is res2 else "pass1"
    (OUT / "layout.json").write_text(
        json.dumps(dict(chosen=chosen, cells=cells_record(cfg, final[0], final[2])), indent=1)
    )
    lay = cfg.layout
    summary = dict(
        inputs=dict(
            freq_hz=FREQ_HZ,
            d_flats_mm=lay.d_flats_m * 1e3,
            pitch_mm=lay.pitch_m * 1e3,
            focal_mm=lay.focal_m * 1e3,
            beam_deg=BEAM_DEG,
            edge_taper_db=EDGE_TAPER_DB,
            q=final[0]["syn"]["q"],
            plate_mm=PLATE_M * 1e3,
            total_mm=TOTAL_M * 1e3,
            through_hole=True,
            d_seam_deg=D_SEAM_DEG,
            cap_fill=CAP_FILL,
            n_harm=N_HARM,
            eps_r=EPS,
            tan_delta=TAND,
        ),
        pass1=summarize(cfg, res1, m1, "pass1"),
        pass2=s2,
        chosen=chosen,
        hfss_comparison="NOT_AVAILABLE_FROM_CURRENT_DATA",
    )
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print("done", chosen, flush=True)


if __name__ == "__main__":
    main()
