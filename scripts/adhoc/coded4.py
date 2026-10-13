import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
spec = importlib.util.spec_from_file_location("pb", "examples/make_paper_board_26ghz.py")
pb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pb)
from src.app.routers.hmsa import (  # noqa: E402
    HmsaSynthesizePortsRequest,
    TwoStateSourceRequest,
    _annular_library_base_config,
    _run_synthesize_ports_job,
)
from src.export.array_export import export_multi_port_array  # noqa: E402
from src.hologram.two_state_cells import LIBRARY_REFERENCE_OFFSET_M, patch_plane_gamma  # noqa: E402
from src.library_sweep import build_annular_slot_library_cells  # noqa: E402

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
GAP = float(sys.argv[2])
POL = sys.argv[3]
TH = float(sys.argv[4])
BMODE = sys.argv[5] if len(sys.argv) > 5 else "model"
F = 10.54e9
P = 6.4
b = 5.56
R = 1.922
r = 1.179
LAM = 299792458 / F * 1e3
S = float(sys.argv[6]) if len(sys.argv) > 6 else round(5 * LAM, 1)  # sub-aperture side, mm
pb.F_HZ = F
pb.PERIOD_MM = P
pb.PERIOD_M = P * 1e-3
BEAMS = [(TH, 45.0), (TH, 135.0), (TH, 225.0), (TH, 315.0)]
pb.BEAMS = BEAMS
ge = patch_plane_gamma(pb._gamma_top_empty(F), F, LIBRARY_REFERENCE_OFFSET_M)
body = dict(
    period_mm=P,
    eps_r=2.2,
    tan_delta=0.0009,
    thickness_mm=1.575,
    frequency_hz=F,
    accuracy="accurate",
    a_mm=b,
    b_mm=b,
    r_mm=r,
)
c = build_annular_slot_library_cells(_annular_library_base_config(body, placeholder_R_mm=R), b, b, r, [R], F)[0]
gs = patch_plane_gamma(complex(c["S11"][0]["real"], c["S11"][0]["imag"]), F, LIBRARY_REFERENCE_OFFSET_M)
ge5 = ge
ts = TwoStateSourceRequest(
    source="uploaded",
    uploaded_label_0="empty (bare substrate)",
    uploaded_value_0_real=ge5.real,
    uploaded_value_0_imag=ge5.imag,
    uploaded_label_1=f"{b} mm square, annular slot r {r} / R {R} mm",
    uploaded_value_1_real=gs.real,
    uploaded_value_1_imag=gs.imag,
    uploaded_geometry_0={"kind": "empty"},
    uploaded_geometry_1={"kind": "annular_slot", "a_mm": b, "b_mm": b, "R_mm": R, "r_mm": r},
).model_dump()
g = GAP / 2
regions = [[g, g + S, g, g + S], [-g - S, -g, g, g + S], [-g - S, -g, -g - S, -g], [g, g + S, -g - S, -g]]
beam = lambda t, p: {
    "theta_deg": t,
    "phi_deg": p,
    "weight_real": 1.0,
    "weight_imag": 0.0,
    "initial_phase_deg": 0.0,
    "desired_polarization": POL,
    "oam_charge": 0,
}
req = HmsaSynthesizePortsRequest(
    formulation="coded",
    frequency_hz=F,
    period_mm=P,
    layout={"preset": "custom", "regions_mm": regions, "gap_mm": GAP},
    ports=[{"beams": [beam(*BEAMS[i])]} for i in range(4)],
    port_eps_r=2.2,
    port_thickness_mm=1.575,
    port_tan_delta=0.0009,
    two_state=ts,
    beta_mode=BMODE,
).model_dump()
res = _run_synthesize_ports_job(pb._Job(f"coded4_gap{int(GAP)}"), req)
print("synth done", flush=True)
pats = [(p, pb._port_pattern(p, F)[0]) for p in res["ports"]]
mets = [pb._port_metrics(p, pat) for p, pat in pats]
# Fine-grid (0.5 deg theta x 1 deg phi) re-evaluation of each port's own pattern for SLL/XPD/pointing:
# the shared 1 x 5 deg grid is too coarse for a > 10 lambda0 sub-aperture's ~5 deg beam.
from src.pattern.vector_aperture_pattern import (  # noqa: E402
    VectorPatternResult,
    _vector_far_field,
    beam_polarization_metrics,
    circular_components,
)


def _chunked_pattern(x_m, y_m, jx, jy, green, dx, dy, n_theta=181, n_phi=360, rows=4):
    """compute_vector_pattern's own far field (_vector_far_field, unmodified) on a fine grid,
    evaluated a few theta rows at a time so a ~2000-cell aperture fits in memory."""
    theta = np.linspace(0.0, 90.0, n_theta)
    phi = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    et = np.empty(tg.shape, complex)
    ep = np.empty(tg.shape, complex)
    for i in range(0, n_theta, rows):
        et[i : i + rows], ep[i : i + rows] = _vector_far_field(
            tg[i : i + rows], pg[i : i + rows], x_m, y_m, jx, jy, F, green, dx, dy
        )
    er, el = circular_components(et, ep)
    power = np.abs(et) ** 2 + np.abs(ep) ** 2
    i, j = np.unravel_index(int(np.argmax(power)), power.shape)
    w = np.sin(np.radians(tg)) * np.radians(theta[1] - theta[0]) * (2 * np.pi / n_phi)
    directivity = 10 * np.log10(4 * np.pi * power.max() / np.sum(power * w))
    return VectorPatternResult(
        theta_deg=tg,
        phi_deg=pg,
        e_theta=et,
        e_phi=ep,
        e_r=er,
        e_l=el,
        power_db=10 * np.log10(np.maximum(power / power.max(), 1e-30)),
        main_beam_theta_deg=float(tg[i, j]),
        main_beam_phi_deg=float(pg[i, j]),
        sll_db=None,
        directivity_estimate_dbi=float(directivity),
    )


def _fine_metrics(port, k):
    x_m, y_m = np.asarray(port["x_m"]), np.asarray(port["y_m"])
    jx, jy = pb._port_pattern(port, F)[1]
    dx, dy = float(x_m[0, 1] - x_m[0, 0]), float(y_m[1, 0] - y_m[0, 0])
    green = pb.LayeredGreenFunction(pb.SUBSTRATE, F, backing="pec")
    fp = _chunked_pattern(x_m, y_m, jx, jy, green, dx, dy)
    t, ph = BEAMS[k]
    m = beam_polarization_metrics(fp, POL, t, ph, float(x_m.max() - x_m.min()) + dx, 299792458.0 / F)
    from src.pattern.vector_aperture_pattern import _direction_mask, copol_component, first_null_half_width_deg

    co = np.abs(copol_component(fp.e_theta, fp.e_phi, POL)) ** 2
    excl = first_null_half_width_deg(float(x_m.max() - x_m.min()) + dx, 299792458.0 / F, t)
    outside = np.where(_direction_mask(fp.theta_deg, fp.phi_deg, t, ph, excl), -1.0, co)
    si, sj = np.unravel_index(int(np.argmax(outside)), outside.shape)
    # Cuts through the beam peak along the board's own x and y directions (constant v / constant u
    # in direction-cosine space) -- a square aperture's first sidelobes sit on these lines.
    u0 = np.sin(np.radians(fp.main_beam_theta_deg)) * np.cos(np.radians(fp.main_beam_phi_deg))
    v0 = np.sin(np.radians(fp.main_beam_theta_deg)) * np.sin(np.radians(fp.main_beam_phi_deg))
    from src.pattern.vector_aperture_pattern import ORTHOGONAL_POLARIZATION as _OPOL

    cuts = {}
    for name, (uu, vv) in {
        "x": (np.linspace(-1, 1, 801), np.full(801, v0)),
        "y": (np.full(801, u0), np.linspace(-1, 1, 801)),
    }.items():
        ok = uu**2 + vv**2 < 1.0
        uu, vv = uu[ok], vv[ok]
        tt = np.degrees(np.arcsin(np.sqrt(uu**2 + vv**2)))
        pp = np.degrees(np.arctan2(vv, uu))
        et, ep = _vector_far_field(tt, pp, x_m, y_m, jx, jy, F, green, dx, dy)
        cuts[name] = (
            uu if name == "x" else vv,
            np.abs(copol_component(et, ep, POL)) ** 2,
            np.abs(copol_component(et, ep, _OPOL[POL])) ** 2,
        )
    req = pb._direction_cosines(np.array(t), np.array(ph))
    got = pb._direction_cosines(np.array(fp.main_beam_theta_deg), np.array(fp.main_beam_phi_deg))
    co_pk = float(co.max())
    cut_sll = {}
    for name, (_axis, cco, _) in cuts.items():
        db = 10 * np.log10(np.maximum(cco / co_pk, 1e-12))
        k0 = int(np.argmax(db))
        # walk down from the peak to the first local minimum on each side, then take the max beyond
        lo = k0
        while lo > 0 and db[lo - 1] < db[lo]:
            lo -= 1
        hi = k0
        while hi < len(db) - 1 and db[hi + 1] < db[hi]:
            hi += 1
        side = np.concatenate([db[:lo], db[hi + 1 :]])
        cut_sll[name] = float(side.max()) if side.size else None
    return (
        fp,
        cuts,
        {
            "grid": "0.5 deg theta x 1 deg phi",
            "sll_cut_along_x_db": cut_sll["x"],
            "sll_cut_along_y_db": cut_sll["y"],
            "peak_theta_deg": float(fp.main_beam_theta_deg),
            "peak_phi_deg": float(fp.main_beam_phi_deg),
            "pointing_error_deg": float(np.degrees(np.arccos(np.clip(req @ got, -1, 1)))),
            "directivity_dbi": float(fp.directivity_estimate_dbi),
            "sll_copol_db": m.sll_copol_db,
            "sll_direction_deg": [float(fp.theta_deg[si, sj]), float(fp.phi_deg[si, sj])],
            "main_lobe_exclusion_deg": float(excl),
            "point_xpd_db": m.point_xpd_db,
            "integrated_xpd_db": m.integrated_xpd_db,
            "cross_peak_db": m.cross_peak_db,
        },
    )


fine = [_fine_metrics(p, k) for k, (p, _) in enumerate(pats)]
for met, (_, _, fm) in zip(mets, fine, strict=True):
    met["fine_grid"] = fm


def _cut_figure(path):
    fig, axes = pb.plt.subplots(4, 2, figsize=(11, 13), sharey=True)
    for k, (_, cuts, fm) in enumerate(fine):
        co_pk = max(float(c[1].max()) for c in cuts.values())
        for col, name in enumerate(("x", "y")):
            ax = axes[k, col]
            axis, co, cx = cuts[name]
            ax.plot(axis, 10 * np.log10(np.maximum(co / co_pk, 1e-6)), label=POL.upper() + " (co)")
            ax.plot(axis, 10 * np.log10(np.maximum(cx / co_pk, 1e-6)), label="cross", alpha=0.8)
            sll = fm[f"sll_cut_along_{name}_db"]
            other = "v" if name == "x" else "u"
            ax.set_title(
                f"Port {k + 1}: cut along board {name} through the beam ({other} fixed), SLL {sll:.1f} dB", fontsize=9
            )
            ax.set_xlim(-1, 1)
            ax.set_ylim(-45, 2)
            ax.grid(alpha=0.3)
            ax.set_xlabel("u = sin(theta) cos(phi)" if name == "x" else "v = sin(theta) sin(phi)")
            ax.legend(fontsize=7)
        axes[k, 0].set_ylabel("dB rel. co-pol peak")
    fig.suptitle(f"4-port 1-bit HMSA, {S:g} mm sub-apertures, gap {GAP:g} mm, 10.54 GHz (model prediction)")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    pb.plt.close(fig)


_cut_figure(OUT / "port_cuts.png")
pb.figure_patterns(pats, OUT / "port_patterns.png")
pb.figure_layout(res, OUT / "board_layout.png")
feeds = [(p["feed_x_mm"], p["feed_y_mm"]) for p in res["ports"]]
rep = {
    "frequency_GHz": 10.54,
    "cell_mm": P,
    "port_side_mm": S,
    "gap_mm": GAP,
    "polarization": POL,
    "beta_mode": BMODE,
    "beams": BEAMS,
    "cell_pair": {
        "empty_patch_plane_deg": float(np.degrees(np.angle(ge))),
        "slotted_patch_plane_deg": float(np.degrees(np.angle(gs))),
        "slotted": {"b_mm": b, "R_mm": R, "r_mm": r},
    },
    "two_state": res.get("two_state"),
    "isolation": res["isolation"],
    "board_mm": [res["board_width_mm"], res["board_height_mm"]],
    "cells": [len(res["lw_grid_m"][0]), len(res["lw_grid_m"])],
    "feeds_mm": feeds,
    "ports": mets,
}
json.dump(rep, open(OUT / "design_report.json", "w"), indent=1, default=float)
print(json.dumps(rep, indent=1, default=float)[:5000], flush=True)
exp = export_multi_port_array(
    res["lw_grid_m"],
    P * 1e-3,
    P * 1e-3,
    pb.SUBSTRATE,
    pb.METAL,
    OUT,
    include_step=True,
    feeds_mm=feeds,
    ports_info=[(k + 1, x, y) for k, (x, y) in enumerate(feeds)],
)
print("step", exp.file_size_bytes, flush=True)
