import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.app.routers.hmsa import _l_fractions_of_period, _dispersion_summary, _port_cell_to_m
from src.basis_functions import MeshSettings
from src.dispersion.mom_scaffold import square_patch_unit_cell
from src.dispersion.sheet_model import build_dispersion_table_from_sheet_model
from src.dispersion.analytic_tm0 import tm0_beta
from src.materials import Metal, Substrate
from src.phase8b.object_wave import Beam
from src.hologram.hmsa_synthesis import synthesize_formulation_a
from src.hologram.multi_port_synthesis import synthesize_ports, combined_board_lw_grid_m
from src.hologram.ports import layout_preset_1xn
from src.pattern.vector_aperture_pattern import aperture_current_from_scalar, compute_vector_pattern
from src.green_function import LayeredGreenFunction
from src.export.array_export import export_multi_port_array
from src.app.plots import patch_layout_png

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
F = 10.54e9
P = 4e-3
PORT = float(sys.argv[2]) * 1e-3
GAP = float(sys.argv[3]) * 1e-3
TH = 25.0
BEAMS = [(TH, 180.0), (TH, 90.0), (TH, 0.0)]
sub = Substrate("user", 2.2, 0.0009, 1.575e-3)
metal = Metal()
t = time.time()
mc = 32
ls = [P * f for f in _l_fractions_of_period(mc)]
table = build_dispersion_table_from_sheet_model(
    lambda l: square_patch_unit_cell(P, sub, metal, l),
    ls,
    F,
    MeshSettings(cells_x=mc, cells_y=mc, floquet_order=mc),
    period_m=P,
)
S = _dispersion_summary(table)
print("dispersion", round(time.time() - t, 1), "s", S, flush=True)
regions = layout_preset_1xn(3, PORT + GAP, gap_m=GAP)
pb = {
    r.port_id: [
        Beam(theta_deg=BEAMS[r.port_id][0], phi_deg=BEAMS[r.port_id][1], weight=1 + 0j, desired_polarization="rhcp")
    ]
    for r in regions
}
k_swg = tm0_beta(F, 2.2, 1.575e-3)
multi = synthesize_ports(
    regions,
    pb,
    lambda ref, beams, feed: synthesize_formulation_a(
        ref, beams, F, S["x_avg_ohm"], S["modulation_index"], table, P, rescale_to_beta_bar=True, feed_xy_m=feed
    ),
    P,
    k_swg=k_swg,
    min_lambda_sw=1.0,
)
green = LayeredGreenFunction(sub, F, backing="pec")
k0 = 2 * np.pi * F / 299792458.0


def ff(x, y, jx, jy, th, ph):
    th = np.radians(th)
    ph = np.radians(ph)
    kx = k0 * np.sin(th) * np.cos(ph)
    ky = k0 * np.sin(th) * np.sin(ph)
    Pm = np.exp(1j * (np.multiply.outer(kx, x.ravel()) + np.multiply.outer(ky, y.ravel())))
    Jx = Pm @ jx.ravel()
    Jy = Pm @ jy.ravel()
    Et = np.cos(th) * (Jx * np.cos(ph) + Jy * np.sin(ph))
    Ep = -Jx * np.sin(ph) + Jy * np.cos(ph)
    return (Et + 1j * Ep) / np.sqrt(2), (Et - 1j * Ep) / np.sqrt(2)


rep = {
    "frequency_GHz": 10.54,
    "period_mm": 4,
    "port_mm": PORT * 1e3,
    "gap_mm": GAP * 1e3,
    "x_avg_ohm": S["x_avg_ohm"],
    "modulation_index": S["modulation_index"],
    "ports": [],
    "isolation": [
        {"ports": [f.port_a + 1, f.port_b + 1], "gap_lambda_sw": round(f.gap_lambda_sw, 2), "flagged": f.flagged}
        for f in multi.isolation_flags
    ],
}
fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, port in zip(axes, multi.ports):
    d = port.design
    jx, jy = aperture_current_from_scalar(d.x_m, d.y_m, d.aperture_field, feed_xy_m=port.feed_xy_m)
    pat = compute_vector_pattern(d.x_m, d.y_m, jx, jy, F, green, P, P)
    THg, PHg = np.meshgrid(np.arange(0, 90.1, 0.5), np.arange(0, 360, 1.0), indexing="ij")
    R, L = ff(d.x_m, d.y_m, jx, jy, THg, PHg)
    Pw = abs(R) ** 2 + abs(L) ** 2
    i = np.unravel_index(Pw.argmax(), Pw.shape)
    t0, p0 = BEAMS[port.port_id]
    r, l = ff(d.x_m, d.y_m, jx, jy, np.array([t0]), np.array([p0]))
    rep["ports"].append(
        {
            "port": port.port_id + 1,
            "feed_mm": [round(port.feed_xy_m[0] * 1e3, 3), round(port.feed_xy_m[1] * 1e3, 3)],
            "region_mm": [
                round(v * 1e3, 3)
                for v in (port.region.x_min_m, port.region.x_max_m, port.region.y_min_m, port.region.y_max_m)
            ],
            "cells": list(d.x_m.shape),
            "requested": [t0, p0],
            "unreachable_pct": round(100 * d.unreachable_fraction, 1),
            "model_peak": [round(float(pat.main_beam_theta_deg), 1), round(float(pat.main_beam_phi_deg), 1)],
            "model_directivity_dBi": round(float(pat.directivity_estimate_dbi), 1),
            "freespace_peak": [float(THg[i]), float(PHg[i])],
            "freespace_XPD_dB_at_request": round(float(20 * np.log10(abs(r[0]) / abs(l[0]))), 1),
        }
    )
    u = np.sin(np.radians(pat.theta_deg)) * np.cos(np.radians(pat.phi_deg))
    v = np.sin(np.radians(pat.theta_deg)) * np.sin(np.radians(pat.phi_deg))
    cs = ax.contourf(u, v, np.maximum(pat.power_db, -40), levels=np.linspace(-40, 0, 41), cmap="viridis")
    ax.plot(
        [np.sin(np.radians(t0)) * np.cos(np.radians(p0))],
        [np.sin(np.radians(t0)) * np.sin(np.radians(p0))],
        "r+",
        ms=12,
        mew=2,
    )
    ax.set_aspect("equal")
    ax.set_title(
        f"Port {port.port_id + 1}: RHCP beam ({t0:g}°, {p0:g}°)\npeak ({pat.main_beam_theta_deg:.0f}°, {pat.main_beam_phi_deg:.0f}°), D = {pat.directivity_estimate_dbi:.1f} dBi",
        fontsize=9,
    )
    ax.set_xlabel("u")
    ax.set_ylabel("v")
fig.colorbar(cs, ax=axes, label="dB")
fig.savefig(OUT / "port_patterns_uv.png", dpi=150)
plt.close(fig)
board = combined_board_lw_grid_m(multi, _port_cell_to_m("A", P * 1e3))
FEEDS = [(p.feed_xy_m[0], p.feed_xy_m[1]) for p in multi.ports]
cleared = 0
for i, row in enumerate(board.lw_grid_m):
    for j, (lm, wm) in enumerate(row):
        cx, cy = board.cell_center_m(i, j)
        if lm > 0 and any(abs(cx - fx) < 4e-3 and abs(cy - fy) < 4e-3 for fx, fy in FEEDS):
            row[j] = (0.0, 0.0)
            cleared += 1
print("feed keep-out cleared", cleared, flush=True)
cells = []
for i, row in enumerate(board.lw_grid_m):
    for j, (lm, wm) in enumerate(row):
        if lm > 0 and wm > 0:
            cx, cy = board.cell_center_m(i, j)
            cells.append((cx * 1e3, cy * 1e3, lm * 1e3, wm * 1e3))
ny, nx = len(board.lw_grid_m), len(board.lw_grid_m[0])
patch_layout_png(
    cells, P * 1e3, nx, ny, OUT / "patch_layout.png", title=f"3-port HMSA, {nx} x {ny} cells, 4 mm", feed_xy_mm=None
)
feeds = [(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports]
export_multi_port_array(
    board.lw_grid_m,
    board.period_x_m,
    board.period_y_m,
    sub,
    metal,
    OUT,
    include_step=True,
    feeds_mm=feeds,
    ports_info=[(k + 1, x, y) for k, (x, y) in enumerate(feeds)],
)
rep.update(
    board_mm=[round(board.width_m * 1e3, 3), round(board.height_m * 1e3, 3)],
    board_cells=[nx, ny],
    board_origin_mm=[round(board.x_min_m * 1e3, 3), round(board.y_min_m * 1e3, 3)],
    pitch_mm=[round(board.period_x_m * 1e3, 4), round(board.period_y_m * 1e3, 4)],
    n_patches=len(cells),
    feed_keepout_mm="8 x 8 (4 cells cleared per feed)",
    patch_sizes_mm=sorted({round(c[2], 3) for c in cells}),
)
json.dump(rep, open(OUT / "design_report.json", "w"), indent=1)
print(json.dumps(rep, indent=1))
