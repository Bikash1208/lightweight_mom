"""Submission figures for Phase 15: vortex vs pencil (normal vs compensated beam), an intentional
OAM vortex, a Mishra-style multi-port board, and its STEP export with HFSS README.

Run from the project root (needs the `cad` extra for the STEP file):
    .venv/bin/python examples/make_submission_figures.py

Writes into docs/submission/: one PNG per figure, a `<name>.caption.md` listing the inputs and the
measured numbers, `submission_metrics.json`, and `multiport_step/` (STEP + HFSS_README.md).

Inputs (all listed again in each caption): the tracked HFSS phase-only field
`inputs/reference_fields/efl_phase_EzZ_example.fld` at 10.54 GHz; 8 mm cells; RT/Duroid 5880
(eps_r 2.2, tan_delta 0.0009, 1.575 mm) on a PEC backing; a fast-mesh square-patch library built
here with the app's own settings.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.app.plots import vortex_phase_ring_png  # noqa: E402
from src.app.routers.hmsa import AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY, _l_fractions_of_period, _mask_feed_region  # noqa: E402
from src.app.library_query import library_candidates_at_frequency  # noqa: E402
from src.dispersion.analytic_tm0 import tm0_beta  # noqa: E402
from src.export.array_export import export_multi_port_array  # noqa: E402
from src.fields.fld_parser import read_fld  # noqa: E402
from src.fields.reference_field_grid import (  # noqa: E402
    feed_position_from_nan_mask,
    fld_to_reference_field_grid,
    resample_to_cell_grid,
)
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.hologram.hmsa_synthesis import synthesize_formulation_b  # noqa: E402
from src.hologram.multi_port_synthesis import combined_board_lw_grid_m, port_reference_field, synthesize_ports  # noqa: E402
from src.hologram.ports import layout_preset_2x2_plus_1_large  # noqa: E402
from src.library_sweep import build_square_library_cells, write_square_library_json  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.pattern.vector_aperture_pattern import (  # noqa: E402
    aperture_current_from_scalar,
    aperture_efficiency_estimate,
    compute_vector_pattern,
    copol_label,
    measure_vortex_winding,
)
from src.phase8b.object_wave import Beam  # noqa: E402
from src.unit_cell_library import UnitCellLibrary  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "submission"
FLD = ROOT / "inputs" / "reference_fields" / "efl_phase_EzZ_example.fld"
F_HZ = 10.54e9
PERIOD_M = 8.0e-3
EPS_R, TAN_D, THICK_M = 2.2, 0.0009, 1.575e-3
BEAM_THETA, BEAM_PHI = 20.0, 0.0
MULTIPORT_BEAMS = [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
SUBSTRATE = Substrate(name="user", eps_r=EPS_R, tan_delta=TAN_D, thickness=THICK_M)


def _library_candidates(work_dir):
    accuracy = "fast"
    mesh_cells = AUTO_LIBRARY_MESH_CELLS_BY_ACCURACY[accuracy]
    l_values = tuple(round(PERIOD_M * 1e3 * f, 6) for f in _l_fractions_of_period(mesh_cells))
    base = {
        "period_x": PERIOD_M, "period_y": PERIOD_M,
        "substrate": {"name": "user", "eps_r": EPS_R, "tan_delta": TAN_D, "thickness": THICK_M, "backing": "pec"},
        "geometry": {"type": "rectangle", "length": l_values[0] * 1e-3, "width": l_values[0] * 1e-3},
        "excitation": {"frequency": F_HZ},
        "mesh": {"cells_x": mesh_cells, "cells_y": mesh_cells, "floquet_order": mesh_cells},
        "reference_planes": {"top_air_thickness": 5.0e-3},
    }
    cells = build_square_library_cells(base, l_values, F_HZ)
    path = work_dir / "square_library.json"
    write_square_library_json(cells, path)
    lib = UnitCellLibrary.from_json(path)
    return library_candidates_at_frequency(lib, "submission_square_library", F_HZ), l_values


def _reference():
    raw = fld_to_reference_field_grid(read_fld(FLD))
    raw_nan = np.isnan(raw.phase())
    masked, _ = _mask_feed_region(raw)
    x0, y0 = feed_position_from_nan_mask(raw.x_m, raw.y_m, raw_nan)
    return resample_to_cell_grid(masked, PERIOD_M), (float(x0), float(y0))


def _pattern(x_m, y_m, aperture_field, feed_xy_m, green):
    jx, jy = aperture_current_from_scalar(x_m, y_m, aperture_field, feed_xy_m=feed_xy_m)
    pat = compute_vector_pattern(x_m, y_m, jx, jy, F_HZ, green, PERIOD_M, PERIOD_M)
    return pat, jx, jy


def _nearest(pat, theta_deg, phi_deg):
    i = int(np.argmin(np.abs(pat.theta_deg[:, 0] - theta_deg)))
    j = int(np.argmin(np.abs(((pat.phi_deg[0, :] - phi_deg + 180.0) % 360.0) - 180.0)))
    return i, j


def _axial_ratio_db(e_r, e_l):
    major = np.abs(e_r) + np.abs(e_l)
    minor = np.maximum(np.abs(np.abs(e_r) - np.abs(e_l)), 1e-12)
    return 20.0 * np.log10(major / minor)


def _measure(pat, jx, jy, x_m, y_m, green, desired):
    peak = float(np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2))
    if desired == "natural":
        # No copol component exists for a natural beam: measure both circular components and keep the
        # one with the larger winding (the vortex's own handedness), reporting both below.
        both = {h: measure_vortex_winding(x_m, y_m, jx, jy, F_HZ, green, PERIOD_M, PERIOD_M, BEAM_THETA, BEAM_PHI, h,
                                          pattern_peak_power=peak) for h in ("rhcp", "lhcp")}
        vortex = max(both.values(), key=lambda m: abs(m.winding_number))
    else:
        vortex = measure_vortex_winding(
            x_m, y_m, jx, jy, F_HZ, green, PERIOD_M, PERIOD_M, BEAM_THETA, BEAM_PHI, desired,
            pattern_peak_power=peak,
        )
    i, j = _nearest(pat, BEAM_THETA, BEAM_PHI)
    e_r, e_l = pat.e_r[i, j], pat.e_l[i, j]
    ar_db = float(_axial_ratio_db(np.array([e_r]), np.array([e_l]))[0])
    result = {
        "null_depth_db": float(vortex.null_depth_db),
        "winding_number": float(vortex.winding_number),
        "axial_ratio_db_at_beam": ar_db,
        "e_r_db_at_beam": float(20.0 * np.log10(max(abs(e_r), 1e-12) / np.sqrt(peak))),
        "e_l_db_at_beam": float(20.0 * np.log10(max(abs(e_l), 1e-12) / np.sqrt(peak))),
        "main_beam_theta_deg": float(pat.main_beam_theta_deg),
        "main_beam_phi_deg": float(pat.main_beam_phi_deg),
        "directivity_dbi": float(pat.directivity_estimate_dbi),
    }
    if desired == "natural":
        result["winding_rhcp_lhcp"] = [float(both["rhcp"].winding_number), float(both["lhcp"].winding_number)]
    if desired in ("rhcp", "lhcp"):
        co_key = copol_label(desired)
        co = pat.e_r[i, j] if co_key == "E_R" else pat.e_l[i, j]
        cross = pat.e_l[i, j] if co_key == "E_R" else pat.e_r[i, j]
        result["copol_component"] = co_key
        result["xpd_db_at_beam"] = float(20.0 * np.log10(max(abs(co), 1e-12) / max(abs(cross), 1e-12)))
    return result, vortex


def _uv_panel(ax, pat, title, marker=True):
    u = np.sin(np.radians(pat.theta_deg)) * np.cos(np.radians(pat.phi_deg))
    v = np.sin(np.radians(pat.theta_deg)) * np.sin(np.radians(pat.phi_deg))
    levels = np.linspace(-40.0, 0.0, 41)
    cs = ax.contourf(u, v, np.maximum(pat.power_db, -40.0), levels=levels, cmap="viridis")
    if marker:
        bu = np.sin(np.radians(BEAM_THETA)) * np.cos(np.radians(BEAM_PHI))
        bv = np.sin(np.radians(BEAM_THETA)) * np.sin(np.radians(BEAM_PHI))
        ax.plot([bu], [bv], "r+", ms=12, mew=2)
    ax.set_aspect("equal")
    ax.set_xlabel("u")
    ax.set_ylabel("v")
    ax.set_title(title, fontsize=10)
    return cs


def _circular_cut_figure(pat, title, path):
    i0, _ = _nearest(pat, 0.0, 0.0)
    theta = pat.theta_deg[:, 0]
    peak = float(np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2))
    e_r = pat.e_r[:, 0]
    e_l = pat.e_l[:, 0]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 6.0), sharex=True)
    ax1.plot(theta, 10 * np.log10(np.maximum(np.abs(e_r) ** 2 / peak, 1e-12)), label="E_R")
    ax1.plot(theta, 10 * np.log10(np.maximum(np.abs(e_l) ** 2 / peak, 1e-12)), label="E_L")
    ax1.axvline(BEAM_THETA, color="r", ls=":", lw=1)
    ax1.set_ylabel("power (dB, peak = 0)")
    ax1.set_ylim(-45, 2)
    ax1.legend()
    ax1.set_title(f"{title}: circular cuts at phi = 0 deg", fontsize=10)
    ax2.plot(theta, _axial_ratio_db(e_r, e_l), color="k")
    ax2.axvline(BEAM_THETA, color="r", ls=":", lw=1)
    ax2.set_ylabel("axial ratio (dB)")
    ax2.set_xlabel("theta (deg)")
    ax2.set_ylim(0, 40)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _write_caption(name, text, inputs, numbers):
    lines = [f"# {name}", "", text, "", "## Inputs", ""] + [f"- {k}: {v}" for k, v in inputs.items()]
    lines += ["", "## Measured", ""] + [f"- {k}: {v}" for k, v in numbers.items()]
    (OUT / f"{name}.caption.md").write_text("\n".join(lines) + "\n")


def _single_beam_figures(reference, feed_xy_m, candidates, green, metrics):
    common_inputs = {
        "field": "inputs/reference_fields/efl_phase_EzZ_example.fld (phase-only HFSS field, 10.54 GHz)",
        "feed (x0, y0) mm": f"({feed_xy_m[0] * 1e3:.1f}, {feed_xy_m[1] * 1e3:.1f})",
        "cells": f"{PERIOD_M * 1e3:g} mm pitch, {reference.x_m.shape[1]} x {reference.x_m.shape[0]} reference grid",
        "substrate": f"RT/Duroid 5880, eps_r {EPS_R}, tan_delta {TAN_D}, {THICK_M * 1e3:g} mm, PEC backing",
        "method": "Formulation B (phase-difference hologram) with the square-patch library (fast mesh)",
        "scripts": "examples/make_submission_figures.py",
    }
    cases = [
        ("vortex_natural", Beam(theta_deg=BEAM_THETA, phi_deg=BEAM_PHI, weight=1.0 + 0j, desired_polarization="natural"),
         "Normal (natural, no compensation) beam at theta=20 deg, phi=0 deg: the monopole's radial surface current "
         "is not compensated, so the aperture radiates a vortex (null on boresight, nonzero OAM winding)."),
        ("pencil_rhcp", Beam(theta_deg=BEAM_THETA, phi_deg=BEAM_PHI, weight=1.0 + 0j, desired_polarization="rhcp"),
         "The same beam with per-beam RHCP compensation: a pencil beam, co-polar with E_R at boresight."),
        ("vortex_l1", Beam(theta_deg=BEAM_THETA, phi_deg=BEAM_PHI, weight=1.0 + 0j, desired_polarization="rhcp", oam_charge=1),
         "An intentionally requested OAM vortex, charge l=+1 (RHCP request). It reproduces the natural beam: the "
         "uncompensated monopole current already carries winding +1, so the natural beam is an l=+1 vortex."),
        ("vortex_l2_limitation", Beam(theta_deg=BEAM_THETA, phi_deg=BEAM_PHI, weight=1.0 + 0j, desired_polarization="rhcp", oam_charge=2),
         "LIMITATION, not a clean charge-2 vortex: a requested l=+2 (RHCP) puts the boresight power into the opposite "
         "(cross) circular component with zero winding, so the charge is not realised as requested in this model. "
         "Reported as measured."),
    ]
    for name, beam, text in cases:
        design = synthesize_formulation_b(reference, [beam], F_HZ, candidates, PERIOD_M, feed_xy_m=feed_xy_m)
        pat, jx, jy = _pattern(design.x_m, design.y_m, design.aperture_field, feed_xy_m, green)
        desired = beam.desired_polarization
        numbers, vortex = _measure(pat, jx, jy, design.x_m, design.y_m, green, desired)
        metrics[name] = numbers
        uv_path = OUT / f"{name}_uv.png"
        fig, ax = plt.subplots(figsize=(5.2, 4.6))
        cs = _uv_panel(ax, pat, f"{name}: co-pol power (dB)")
        fig.colorbar(cs, ax=ax, label="dB")
        fig.tight_layout()
        fig.savefig(uv_path, dpi=150)
        plt.close(fig)
        _circular_cut_figure(pat, name, OUT / f"{name}_circular_cuts.png")
        vortex_phase_ring_png(vortex.ring_alpha_deg, vortex.ring_phase_deg, vortex.winding_number,
                              OUT / f"{name}_phase_ring.png", title=f"{name}: phase around the beam axis")
        _write_caption(
            name, text,
            {**common_inputs, "beam": f"theta {BEAM_THETA:g} deg, phi {BEAM_PHI:g} deg, desired_polarization {desired}, "
                                       f"oam_charge {beam.oam_charge}"},
            {k: (round(v, 3) if isinstance(v, float) else v) for k, v in numbers.items()},
        )


def _multiport_figures(candidates, green, metrics, work_dir):
    regions = layout_preset_2x2_plus_1_large(60e-3, gap_m=0.0)
    port_beams = {}
    for region in regions:
        if region.port_id < 4:
            theta, phi = MULTIPORT_BEAMS[region.port_id]
            port_beams[region.port_id] = [Beam(theta_deg=theta, phi_deg=phi, weight=1.0 + 0j, desired_polarization="rhcp")]
        else:
            port_beams[region.port_id] = [
                Beam(theta_deg=t, phi_deg=p, weight=1.0 + 0j, desired_polarization="rhcp") for t, p in MULTIPORT_BEAMS
            ]
    k_swg = tm0_beta(F_HZ, EPS_R, THICK_M)

    def per_port(ref, beams, feed_xy_m):
        return synthesize_formulation_b(ref, beams, F_HZ, candidates, PERIOD_M, feed_xy_m=feed_xy_m)

    multi = synthesize_ports(regions, port_beams, per_port, PERIOD_M, k_swg=k_swg, min_lambda_sw=1.0)
    board = combined_board_lw_grid_m(
        multi,
        lambda c: None if c.geometry_param is None else (
            c.geometry_param * 1e-3, (c.geometry_param_w if c.geometry_param_w is not None else c.geometry_param) * 1e-3
        ),
    )

    per_port_metrics = []
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, port in zip(axes.ravel(), multi.ports, strict=False):
        region = port.region
        ref_x, ref_y = port.design.x_m, port.design.y_m
        pat, _, _ = _pattern(ref_x, ref_y, port.design.aperture_field, port.feed_xy_m, green)
        eff, dmax = aperture_efficiency_estimate(pat.directivity_estimate_dbi, region.area_m2, F_HZ)
        per_port_metrics.append({
            "port": port.port_id + 1,
            "beams": MULTIPORT_BEAMS[port.port_id] if port.port_id < 4 else MULTIPORT_BEAMS,
            "feed_mm": [round(port.feed_xy_m[0] * 1e3, 1), round(port.feed_xy_m[1] * 1e3, 1)],
            "aperture_m2": round(region.area_m2, 6),
            "directivity_dbi": round(pat.directivity_estimate_dbi, 2),
            "gain_efficiency_fraction": round(eff, 3),
            "max_directivity_dbi": round(dmax, 2),
            "main_beam_theta_deg": round(pat.main_beam_theta_deg, 1),
            "main_beam_phi_deg": round(pat.main_beam_phi_deg, 1),
        })
        u = np.sin(np.radians(pat.theta_deg)) * np.cos(np.radians(pat.phi_deg))
        v = np.sin(np.radians(pat.theta_deg)) * np.sin(np.radians(pat.phi_deg))
        ax.contourf(u, v, np.maximum(pat.power_db, -40.0), levels=np.linspace(-40, 0, 41), cmap="viridis")
        for t, p in (MULTIPORT_BEAMS if port.port_id == 4 else [MULTIPORT_BEAMS[port.port_id]]):
            ax.plot([np.sin(np.radians(t)) * np.cos(np.radians(p))], [np.sin(np.radians(t)) * np.sin(np.radians(p))], "r+", ms=10)
        ax.set_aspect("equal")
        ax.set_title(f"port {port.port_id + 1}: D = {pat.directivity_estimate_dbi:.1f} dBi, "
                     f"eff {eff:.2f}", fontsize=9)
    axes.ravel()[-1].axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "multiport_patterns.png", dpi=150)
    plt.close(fig)

    lw = board.lw_grid_m
    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    for i, row in enumerate(lw):
        for j, cell in enumerate(row):
            if cell is None or cell[0] <= 0 or cell[1] <= 0:
                continue
            cx, cy = board.cell_center_m(i, j)
            ax.add_patch(plt.Rectangle((cx * 1e3 - cell[0] * 5e2, cy * 1e3 - cell[1] * 5e2), cell[0] * 1e3, cell[1] * 1e3,
                                       fc="#b87333", ec="none"))
    colors = ["tab:blue", "tab:green", "tab:purple", "tab:orange", "tab:red"]
    flagged = {f.port_a for f in multi.isolation_flags if f.flagged} | {f.port_b for f in multi.isolation_flags if f.flagged}
    for port in multi.ports:
        fx, fy = port.feed_xy_m[0] * 1e3, port.feed_xy_m[1] * 1e3
        ax.plot([fx], [fy], marker="o", ms=9, color=colors[port.port_id], mec="k")
        ax.annotate(f"port {port.port_id + 1}", (fx, fy), xytext=(6, 6), textcoords="offset points", fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    n_flagged = sum(1 for f in multi.isolation_flags if f.flagged)
    ax.set_title(f"Multi-port board, 8 mm cells, feed markers\n"
                 f"isolation flags: {n_flagged} of {len(multi.isolation_flags)} port pairs (edge gap under 1 lambda_sw)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "multiport_patch_layout.png", dpi=150)
    plt.close(fig)

    isolation = [
        {"ports": [f.port_a + 1, f.port_b + 1], "gap_lambda_sw": round(f.gap_lambda_sw, 3), "flagged": f.flagged}
        for f in multi.isolation_flags
    ]
    metrics["multiport"] = {"ports": per_port_metrics, "isolation": isolation}

    substrate, metal = SUBSTRATE, Metal()
    feeds_mm = [(p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports]
    ports_info = [(p.port_id + 1, p.feed_xy_m[0] * 1e3, p.feed_xy_m[1] * 1e3) for p in multi.ports]
    step_dir = OUT / "multiport_step"
    export_multi_port_array(board.lw_grid_m, PERIOD_M, PERIOD_M, substrate, metal, step_dir, include_step=True,
                            feeds_mm=feeds_mm, ports_info=ports_info)

    inputs = {
        "field": "analytic TM0 monopole, phase-only, per port (no per-port .fld in this app version)",
        "k_swg": f"tm0_beta = {k_swg:.2f} rad/m",
        "layout": "2x2 small ports (60 mm) plus one large port, gap 0 mm",
        "beams": "ports 1-4: one beam each (15,0), (30,90), (45,180), (60,270); port 5: all four beams",
        "polarization": "RHCP compensation on every beam (pencil beams)",
        "cells": f"{PERIOD_M * 1e3:g} mm pitch; formulation B with the square-patch library",
        "scripts": "examples/make_submission_figures.py",
    }
    text = ("Multi-port board: four single-beam sub-apertures and one quad-beam sub-aperture, each with its own "
            "feed. Per-port patterns use the vector model; gain is reported as directivity and as aperture "
            "efficiency against the port's own uniform-aperture maximum. The isolation flag is a first-order "
            "edge-distance check in lambda_sw, not a full-wave result.")
    numbers = {}
    for pm in per_port_metrics:
        numbers[f"port {pm['port']} requested beams (theta, phi deg)"] = pm["beams"]
        numbers[f"port {pm['port']} main beam realised (theta, phi deg)"] = (pm["main_beam_theta_deg"], pm["main_beam_phi_deg"])
        numbers[f"port {pm['port']} directivity (dBi) / aperture efficiency"] = (
            f"{pm['directivity_dbi']:.2f} / {pm['gain_efficiency_fraction']:.3f} "
            f"(port maximum {pm['max_directivity_dbi']:.2f} dBi, aperture {pm['aperture_m2'] * 1e4:.0f} cm^2)"
        )
    numbers["isolation flags (first-order, edge gap in lambda_sw)"] = (
        f"{sum(1 for x in isolation if x['flagged'])} of {len(isolation)} port pairs flagged"
    )
    _write_caption("multiport", text, inputs, numbers)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    work = ROOT / "results" / "submission_work"
    work.mkdir(parents=True, exist_ok=True)
    candidates, l_values = _library_candidates(work)
    reference, feed_xy_m = _reference()
    green = LayeredGreenFunction(SUBSTRATE, F_HZ, backing="pec")
    metrics = {"library_L_mm": list(l_values), "candidates": len(candidates)}
    _single_beam_figures(reference, feed_xy_m, candidates, green, metrics)
    _multiport_figures(candidates, green, metrics, work)
    (OUT / "submission_metrics.json").write_text(json.dumps(metrics, indent=2, default=float) + "\n")
    print(json.dumps({k: metrics[k] for k in ("vortex_natural", "pencil_rhcp", "vortex_l1")}, indent=2, default=float))


if __name__ == "__main__":
    main()
