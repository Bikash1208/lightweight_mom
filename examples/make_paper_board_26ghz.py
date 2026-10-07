"""Mishra et al. (IEEE TVT 2026) 1-bit coded 5-port board at 26 GHz, run through the app's own coded
ports path, with both beta modes, the paper's beams (LHCP) and a comparison against the paper's
reported values.

Run from the project root (needs the `cad` extra for the STEP file):
    .venv/bin/python examples/make_paper_board_26ghz.py

Cell pair (read from the paper's Fig. 2, not its text): period 2.57 mm; state 0 is EMPTY bare
RT/Duroid 5880 (1.575 mm, eps_r 2.2) with no metal; state 1 is a 2.17 mm square patch with an
annular slot (r 0.46 mm, R 0.75 mm). Both states are solved here with the app's periodic MoM
(state 1) and the bare-stack reflection (state 0, no currents). Their reflections are referred to
the patch plane, then passed to the coded ports path as uploaded values.

Writes into docs/submission/paper26/: PNG figures, `comparison.md`, `metrics.json`, and
`step/` (STEP + HFSS_README.md for the paper-mode board).
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from src.app.routers.hmsa import (  # noqa: E402
    HmsaSynthesizePortsRequest,
    TwoStateSourceRequest,
    _annular_library_base_config,
    _run_synthesize_ports_job,
)
from src.dispersion.sheet_model import find_beta_from_gamma  # noqa: E402
from src.excitation import ReferencePlanes  # noqa: E402
from src.export.array_export import export_multi_port_array  # noqa: E402
from src.green_function import LayeredGreenFunction  # noqa: E402
from src.hologram.two_state_cells import LIBRARY_REFERENCE_OFFSET_M, patch_plane_gamma  # noqa: E402
from src.library_sweep import build_annular_slot_library_cells  # noqa: E402
from src.materials import Metal, Substrate  # noqa: E402
from src.pattern.vector_aperture_pattern import (  # noqa: E402
    aperture_current_from_scalar,
    aperture_efficiency_estimate,
    circular_components,
    compute_vector_pattern,
    copol_label,
)
from src.sparameters import scattering_from_fields  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "submission" / "paper26"
C0 = 299792458.0
F_HZ = 26e9
PERIOD_MM = 2.57
EPS_R, TAN_D, THICK_MM = 2.2, 0.0009, 1.575
PATCH_MM, R_MM, R_SLOT_MM = 2.17, 0.75, 0.46
LAMBDA0_MM = C0 / F_HZ * 1e3
SMALL_PORT_MM = 5.0 * LAMBDA0_MM
SUBSTRATE = Substrate("user", EPS_R, TAN_D, THICK_MM * 1e-3)
METAL = Metal()
PERIOD_M = PERIOD_MM * 1e-3
BEAMS = [(15.0, 0.0), (30.0, 90.0), (45.0, 180.0), (60.0, 270.0)]
PAPER_LHCP_GAIN_DBIC = {"port 1": 16.8, "port 2": 16.5, "port 3": 15.5, "port 4": 14.8, "port 5": 15.8}
PAPER_RHCP_SUPPRESSION_DB = "10-15"
PAPER_ISOLATION_DB = "< -30"
PAPER_S11_DB = "< -10 over 25-27 GHz"
EQUALISE_ITERATIONS = 2
BEAM_RADIUS_DEG = 12.0
EDGE_RHO_MM = SMALL_PORT_MM / 2.0


def _gamma_top_empty(frequency_hz):
    green = LayeredGreenFunction(SUBSTRATE, frequency_hz, backing="pec")
    planes = ReferencePlanes(top_air_thickness=5.0e-3, bottom_air_thickness=None, substrate_thickness=THICK_MM * 1e-3)
    return complex(scattering_from_fields(green, planes, np.zeros((1, 2, 2), dtype=complex)).s11[0, 0])


def _gamma_top_slotted(frequency_hz):
    body = {
        "period_mm": PERIOD_MM,
        "eps_r": EPS_R,
        "tan_delta": TAN_D,
        "thickness_mm": THICK_MM,
        "frequency_hz": frequency_hz,
        "accuracy": "fast",
        "a_mm": PATCH_MM,
        "b_mm": PATCH_MM,
        "r_mm": R_SLOT_MM,
    }
    base = _annular_library_base_config(body, placeholder_R_mm=R_MM)
    cell = build_annular_slot_library_cells(base, PATCH_MM, PATCH_MM, R_SLOT_MM, [R_MM], frequency_hz)[0]
    if cell["solver_status"] != "OK":
        raise RuntimeError(f"annular cell solve failed at {frequency_hz / 1e9} GHz: {cell['solver_status']}")
    entry = cell["S11"][0]
    return complex(entry["real"], entry["imag"])


def cell_pair_table():
    rows = []
    for f in (25e9, 26e9, 27e9):
        g_empty = patch_plane_gamma(_gamma_top_empty(f), f, LIBRARY_REFERENCE_OFFSET_M)
        g_slot = patch_plane_gamma(_gamma_top_slotted(f), f, LIBRARY_REFERENCE_OFFSET_M)
        beta = {}
        for name, g in (("empty", g_empty), ("slotted", g_slot)):
            mode = find_beta_from_gamma(g, SUBSTRATE, f, period_m=PERIOD_M)
            beta[name] = mode.beta_rad_per_m
        sep = abs((np.degrees(np.angle(g_slot)) - np.degrees(np.angle(g_empty)) + 180.0) % 360.0 - 180.0)
        rows.append(
            {
                "frequency_GHz": f / 1e9,
                "empty_patch_plane_phase_deg": float(np.degrees(np.angle(g_empty))),
                "slotted_patch_plane_phase_deg": float(np.degrees(np.angle(g_slot))),
                "empty_magnitude_db": float(20 * np.log10(abs(g_empty))),
                "slotted_magnitude_db": float(20 * np.log10(abs(g_slot))),
                "empty_beta_rad_per_m": float(beta["empty"]),
                "slotted_beta_rad_per_m": float(beta["slotted"]),
                "separation_deg": float(sep),
                "magnitude_difference_db": float(20 * np.log10(abs(g_slot)) - 20 * np.log10(abs(g_empty))),
                "_gamma_empty": g_empty,
                "_gamma_slotted": g_slot,
            }
        )
    return rows


def _two_state_request(gamma_empty, gamma_slot):
    return TwoStateSourceRequest(
        source="uploaded",
        uploaded_label_0="empty (bare substrate)",
        uploaded_value_0_real=gamma_empty.real,
        uploaded_value_0_imag=gamma_empty.imag,
        uploaded_label_1="2.17 mm square, annular slot r 0.46 / R 0.75 mm",
        uploaded_value_1_real=gamma_slot.real,
        uploaded_value_1_imag=gamma_slot.imag,
        uploaded_geometry_0={"kind": "empty"},
        uploaded_geometry_1={
            "kind": "annular_slot",
            "a_mm": PATCH_MM,
            "b_mm": PATCH_MM,
            "R_mm": R_MM,
            "r_mm": R_SLOT_MM,
        },
    ).model_dump()


def _port_beam_dict(theta, phi, weight):
    return {
        "theta_deg": theta,
        "phi_deg": phi,
        "weight_real": weight.real,
        "weight_imag": weight.imag,
        "initial_phase_deg": 0.0,
        "desired_polarization": "lhcp",
        "oam_charge": 0,
    }


def _board_body(beta_mode, gamma_empty, gamma_slot, port5_weights):
    ports = [{"beams": [_port_beam_dict(*BEAMS[i], 1.0 + 0j)]} for i in range(4)]
    ports.append({"beams": [_port_beam_dict(t, p, w) for (t, p), w in zip(BEAMS, port5_weights, strict=True)]})
    return HmsaSynthesizePortsRequest(
        formulation="coded",
        frequency_hz=F_HZ,
        period_mm=PERIOD_MM,
        layout={"preset": "2x2_plus_1_large", "small_side_mm": SMALL_PORT_MM, "gap_mm": 0.0},
        ports=ports,
        port_eps_r=EPS_R,
        port_thickness_mm=THICK_MM,
        port_tan_delta=TAN_D,
        two_state=_two_state_request(gamma_empty, gamma_slot),
        beta_mode=beta_mode,
    ).model_dump()


class _Job:
    def __init__(self, job_id):
        self.id = job_id
        self.progress = 0.0
        self.progress_message = ""


def run_board(beta_mode, gamma_empty, gamma_slot):
    weights = [1.0 + 0j] * 4
    history = []
    for iteration in range(EQUALISE_ITERATIONS + 1):
        result = _run_synthesize_ports_job(
            _Job(f"paper26_{beta_mode}_{iteration}"), _board_body(beta_mode, gamma_empty, gamma_slot, weights)
        )
        port5 = result["ports"][4]
        pat, _ = _port_pattern(port5, F_HZ)
        levels = [_beam_level(pat, t, p)["level_db"] for t, p in BEAMS]
        history.append([round(v, 2) for v in levels])
        if iteration == EQUALISE_ITERATIONS:
            break
        weights = [w * 10.0 ** (-level / 20.0) for w, level in zip(weights, levels, strict=True)]
    return result, history


def _port_pattern(port, frequency_hz):
    x_m, y_m = np.asarray(port["x_m"]), np.asarray(port["y_m"])
    field = np.asarray(port["aperture_field_real"]) + 1j * np.asarray(port["aperture_field_imag"])
    feed = (port["feed_x_mm"] * 1e-3, port["feed_y_mm"] * 1e-3)
    green = LayeredGreenFunction(SUBSTRATE, frequency_hz, backing="pec")
    jx, jy = aperture_current_from_scalar(x_m, y_m, field, feed_xy_m=feed)
    dx = float(x_m[0, 1] - x_m[0, 0])
    dy = float(y_m[1, 0] - y_m[0, 0])
    return compute_vector_pattern(x_m, y_m, jx, jy, frequency_hz, green, dx, dy), (jx, jy)


def _direction_cosines(theta_deg, phi_deg):
    t, p = np.radians(theta_deg), np.radians(phi_deg)
    return np.stack([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)], axis=-1)


def _beam_level(pat, theta_deg, phi_deg):
    """Local peak within BEAM_RADIUS_DEG of the requested direction: LHCP (co) and RHCP (cross)
    levels there, and the pointing error. The copol component follows `copol_label` (the
    project's own handedness rule), not the basis letter."""
    req = _direction_cosines(np.array(theta_deg), np.array(phi_deg))
    grid = _direction_cosines(pat.theta_deg, pat.phi_deg)
    near = (grid @ req) >= np.cos(np.radians(BEAM_RADIUS_DEG))
    power = np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2
    peak_power = float(np.max(power))
    masked = np.where(near, power, -1.0)
    i, j = np.unravel_index(int(np.argmax(masked)), masked.shape)
    co_key = copol_label("lhcp").split("_")[1]
    co = pat.e_r[i, j] if co_key == "R" else pat.e_l[i, j]
    cross = pat.e_l[i, j] if co_key == "R" else pat.e_r[i, j]
    realised = _direction_cosines(pat.theta_deg[i, j], pat.phi_deg[i, j])
    pointing_error = float(np.degrees(np.arccos(np.clip(realised @ req, -1.0, 1.0))))
    return {
        "theta_deg": float(pat.theta_deg[i, j]),
        "phi_deg": float(pat.phi_deg[i, j]),
        "pointing_error_deg": pointing_error,
        "level_db": float(10.0 * np.log10(max(power[i, j], 1e-30) / peak_power)),
        "lhcp_db": float(20.0 * np.log10(max(abs(co), 1e-12) / np.sqrt(peak_power))),
        "rhcp_suppression_db": float(20.0 * np.log10(max(abs(co), 1e-12) / max(abs(cross), 1e-12))),
    }


def _port_metrics(port, pat):
    beams = [(t, p) for t, p in BEAMS] if port["port_id"] == 4 else [BEAMS[port["port_id"]]]
    area = (port["x_max_mm"] - port["x_min_mm"]) * (port["y_max_mm"] - port["y_min_mm"]) * 1e-6
    eff, _ = aperture_efficiency_estimate(pat.directivity_estimate_dbi, area, F_HZ)
    per_beam = []
    for t, p in beams:
        b = _beam_level(pat, t, p)
        b["requested_theta_deg"], b["requested_phi_deg"] = t, p
        b["directivity_dbi"] = float(pat.directivity_estimate_dbi + b["level_db"])
        per_beam.append(b)
    return {
        "port": port["port_id"] + 1,
        "feed_mm": [round(port["feed_x_mm"], 2), round(port["feed_y_mm"], 2)],
        "cells": int(port["n_cells"]),
        "unreachable_percent": round(100.0 * port["unreachable_fraction"], 2),
        "peak_directivity_dbi": float(pat.directivity_estimate_dbi),
        "aperture_efficiency_upper_bound": float(eff),
        "beams": per_beam,
    }


def _uv_panel(ax, pat, beams, title):
    u = np.sin(np.radians(pat.theta_deg)) * np.cos(np.radians(pat.phi_deg))
    v = np.sin(np.radians(pat.theta_deg)) * np.sin(np.radians(pat.phi_deg))
    ax.contourf(
        u,
        v,
        np.maximum(
            10
            * np.log10(
                np.maximum(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2, 1e-30)
                / np.max(np.abs(pat.e_theta) ** 2 + np.abs(pat.e_phi) ** 2)
            ),
            -40.0,
        ),
        levels=np.linspace(-40, 0, 41),
        cmap="viridis",
    )
    for t, p in beams:
        ax.plot(
            [np.sin(np.radians(t)) * np.cos(np.radians(p))],
            [np.sin(np.radians(t)) * np.sin(np.radians(p))],
            "r+",
            ms=10,
        )
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8)


def figure_patterns(port_patterns, path):
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, (port, pat) in zip(axes.ravel(), port_patterns, strict=False):
        beams = BEAMS if port["port_id"] == 4 else [BEAMS[port["port_id"]]]
        _uv_panel(
            ax, pat, beams, f"port {port['port_id'] + 1}: LHCP power, peak {pat.directivity_estimate_dbi:.1f} dBi"
        )
    axes.ravel()[-1].axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def figure_cuts(port_patterns, path):
    fig, axes = plt.subplots(len(BEAMS), 1, figsize=(7.0, 2.2 * len(BEAMS)), sharex=True)
    port5_pat = next(pat for port, pat in port_patterns if port["port_id"] == 4)
    theta = port5_pat.theta_deg[:, 0]
    for ax, (t, p) in zip(axes, BEAMS, strict=True):
        i_cut = int(np.argmin(np.abs(((port5_pat.phi_deg[0, :] - p + 180.0) % 360.0) - 180.0)))
        e_theta = port5_pat.e_theta[:, i_cut]
        e_phi = port5_pat.e_phi[:, i_cut]
        peak = float(np.max(np.abs(port5_pat.e_theta) ** 2 + np.abs(port5_pat.e_phi) ** 2))
        e_r, e_l = circular_components(e_theta, e_phi)
        co_key = copol_label("lhcp").split("_")[1]
        co, cross = (e_r, e_l) if co_key == "R" else (e_l, e_r)
        ax.plot(theta, 10 * np.log10(np.maximum(np.abs(co) ** 2 / peak, 1e-12)), lw=2.0, label="LHCP (co)")
        ax.plot(
            theta, 10 * np.log10(np.maximum(np.abs(cross) ** 2 / peak, 1e-12)), lw=1.2, ls="--", label="RHCP (cross)"
        )
        ax.axvline(t, color="r", ls=":", lw=1)
        ax.set_ylim(-45, 2)
        ax.set_ylabel(f"phi={p:g}\n(dB)", fontsize=8)
    axes[0].legend(fontsize=8)
    axes[-1].set_xlabel("theta (deg), port 5 (all four beams)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _board_feeds_mm(result):
    """Feed positions in the board-corner frame the combined cell grid uses (the ports report
    global coordinates; the layout preset puts the board's lower-left corner at negative x and y)."""
    x0 = min(p["x_min_mm"] for p in result["ports"])
    y0 = min(p["y_min_mm"] for p in result["ports"])
    return [(p["feed_x_mm"] - x0, p["feed_y_mm"] - y0) for p in result["ports"]]


def figure_layout(result, path):
    board = result["lw_grid_m"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    n_y, n_x = len(board), len(board[0])
    for i in range(n_y):
        for j in range(n_x):
            cell = board[i][j]
            cx = (j + 0.5) * PERIOD_MM
            cy = (i + 0.5) * PERIOD_MM
            if isinstance(cell, dict):
                ax.add_patch(
                    plt.Rectangle(
                        (cx - PERIOD_MM / 2, cy - PERIOD_MM / 2), PERIOD_MM, PERIOD_MM, fc="#f2f2f2", ec="none"
                    )
                )
                ax.add_patch(
                    plt.Rectangle((cx - PATCH_MM / 2, cy - PATCH_MM / 2), PATCH_MM, PATCH_MM, fc="#b87333", ec="none")
                )
                ax.add_patch(plt.Circle((cx, cy), R_MM, fc="#f2f2f2", ec="none"))
                ax.add_patch(plt.Circle((cx, cy), R_SLOT_MM, fc="#b87333", ec="none"))
            elif cell[0] > 0 and cell[1] > 0:
                ax.add_patch(
                    plt.Rectangle(
                        (cx - cell[0] * 500, cy - cell[1] * 500), cell[0] * 1e3, cell[1] * 1e3, fc="#b87333", ec="none"
                    )
                )
    colours = ["tab:blue", "tab:green", "tab:purple", "tab:orange", "tab:red"]
    for (fx, fy), port in zip(_board_feeds_mm(result), result["ports"], strict=True):
        ax.plot([fx], [fy], "o", ms=8, color=colours[port["port_id"]], mec="k")
    ax.set_xlim(0, n_x * PERIOD_MM)
    ax.set_ylim(0, n_y * PERIOD_MM)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(
        f"{F_HZ / 1e9:g} GHz board, {n_x} x {n_y} cells at {PERIOD_MM} mm (copper = annular-slot patch, grey = empty)",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_comparison(metrics, path):
    lines = [
        "# 26 GHz 1-bit coded 5-port board: this model vs Mishra et al. (TVT 2026)",
        "",
        "Cell pair read from the paper's Fig. 2 (cell 1 = bare substrate, cell 2 = square patch with annular slot), not from its text.",
        "",
        "## Cell pair, 2.57 mm period",
        "",
        "| f (GHz) | empty phase (deg) | slotted phase (deg) | separation (deg) | empty |S11| (dB) | slotted |S11| (dB) | |S11| difference (dB) | beta empty (rad/m) | beta slotted (rad/m) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for c in metrics["cells"]:
        lines.append(
            f"| {c['frequency_GHz']:g} | {c['empty_patch_plane_phase_deg']:.2f} | {c['slotted_patch_plane_phase_deg']:.2f} | "
            f"{c['separation_deg']:.2f} | {c['empty_magnitude_db']:.3f} | {c['slotted_magnitude_db']:.3f} | "
            f"{c['magnitude_difference_db']:.3f} | {c['empty_beta_rad_per_m']:.2f} | {c['slotted_beta_rad_per_m']:.2f} |"
        )
    b = metrics["beta_modes"]
    lines += [
        "",
        "## Surface-wave beta modes at 26 GHz",
        "",
        f"- Paper mode, beta = k0 sqrt((eps_r+1)/2): {b['paper_beta_rad_per_m']:.2f} rad/m.",
        f"- Model mode, beta of the mixed 1-bit surface (equal-area mean of the two sheet admittances): {b['model_beta_rad_per_m']:.2f} rad/m.",
        f"- Aperture-edge phase difference between the modes at the 28.8 mm edge of a 5 lambda_0 sub-aperture: "
        f"{b['edge_phase_difference_deg_at_5lambda0_edge']:.1f} deg.",
        "",
    ]
    lines += ["## Board", "", "| | Paper | This model, paper mode | This model, model mode |", "|---|---|---|---|"]
    bd = metrics["board"]
    lines.append(
        f"| Board | 20 lambda_0 x 10 lambda_0 | {bd['paper']['width_mm']:.1f} x {bd['paper']['height_mm']:.1f} mm | "
        f"{bd['model']['width_mm']:.1f} x {bd['model']['height_mm']:.1f} mm |"
    )
    lines.append(
        f"| Cells | 92 x 46 | {bd['paper']['cells_x']} x {bd['paper']['cells_y']} | {bd['model']['cells_x']} x {bd['model']['cells_y']} |"
    )
    lines += [
        "",
        "## Beams: pointing and LHCP/RHCP (per requested beam)",
        "",
        "The paper gives the target directions, not measured beam angles in its text (measured patterns are in its Fig. 5 only). "
        "Levels are relative to each port's own peak; the paper's gains are realised dBic, ours are lossless directivity (dBi).",
        "",
        "| port | requested (theta, phi) | realised (theta, phi) | pointing error (deg) | LHCP level (dB rel. port peak) | RHCP suppression (dB) | directivity (dBi, lossless) | paper realised LHCP gain (dBic) | mode |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for mode in ("paper", "model"):
        for port in metrics["ports"][mode]:
            paper_gain = PAPER_LHCP_GAIN_DBIC[f"port {port['port']}"]
            for beam in port["beams"]:
                lines.append(
                    f"| {port['port']} | ({beam['requested_theta_deg']:g}, {beam['requested_phi_deg']:g}) | "
                    f"({beam['theta_deg']:.1f}, {beam['phi_deg']:.1f}) | {beam['pointing_error_deg']:.1f} | "
                    f"{beam['level_db']:.2f} | {beam['rhcp_suppression_db']:.1f} | {beam['directivity_dbi']:.2f} | "
                    f"{paper_gain} | {mode} |"
                )
    lines += [
        "",
        "RHCP suppression is shown as the RHCP level below LHCP at the beam, so a positive number means RHCP is that many dB down. "
        f"Paper: {PAPER_RHCP_SUPPRESSION_DB} dB.",
        "",
    ]
    lines += [
        "## Other paper targets",
        "",
        "| Quantity | Paper | This model |",
        "|---|---|---|",
        f"| Isolation between ports | {PAPER_ISOLATION_DB} dB (measured) | not computed: first-order edge-gap flags only; check in HFSS with the STEP |",
        f"| |S11| | {PAPER_S11_DB} | not comparable: a reflectarray cell reflects almost all power (|S11| about 0 dB by construction); the paper's |S11| is a feed-match quantity |",
        "",
    ]
    lines += ["## Port 5 equalisation", "", "Beam levels (dB relative to port 5's peak) per iteration:", ""]
    for mode in ("paper", "model"):
        lines.append(f"- {mode} mode: " + "; ".join(str(h) for h in metrics["port5_equalisation_db"][mode]))
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells = cell_pair_table()
    gamma_empty = cells[1]["_gamma_empty"]
    gamma_slot = cells[1]["_gamma_slotted"]
    metrics = {"cells": [{k: v for k, v in row.items() if not k.startswith("_")} for row in cells]}

    boards = {}
    for beta_mode in ("paper", "model"):
        result, history = run_board(beta_mode, gamma_empty, gamma_slot)
        patterns = [(port, _port_pattern(port, F_HZ)[0]) for port in result["ports"]]
        boards[beta_mode] = {
            "result": result,
            "history": history,
            "patterns": patterns,
            "ports": [_port_metrics(port, pat) for port, pat in patterns],
        }
    paper_result = boards["paper"]["result"]
    figure_patterns(boards["paper"]["patterns"], OUT / "patterns_paper_mode.png")
    figure_cuts(boards["paper"]["patterns"], OUT / "port5_circular_cuts_paper_mode.png")
    figure_layout(paper_result, OUT / "board_layout.png")

    edge = {
        "paper_beta_rad_per_m": paper_result["two_state"]["beta_paper_rad_per_m"],
        "model_beta_rad_per_m": paper_result["two_state"]["beta_model_rad_per_m"],
    }
    edge["beta_difference_rad_per_m"] = edge["paper_beta_rad_per_m"] - edge["model_beta_rad_per_m"]
    edge["edge_phase_difference_deg_at_5lambda0_edge"] = float(
        np.degrees(edge["beta_difference_rad_per_m"] * EDGE_RHO_MM * 1e-3)
    )
    metrics["beta_modes"] = edge
    metrics["separation_at_26ghz_deg"] = cells[1]["separation_deg"]
    metrics["two_state"] = {m: boards[m]["result"]["two_state"] for m in boards}
    metrics["ports"] = {m: boards[m]["ports"] for m in boards}
    metrics["port5_equalisation_db"] = {m: boards[m]["history"] for m in boards}
    metrics["board"] = {
        m: {
            "width_mm": boards[m]["result"]["board_width_mm"],
            "height_mm": boards[m]["result"]["board_height_mm"],
            "cells_x": len(boards[m]["result"]["lw_grid_m"][0]),
            "cells_y": len(boards[m]["result"]["lw_grid_m"]),
        }
        for m in boards
    }

    step_dir = OUT / "step"
    board_grid = paper_result["lw_grid_m"]
    feeds_mm = paper_result["export_feeds_mm"]
    ports_info = [(p["port_id"] + 1, fx, fy) for p, (fx, fy) in zip(paper_result["ports"], feeds_mm, strict=True)]
    ox, oy = paper_result["export_offset_mm"]
    marker_check = [
        {
            "port": p["port_id"] + 1,
            "feed_mm": [round(fx, 3), round(fy, 3)],
            "inside_own_port": bool(
                p["x_min_mm"] - ox <= fx <= p["x_max_mm"] - ox and p["y_min_mm"] - oy <= fy <= p["y_max_mm"] - oy
            ),
        }
        for p, (fx, fy) in zip(paper_result["ports"], feeds_mm, strict=True)
    ]
    metrics["feed_markers_in_step_frame"] = marker_check
    if not all(m["inside_own_port"] for m in marker_check):
        raise RuntimeError(f"a feed marker lies outside its own port in the STEP frame: {marker_check}")
    export = export_multi_port_array(
        board_grid,
        PERIOD_M,
        PERIOD_M,
        SUBSTRATE,
        METAL,
        step_dir,
        include_step=True,
        feeds_mm=feeds_mm,
        ports_info=ports_info,
    )
    metrics["step"] = {
        "step_path": str(Path(export.step_path).relative_to(ROOT)) if export.step_path else None,
        "file_size_bytes": export.file_size_bytes,
    }

    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2, default=float) + "\n")
    write_comparison(metrics, OUT / "comparison.md")
    print(
        json.dumps({"cells": metrics["cells"], "beta_modes": edge, "board": metrics["board"]}, indent=2, default=float)
    )


if __name__ == "__main__":
    main()
