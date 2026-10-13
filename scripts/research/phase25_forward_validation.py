"""Phase 25 A: validation of the arbitrary-excitation forward model (`src/hologram/butler_forward.py`).

    python scripts/research/phase25_forward_validation.py <out.json> [--quick] [--skip-mom] [--mom-only[=24,20]]

Part I (pattern path, hex4 side 150 mm, Phase 16 pointing-corrected, 10.54 GHz, held-out 0.5 deg x 1 deg grid):
  for the DESIGNED excitations (each monopole alone; all four together) the forward model must reproduce the
  existing synthesis prediction (`sll_board.evaluate`, which evaluates the stored `aperture_field`):
    (a) with the design's own incident field (`incident='stored'`): complex far field to rounding and every
        metric identical;
    (b) with the incident field rebuilt analytically at the cell centres (`incident='analytic'`): the
        difference this reconstruction makes (the synthesis interpolates a finer native grid);
    (c) `amplitude_model='physical'` (rho^-1/2 taper): how far the metrics move if the surface-wave
        amplitude taper is kept (the synthesis is phase-only).

Part II (what the repo's full-array MoM CAN check): the only multi-cell MoM here is the periodic N x N
supercell under a NORMAL-INCIDENCE PLANE WAVE (`phase12c.grounded_supercell`, `phase12d.floquet_power`) -- there
is no finite board fed by a surface wave, so the monopole's incident field E_inc is NOT testable by it and is
reported as NOT_AVAILABLE_FROM_CURRENT_DATA. What it does test is the other half of the forward model, the
local-periodicity step `cell radiates Gamma_state x E_inc`: for a 1-bit 6.4 mm-cell supercell (empty / 5.85 mm
square, grounded RT/duroid 5880 1.575 mm, 10.54 GHz) the forward model predicts every Floquet order's
reflected-E amplitude as the Fourier sum of the ISOLATED per-cell reflection phasors (same mesh, patch
plane), and the full N x N MoM solve gives the true value including mutual coupling. Cases: a designed block
(cropped from the hex4 design), a block from another port's pattern (non-designed for this block), a seeded
random block, and the uniform board (known answer: must agree to rounding).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "research"))

import phase23_compare as v1  # noqa: E402

from src.app.routers.hmsa import HmsaSynthesizePortsRequest, _run_synthesize_ports_job  # noqa: E402
from src.hologram.butler_forward import (  # noqa: E402
    FeedFarField,
    board_from_cells,
    board_from_ports_result,
    evaluate_excitation,
)
from src.hologram.sll_board import HeldOutGrid, _FarField, build_problem, evaluate  # noqa: E402

F0, PERIOD, EPS_R, TAND, H_MM = v1.F0, v1.PERIOD, v1.EPS_R, v1.TAND, v1.H_MM
BEAMS = v1.BEAMS


def _synthesize():
    body = {
        "formulation": "coded",
        "frequency_hz": F0,
        "period_mm": PERIOD,
        "layout": {"preset": "hex4", "side_mm": 150.0, "gap_mm": 0.0},
        "ports": [{"beams": [v1._beam(t, p)]} for t, p in BEAMS],
        "port_eps_r": EPS_R,
        "port_thickness_mm": H_MM,
        "port_tan_delta": TAND,
        "two_state": {"source": "auto"},
        "beta_mode": "model",
        "pointing_correction": True,
        "pointing_iterations": 2,
    }
    return _run_synthesize_ports_job(v1._Job("p25_forward_validation"), HmsaSynthesizePortsRequest(**body).model_dump())


_KEYS = (
    "directivity_dbi",
    "pointing_error_deg",
    "theta_deg",
    "phi_deg",
    "sll_global_db",
    "sll_xcut_db",
    "sll_ycut_db",
    "xpd_db",
)


def _compare(ev_ref, ev_new):
    """Largest absolute difference of every beam metric (dB / deg) between two `_scenario_eval` dicts."""
    out = {}
    for k in _KEYS:
        d = []
        for br, bn in zip(ev_ref["beams"], ev_new["beams"], strict=True):
            if br.get(k) is None or bn.get(k) is None:
                continue
            d.append(abs(float(br[k]) - float(bn[k])))
        out[k] = max(d) if d else None
    out["peak_directivity_dbi"] = abs(ev_ref["peak_directivity_dbi"] - ev_new["peak_directivity_dbi"])
    return out


def _row(ev):
    return {
        "directivity_dbi": [round(b["directivity_dbi"], 3) for b in ev["beams"]],
        "pointing_error_deg": [round(b["pointing_error_deg"], 3) for b in ev["beams"]],
        "sll_global_db": [None if b["sll_global_db"] is None else round(b["sll_global_db"], 2) for b in ev["beams"]],
        "xpd_db": [round(b["xpd_db"], 1) for b in ev["beams"]],
    }


def part_one(grid, quick):
    t0 = time.time()
    res = _synthesize()
    prob = build_problem(res, F0, EPS_R, TAND, H_MM, PERIOD, [[b] for b in BEAMS])
    ref = evaluate(prob, grid=grid)
    t_ref = time.time() - t0
    beams = [(t, p, "rhcp") for t, p in BEAMS]
    boards = {
        "stored": board_from_ports_result(res, PERIOD, F0, EPS_R, TAND, H_MM, incident="stored"),
        "analytic": board_from_ports_result(res, PERIOD, F0, EPS_R, TAND, H_MM, incident="analytic"),
    }
    a_board = boards["analytic"]
    boards["physical"] = board_from_cells(
        a_board.x_m,
        a_board.y_m,
        a_board.gamma,
        a_board.feeds_xy_m,
        a_board.period_m,
        a_board.frequency_hz,
        a_board.substrate,
        a_board.k_swg_rad_per_m,
        a_board.polarization,
        "physical",
        footprints=np.abs(a_board.e_inc) > 0,
        label="physical",
    )
    own = a_board.meta["owner_port"]
    scen = {f"port{p + 1}": ([p], np.eye(4)[p]) for p in range(4)}
    scen["all"] = ([0, 1, 2, 3], np.ones(4))
    out = {
        "reference_seconds": round(t_ref, 1),
        "n_cells": int(a_board.n_cells),
        "k_swg_rad_per_m": a_board.k_swg_rad_per_m,
    }
    caches = {k: FeedFarField(b) for k, b in boards.items()}
    for kind, board in boards.items():
        res_kind = {}
        for name, (ids, a) in scen.items():
            sub = [beams[i] for i in ids]
            ev = evaluate_excitation(board, caches[kind], a, sub, cell_port=own, beam_port=ids, grid=grid, name=name)
            res_kind[name] = {"vs_existing_max_abs_diff": _compare(ref[name], ev), "metrics": _row(ev)}
        out[kind] = res_kind
    out["existing"] = {name: _row(ref[name]) for name in scen}
    # complex far field, designed excitations, on a seeded direction set (stored incident field)
    rng = np.random.default_rng(25)
    u = rng.uniform(0.0, 1.0, 4000)
    theta = np.degrees(np.arccos(u))  # uniform on the hemisphere
    phi = rng.uniform(0.0, 360.0, 4000)
    levels = {p.port_id: p.state for p in prob.opt_ports}
    alpha = {p.port_id: 0.0 for p in prob.opt_ports}
    worst = {}
    for name, (ids, a) in scen.items():
        et_r, ep_r = _FarField(prob, levels, alpha, F0, set(ids))(theta, phi)
        et_m, ep_m = caches["stored"].fields(a, theta, phi)
        scale = max(np.abs(et_r).max(), np.abs(ep_r).max())
        worst[name] = float(max(np.abs(et_m - et_r).max(), np.abs(ep_m - ep_r).max()) / scale)
    out["complex_far_field_max_rel_diff_stored_incident"] = worst
    # the incident-field reconstruction itself
    s = a_board.meta["stored_aperture_field"]
    mine = (np.ones(4) @ a_board.e_inc) * a_board.gamma
    ph = np.degrees(np.angle(mine / s))
    out["analytic_vs_stored_aperture_field"] = {
        "median_abs_phase_deg": float(np.median(np.abs(ph))),
        "p99_abs_phase_deg": float(np.percentile(np.abs(ph), 99)),
        "max_abs_phase_deg": float(np.abs(ph).max()),
        "cells_over_1deg": int((np.abs(ph) > 1.0).sum()),
        "mean_stored_modulus": float(np.abs(s).mean()),
        "rms_complex_difference": float(np.sqrt(np.mean(np.abs(mine - s) ** 2))),
    }
    return out


# ------------------------------------------------------------------------------------------ part II


def _mom_supercell(lw_mm_grid, mesh_per_cell, substrate, metal, frequency_hz, period_m):
    """Floquet reflected-E_x amplitudes of every order of the N x N grounded supercell (patch plane)."""
    from src import solver as solver_mod
    from src.basis_functions import MeshSettings
    from src.excitation import PlaneWaveExcitation
    from src.green_function import LayeredGreenFunction
    from src.phase12c.grounded_supercell import build_grounded_supercell
    from src.phase12d.aperture_test import aperture_floquet_order
    from src.phase12d.floquet_power import floquet_order_powers

    n = len(lw_mm_grid)
    grid_m = [[(a * 1e-3, b * 1e-3) for a, b in row] for row in lw_mm_grid]
    sc = build_grounded_supercell(n, period_m, period_m, substrate, metal, grid_m)
    total = n * mesh_per_cell
    order = aperture_floquet_order(total)
    exc = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(total, total, order)
    system = solver_mod.assemble_system(sc, exc, mesh)
    cur = solver_mod.solve(system.z, system.v)
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    pw = floquet_order_powers(system.basis, cur, green, sc.lattice, max_order=2)
    return pw, int(system.basis.n_unknowns)


def _isolated_gamma(l_mm, mesh_per_cell, substrate, metal, frequency_hz, period_m, floquet_order):
    """Reflected E_x of the (0,0) order of ONE periodic cell at the patch plane: the same mesh size as the
    supercell (`mesh_per_cell` cells per period) and the SAME harmonic cut-off per cell lattice
    (`floquet_order` = the supercell's order // N), so the uniform board is the known answer to rounding."""
    from src import solver as solver_mod
    from src.basis_functions import MeshSettings
    from src.excitation import PlaneWaveExcitation
    from src.geometry import Geometry, UnitCell
    from src.green_function import LayeredGreenFunction
    from src.periodic_bc import RectangularLattice
    from src.phase12d.floquet_power import floquet_order_powers

    lat = RectangularLattice(period_m, period_m)
    cell = UnitCell(
        lat, Geometry(kind="rectangle", length=l_mm * 1e-3, width=l_mm * 1e-3), substrate, metal, backing="pec"
    )
    exc = PlaneWaveExcitation(frequency=frequency_hz, polarization="x")
    mesh = MeshSettings(mesh_per_cell, mesh_per_cell, floquet_order)
    system = solver_mod.assemble_system(cell, exc, mesh)
    cur = solver_mod.solve(system.z, system.v)
    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    pw = floquet_order_powers(system.basis, cur, green, lat, max_order=2)
    return complex(pw.e_reflected_x[2, 2])


def _empty_gamma(substrate, frequency_hz):
    """Reflected E_x (0,0) of the bare grounded slab at the patch plane (zero current)."""
    from src.constants import ETA0
    from src.green_function import LayeredGreenFunction
    from src.sparameters import decompose_air_fields

    green = LayeredGreenFunction(substrate, frequency_hz, backing="pec")
    e_stack, h_stack = green.stack_fields(0.0, 0.0, np.array([1.0, 0.0], dtype=complex), 0.0, port="top")
    up, _ = decompose_air_fields(e_stack, h_stack)
    del ETA0
    return complex(up[0])


# The MoM's rooftop mesh represents a patch as a whole number of mesh cells (centre-sample rule), symmetric
# about the cell centre, i.e. an even count. Mesh 6 / 8 would fill the whole 6.4 mm cell (neighbouring patches
# merge), so only meshes whose even cell count represents the 5.85 mm square without merging are used:
# cells-per-patch-side by cells-per-period.
ALIGNED_PATCH_CELLS = {12: 10, 16: 14, 20: 18, 24: 22}
# (mesh cells per 6.4 mm cell, N x N supercell). N = 5 (D = 1.12 lambda) is the smallest supercell with
# propagating orders besides the specular one at 10.54 GHz (N = 4 radiates everything into (0,0) by energy
# conservation, so it cannot test the forward model).
MOM_CASES = ((12, 5), (16, 5), (20, 5), (24, 5))


def part_two(res, mesh_per_cell, n_cells, seed=25, period_mm=None, patch_cells=None, tag=""):
    from src.materials import Metal, Substrate

    substrate = Substrate("user", EPS_R, TAND, H_MM * 1e-3)
    metal = Metal(model="pec")
    period_mm = PERIOD if period_mm is None else period_mm
    p_m = period_mm * 1e-3
    side_design = float(res["two_state"]["states"][1]["geometry"]["L_mm"])
    h_mm = period_mm / mesh_per_cell
    side = (patch_cells or ALIGNED_PATCH_CELLS[mesh_per_cell]) * h_mm  # the side the mesh actually represents (aligned)
    tiny = 1e-6  # mm: "no patch" (a zero-size patch is rejected by `Geometry`; it covers no mesh cell)
    t0 = time.time()
    g_empty = _empty_gamma(substrate, F0)
    from src.phase12d.aperture_test import aperture_floquet_order

    f_matched = aperture_floquet_order(n_cells * mesh_per_cell) // n_cells
    g_patch = _isolated_gamma(side, mesh_per_cell, substrate, metal, F0, p_m, f_matched)
    g_patch_full = _isolated_gamma(side, mesh_per_cell, substrate, metal, F0, p_m, 2 * mesh_per_cell)
    t_iso = time.time() - t0

    # maps: designed block (port 1), block of another port's pattern (port 3), seeded random, uniform
    def block(port, i0, j0):
        lev = np.zeros(np.asarray(res["ports"][port]["x_m"]).shape, dtype=int)
        for c in res["ports"][port]["patches_mm"]["cells"]:
            lev[c["i"], c["j"]] = 1
        return lev[i0 : i0 + n_cells, j0 : j0 + n_cells]

    lev1 = block(0, 20, 20)
    lev3 = block(2, 20, 20)
    rng = np.random.default_rng(seed)
    cases = {
        "designed_block_port1": lev1,
        "other_port_block_port3": lev3,
        "random_seeded": rng.integers(0, 2, size=(n_cells, n_cells)),
        "uniform_patch": np.ones((n_cells, n_cells), dtype=int),
    }
    # supercell orders (m, n) of the N-cell block: kx = 2 pi m / (N p), ky likewise
    out = {
        "mesh_cells_per_6p4mm_cell": mesh_per_cell,
        "n_cells": n_cells,
        "design_patch_side_mm": side_design,
        "patch_side_mm_represented": side,
        "gap_mm_represented": period_mm - side,
        "period_mm": period_mm,
        "gap_mm_design": PERIOD - side_design,
        "isolated_gamma_empty": [abs(g_empty), float(np.degrees(np.angle(g_empty)))],
        "isolated_gamma_square": [abs(g_patch), float(np.degrees(np.angle(g_patch)))],
        "isolated_floquet_order_matched": f_matched,
        "isolated_gamma_square_full_truncation_2x_mesh": [abs(g_patch_full), float(np.degrees(np.angle(g_patch_full)))],
        "isolated_solve_seconds": round(t_iso, 1),
        "cases": {},
    }
    k0 = 2.0 * np.pi * F0 / 299792458.0
    idx = [(m, n) for m in range(-2, 3) for n in range(-2, 3)]
    for name, lev in cases.items():
        t1 = time.time()
        lw = [[(side, side) if lev[i, j] else (tiny, tiny) for j in range(n_cells)] for i in range(n_cells)]
        pw, n_unk = _mom_supercell(lw, mesh_per_cell, substrate, metal, F0, p_m)
        gamma = np.where(lev == 1, g_patch, g_empty)
        xs = -n_cells * p_m / 2 + (np.arange(n_cells) + 0.5) * p_m
        xx, yy = np.meshgrid(xs, xs, indexing="xy")
        rows = {}
        for m, n in idx:
            kx, ky = 2 * np.pi * m / (n_cells * p_m), 2 * np.pi * n / (n_cells * p_m)
            kt2 = kx * kx + ky * ky
            if kt2 >= k0 * k0:
                continue  # evanescent: carries no radiated power
            kz = np.sqrt(k0 * k0 - kt2)
            amp = np.mean(gamma * np.exp(1j * (kx * xx + ky * yy)))  # point samples (the pattern code's rule)
            mom = complex(pw.e_reflected_x[2 + m, 2 + n])
            # an E_x-only plane wave in order (m, n) carries (1 - (ky/k0)^2)/(kz/k0) times the power of
            # the same amplitude at normal incidence (S_z = |E_x|^2 (k^2 - ky^2) / (2 w mu kz))
            w = (1.0 - (ky / k0) ** 2) / (kz / k0)
            rows[f"{m},{n}"] = {
                "mom_power_fraction": float(pw.fraction[2 + m, 2 + n]),
                "forward_power_fraction_unnormalised": float(w * abs(amp) ** 2),
                "mom_ex_abs": abs(mom),
                "forward_ex_abs": abs(amp),
            }
        total_mom = sum(r["mom_power_fraction"] for r in rows.values())
        total_fwd = sum(r["forward_power_fraction_unnormalised"] for r in rows.values())
        for r in rows.values():
            # the pattern model renormalises to its own radiated power, so the comparable quantity is the
            # forward fraction rescaled to the radiated total the MoM reports
            r["forward_power_fraction_renormalised"] = r["forward_power_fraction_unnormalised"] / total_fwd * total_mom
            r["error_percentage_points"] = 100.0 * abs(
                r["forward_power_fraction_renormalised"] - r["mom_power_fraction"]
            )
        out["cases"][name] = {
            "states": lev.tolist(),
            "n_unknowns": n_unk,
            "seconds": round(time.time() - t1, 1),
            "propagating_orders": rows,
            "mom_radiated_fraction": total_mom,
            "forward_radiated_fraction_before_renormalising": total_fwd,
            "worst_order_error_percentage_points": max(r["error_percentage_points"] for r in rows.values()),
            "forward_power_in_evanescent_orders_percent": 100.0 * (1.0 - total_fwd / total_mom),
        }
        print(
            f"[mom {tag}{name}] {out['cases'][name]['seconds']}s worst order error "
            f"{out['cases'][name]['worst_order_error_percentage_points']:.1f} pp, forward model radiates "
            f"{total_fwd:.3f} of the incident power (MoM {total_mom:.3f})",
            flush=True,
        )
    return out


def main(out_path, quick=False, skip_mom=False, mom_only=False):
    grid = (
        HeldOutGrid(theta_step_deg=2.0, phi_step_deg=4.0, cut_step_u=0.01, patch_step_deg=0.25)
        if quick
        else HeldOutGrid()
    )
    only = None
    if mom_only:
        result = json.loads(Path(out_path).read_text())
        only = [int(v) for v in mom_only.split(",")] if isinstance(mom_only, str) else None
        if only is None:
            result.pop("part_two_plane_wave_supercell_mom", None)
    else:
        t0 = time.time()
        one = part_one(grid, quick)
        print(f"part I done in {time.time() - t0:.0f}s", flush=True)
        result = {
            "grid": "quick (2x4 deg)" if quick else "held-out 0.5 deg x 1 deg",
            "labels": ["MODEL_ONLY", "NOT_HFSS_VALIDATED"],
            "part_one_pattern_path": one,
        }
    Path(out_path).write_text(json.dumps(result, indent=1))
    if not skip_mom:
        res = _synthesize()
        for mesh, n in ((12, 5),) if quick else MOM_CASES:
            if only is not None and mesh not in only:
                continue
            result.setdefault("part_two_plane_wave_supercell_mom", {})[f"mesh{mesh}"] = part_two(res, mesh, n)
            Path(out_path).write_text(json.dumps(result, indent=1))
    result["part_two_surface_wave_fed_finite_board_mom"] = "NOT_AVAILABLE_FROM_CURRENT_DATA"
    Path(out_path).write_text(json.dumps(result, indent=1))
    print(json.dumps({k: v for k, v in result.items() if k != "part_one_pattern_path"}, indent=1)[:3000])


if __name__ == "__main__":
    main(
        sys.argv[1],
        quick="--quick" in sys.argv,
        skip_mom="--skip-mom" in sys.argv,
        mom_only=next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mom-only=")), "--mom-only" in sys.argv),
    )
