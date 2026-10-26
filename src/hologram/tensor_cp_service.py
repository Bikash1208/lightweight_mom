"""Phase 27 item 3: the "Tensor CP (monopole-fed)" design request as a plain function (the app router wraps it in a job).

`run_design(request, out_dir=None, progress=None, cancelled=None)` -> JSON-able dict. Boards: `single` (one 150 mm hexagon),
`hex4_static` (four hexagons, one beam per port) and `hex4_combined` (four hexagons, one beam, equal amplitude and phase).
Requests outside the registered family's domain (substrate, frequency, period) return the named state
`NOT_AVAILABLE_FROM_CURRENT_DATA` with the reasons -- nothing is extrapolated. Every number carries its source label; the
registry's known biases ride along. No physics lives here: it calls `tensor_cp_design` and `cell_family_registry`.

LHCP: the design uses `hand = -1`; metrics are those of the RHCP code with `E_phi -> -E_phi` (swaps the circular components).
Pin height: the launch fraction scales with the thin-wire model (`dispersion.pin_launch`, effective radius 1.5 mm, no top
hat) relative to the HFSS-calibrated 6.0 mm value (0.62); other heights are MODEL_ONLY.
"""

import numpy as np

from ..dispersion.pin_launch import Launcher, PinLaunch
from ..pattern.vector_aperture_pattern import circular_components, copol_component
from . import cell_family_registry as registry
from . import tensor_cp_design as t
from .tensor_leaky_cp import SupercellTable, cell_design, march_rays

FAMILY = "tensor_slotted_circle_rt5880_6p4"
BOARDS = ("single", "hex4_static", "hex4_combined")
PIN_CALIBRATED_MM = 6.0
HFSS_ENVELOPE_THETA_DEG = 30.0
S_EDGE = (1.0, 1.7)
_R_EFF_M = 1.5e-3


def launch_fraction(pin_top_above_surface_mm):
    """`(launch fraction, label)`: 0.62 at the calibrated 6.0 mm; other heights scaled by the thin-wire TM0 share."""
    if abs(pin_top_above_surface_mm - PIN_CALIBRATED_MM) < 1e-9:
        return t.LAUNCH, "HFSS_CALIBRATED"

    def tm0(h_mm):
        launcher = Launcher(h_mm * 1e-3, _R_EFF_M, 0.0)
        return float(PinLaunch(t.F, 2.2, 1.575e-3, launcher, n_seg=10).split()["launch_fraction"])

    return float(t.LAUNCH * tm0(pin_top_above_surface_mm) / tm0(PIN_CALIBRATED_MM)), "MODEL_ONLY"


def default_request(board="single"):
    beams = {
        "single": [[30.0, 225.0]],
        "hex4_static": [[30.0, 45.0], [30.0, 135.0], [30.0, 225.0], [30.0, 315.0]],
        "hex4_combined": [[30.0, 225.0]],
    }[board]
    return {
        "board": board,
        "beams_deg": beams,
        "handedness": "rhcp",
        "pin_height_mm": PIN_CALIBRATED_MM,
        "frequency_ghz": 10.54,
        "substrate": {"eps_r": 2.2, "tan_d": 0.0009, "h_mm": 1.575},
        "period_mm": 6.4,
        "quick": False,
    }


def normalise_request(req):
    r = default_request(req.get("board", "single"))
    for k, v in req.items():
        if k == "substrate":
            r["substrate"] = {**r["substrate"], **v}
        elif v is not None:
            r[k] = v
    n_need = 4 if r["board"] == "hex4_static" else 1
    if len(r["beams_deg"]) == 1 and r["board"] == "hex4_static":
        raise ValueError("hex4_static needs four beams (one per port)")
    if len(r["beams_deg"]) != n_need:
        raise ValueError(f"board {r['board']!r} needs {n_need} beam(s), got {len(r['beams_deg'])}")
    if r["board"] not in BOARDS:
        raise ValueError(f"board must be one of {BOARDS}")
    if r["handedness"] not in ("rhcp", "lhcp"):
        raise ValueError("handedness must be 'rhcp' or 'lhcp'")
    return r


def _co_flip(et, ep, hand):
    """Fields whose RHCP is the requested co-polar component (LHCP: swap the circular components)."""
    return (et, ep) if hand > 0 else (et, -ep)


def _sll_db(rh, theta, phi, peak):
    tg, pg = np.meshgrid(theta, phi, indexing="ij")
    i = peak
    cd = np.cos(np.radians(tg)) * np.cos(np.radians(tg[i])) + np.sin(np.radians(tg)) * np.sin(
        np.radians(tg[i])
    ) * np.cos(np.radians(pg) - np.radians(pg[i]))
    d = np.degrees(np.arccos(np.clip(cd, -1, 1)))
    return float(10 * np.log10(rh[d > 15.0].max() / rh[i]))


def _cuts(et, ep, theta, phi, beam, peak_gain_scale):
    """Co/cross gain (dBi, anchored) along the phi-cut through the beam over theta -90..90, and the AR along it."""
    ph0 = beam[1]
    ip = int(np.argmin(np.abs(((phi - ph0 + 180.0) % 360.0) - 180.0)))
    ip2 = int(np.argmin(np.abs(((phi - (ph0 + 180.0) + 180.0) % 360.0) - 180.0)))
    er, el = circular_components(et, ep)
    co = np.concatenate(
        [np.abs(copol_component(et, ep, "rhcp")[::-1, ip2]) ** 2, np.abs(copol_component(et, ep, "rhcp")[:, ip]) ** 2]
    )
    cr = np.concatenate(
        [np.abs(copol_component(et, ep, "lhcp")[::-1, ip2]) ** 2, np.abs(copol_component(et, ep, "lhcp")[:, ip]) ** 2]
    )
    th = np.concatenate([-theta[::-1], theta])
    a, b = np.sqrt(co), np.sqrt(cr)
    ar = 20 * np.log10((a + b) / np.maximum(np.abs(a - b), 1e-12))
    g = lambda v: (10 * np.log10(np.maximum(v * peak_gain_scale, 1e-12))).tolist()  # noqa: E731
    return {
        "phi_cut_deg": [float(phi[ip]), float(phi[ip2])],
        "theta_deg": th.tolist(),
        "co_gain_dbi": g(co),
        "cross_gain_dbi": g(cr),
        "ar_db": np.minimum(ar, 60.0).tolist(),
    }


def _save_png(path, fn):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = fn(plt)
        fig.savefig(path, dpi=100)
        plt.close(fig)
        return path.name
    except Exception as exc:  # noqa: BLE001 -- figures are a convenience, never fatal
        return f"NOT_WRITTEN: {type(exc).__name__}: {exc}"


def _cell_map_png(path, cells, feeds, title):
    def fn(plt):
        fig, ax = plt.subplots(figsize=(8, 8))
        col = {1: "#d9a441", 2: "#6a9fd1", 3: "#8cc08c", 4: "#d98a8a"}
        for c in cells:
            ax.add_patch(
                plt.Circle((c["x_mm"], c["y_mm"]), t_D / 2, fc=col.get(c.get("port", 1), "#d9a441"), ec="none")
            )
            a = np.radians(c["slot_angle_deg"])
            dx, dy = 0.5 * c["slot_L_mm"] * np.cos(a), 0.5 * c["slot_L_mm"] * np.sin(a)
            ax.plot([c["x_mm"] - dx, c["x_mm"] + dx], [c["y_mm"] - dy, c["y_mm"] + dy], "w-", lw=0.8)
        for f in feeds:
            ax.plot(*f, "r+", ms=12)
        ax.set_aspect("equal")
        ax.autoscale_view()
        ax.set_title(title)
        return fig

    t_D = 6.0
    return _save_png(path, fn)


def _pattern_png(path, cut, title):
    def fn(plt):
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        a1.plot(cut["theta_deg"], cut["co_gain_dbi"], label="co-pol")
        a1.plot(cut["theta_deg"], cut["cross_gain_dbi"], label="cross-pol")
        a1.set_ylim(-10, max(cut["co_gain_dbi"]) + 3)
        a1.set_ylabel("gain [dBi] (anchored)")
        a1.legend()
        a1.set_title(title)
        a2.plot(cut["theta_deg"], cut["ar_db"])
        a2.set_ylim(0, 20)
        a2.set_ylabel("axial ratio [dB]")
        a2.set_xlabel("theta [deg] along the phi-cut through the beam")
        return fig

    return _save_png(path, fn)


def _row(port, beam, m, gain, label_gain, extra=None):
    row = {
        "port": port,
        "beam_theta_phi_deg": list(beam),
        "peak_theta_phi_deg": [m["peak_theta_deg"], m["peak_phi_deg"]],
        "gain_rhcp_or_lhcp_dbi": gain,
        "ar_peak_db": m["ar_peak_db"],
        "ar_worst_in_minus3db_beam_db": m["ar_worst_in_minus3db_beam_db"],
        "xpd_db": m["xpd_db"],
        "sll_db": m["sll_db"],
        "sources": {
            "gain": label_gain,
            "ar_xpd": "MODEL_ONLY (pessimistic vs one HFSS run)",
            "sll": "MODEL_ONLY (optimistic by about 4 dB vs one HFSS run)",
        },
    }
    row.update(extra or {})
    return row


def _single(r, table, launch, progress):
    beam = r["beams_deg"][0]
    hand = +1 if r["handedness"] == "rhcp" else -1
    th0, ph0 = np.radians(beam[0]), np.radians(beam[1])
    poly = t.hexagon()
    x, y = t.lattice_cells(poly)
    green = t.LayeredGreenFunction(t.SUB, t.F, backing="pec")
    m = march_rays(table, poly, (0.0, 0.0), t.K0, th0, ph0, 0.9, taper_c=0.5, n_psi=360 if r["quick"] else 720)
    d = cell_design(m, table, x, y, (0.0, 0.0), t.K0, th0, ph0, hand=hand)
    theta = np.arange(0.0, 90.0 + 1e-9, 2.0 if r["quick"] else 1.0)
    phi = np.arange(0.0, 360.0, 3.0 if r["quick"] else 2.0)
    progress(0.4, "pattern")
    et, ep = t.pattern(d, x, y, green, theta, phi)
    et, ep = _co_flip(et, ep, hand)
    mt = t.metrics(et, ep, theta, phi, beam=tuple(beam))
    eta = m["eta_conv"]
    gain = mt["directivity_ap_rhcp_dbi"] + 10 * np.log10(launch * eta) + t.OFFSET_DB
    rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
    pk = np.unravel_index(np.argmax(rh), rh.shape)
    mm = {
        "peak_theta_deg": mt["peak_rhcp_theta_deg"],
        "peak_phi_deg": mt["peak_rhcp_phi_deg"],
        "ar_peak_db": mt["ar_at_peak_db"],
        "ar_worst_in_minus3db_beam_db": mt["ar_max_over_minus3db_beam_db"],
        "xpd_db": mt["xpd_at_peak_db"],
        "sll_db": mt["rhcp_sll_db"],
    }
    scale = 4 * np.pi * 10 ** (t.OFFSET_DB / 10) * launch * eta / t.prad(et, ep, theta, phi)
    lf = launch
    budget = {
        "pin_direct_space_wave_assumed_outside_beam": 1 - lf,
        "unconverted_reaches_edge": lf * (1 - eta),
        "aperture_radiated": lf * eta,
        "main_beam_8deg": lf * eta * mt["within_8deg_fraction"],
        "lhcp_in_aperture_pattern": lf * eta * mt["lhcp_fraction"],
        "grazing_plus_back_upper_bound_unmodelled": (1 - lf) + lf * (1 - eta),
        "source": "ACCEPTED-power split: aperture terms MODEL_ONLY; the pin-direct share is the launch model (HFSS single hexagon: non-beam 0.70 vs 0.68 modelled in total)",
    }
    cells = [
        {
            "x_mm": float(x[i] * 1e3),
            "y_mm": float(y[i] * 1e3),
            "slot_L_mm": float(d["slot_L_mm"][i]),
            "slot_W_mm": float(d["slot_W_mm"][i]),
            "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
            "alpha_np_per_m": float(d["alpha"][i]),
            "port": 1,
        }
        for i in range(x.size)
    ]
    row = _row(1, beam, mm, float(gain), "HFSS_ANCHORED" if lf == t.LAUNCH else "MODEL_ONLY", {"eta_conv": eta})
    del pk
    return {
        "ports": [row],
        "budget": budget,
        "cells": cells,
        "feeds": [(0.0, 0.0)],
        "cut": _cuts(et, ep, theta, phi, beam, scale),
        "eta_conv": eta,
    }


def _hex4(r, table, launch, progress):
    hand = +1 if r["handedness"] == "rhcp" else -1
    board = "A" if r["board"] == "hex4_static" else "B"
    beams = (
        {j: tuple(r["beams_deg"][j - 1]) for j in t.FEEDS}
        if board == "A"
        else {j: tuple(r["beams_deg"][0]) for j in t.FEEDS}
    )
    bd = t.design_board(board, table, beams=beams, hand=hand)
    green = t.LayeredGreenFunction(t.SUB, t.F, backing="pec")
    theta = np.arange(0.0, 90.0 + 1e-9, 3.0 if r["quick"] else 1.0)
    phi = np.arange(0.0, 360.0, 6.0 if r["quick"] else 1.0)
    edges = t.crossing_edges(bd)
    rows, out_by_s = [], {}
    scale_launch = launch / t.LAUNCH
    for si, s in enumerate(S_EDGE):
        progress(0.15 + 0.4 * si, f"far fields s={s}")
        ff = t.far_fields(bd, s, green, theta, phi)
        fields = {j: t.normalised(ff[j], theta, phi) for j in t.FEEDS}
        out_by_s[s] = fields
    cells = []
    for j in t.FEEDS:
        d = bd[j]["design"]
        for i in range(bd[j]["x_mm"].size):
            cells.append(
                {
                    "x_mm": float(bd[j]["x_mm"][i]),
                    "y_mm": float(bd[j]["y_mm"][i]),
                    "slot_L_mm": float(d["slot_L_mm"][i]),
                    "slot_W_mm": float(d["slot_W_mm"][i]),
                    "slot_angle_deg": float(np.degrees(d["slot_angle_rad"][i])),
                    "alpha_np_per_m": float(d["alpha"][i]),
                    "port": int(j),
                }
            )
    gain_shift = 10 * np.log10(scale_launch)

    def metrics_of(et, ep, beam):
        et, ep = _co_flip(et, ep, hand)
        m = t.pol_metrics(et, ep, theta, phi, beam)
        rh = np.abs(copol_component(et, ep, "rhcp")) ** 2
        tg, pg = np.meshgrid(theta, phi, indexing="ij")
        near = (
            np.degrees(
                np.arccos(
                    np.clip(
                        np.cos(np.radians(tg)) * np.cos(np.radians(beam[0]))
                        + np.sin(np.radians(tg))
                        * np.sin(np.radians(beam[0]))
                        * np.cos(np.radians(pg) - np.radians(beam[1])),
                        -1,
                        1,
                    )
                )
            )
            <= 12.0
        )
        i = np.unravel_index(np.argmax(np.where(near, rh, -1)), rh.shape)
        m["sll_db"] = _sll_db(rh, theta, phi, i)
        return m, et, ep, i

    cut = None
    if board == "A":
        for j in t.FEEDS:
            ms = {}
            for s in S_EDGE:
                et, ep, rr = out_by_s[s][j]
                ms[s] = metrics_of(et, ep, beams[j])
            m1, m17 = ms[S_EDGE[0]][0], ms[S_EDGE[1]][0]
            g1, g17 = m1["rhcp_gain_dbi"] + gain_shift, m17["rhcp_gain_dbi"] + gain_shift
            rows.append(
                _row(
                    j,
                    beams[j],
                    m1,
                    g1,
                    "HFSS_ANCHORED" if scale_launch == 1.0 else "MODEL_ONLY",
                    {
                        "gain_range_over_edge_power_scale_dbi": [min(g1, g17), max(g1, g17)],
                        "edge_power_scale": list(S_EDGE),
                        "pickup_over_main_power": float(out_by_s[S_EDGE[0]][j][2]),
                    },
                )
            )
        j0 = 3 if 3 in beams else 1
        et, ep, _ = out_by_s[S_EDGE[0]][j0]
        m, e2, p2, i = metrics_of(et, ep, beams[j0])
        sc = 4 * np.pi * 10 ** (t.OFFSET_DB / 10) * scale_launch
        cut = _cuts(e2, p2, theta, phi, beams[j0], sc)
        cut["port"] = j0
    else:
        beam = beams[1]
        for s in S_EDGE:
            et = sum(0.5 * out_by_s[s][j][0] for j in t.FEEDS)
            ep = sum(0.5 * out_by_s[s][j][1] for j in t.FEEDS)
            m, e2, p2, i = metrics_of(et, ep, beam)
            if s == S_EDGE[0]:
                m1, et1, ep1 = m, e2, p2
            g = m["rhcp_gain_dbi"] + gain_shift
            if s == S_EDGE[0]:
                g1 = g
            else:
                g17 = g
        rows.append(
            _row(
                0,
                beam,
                m1,
                g1,
                "HFSS_ANCHORED" if scale_launch == 1.0 else "MODEL_ONLY",
                {
                    "gain_range_over_edge_power_scale_dbi": [min(g1, g17), max(g1, g17)],
                    "edge_power_scale": list(S_EDGE),
                    "ports": "all four, 0.25 W each, equal phase",
                },
            )
        )
        sc = 4 * np.pi * 10 ** (t.OFFSET_DB / 10) * scale_launch
        cut = _cuts(et1, ep1, theta, phi, beam, sc)
        ports_alone = []
        for j in t.FEEDS:
            et, ep, rr = out_by_s[S_EDGE[0]][j]
            m, *_ = metrics_of(et, ep, beam)
            ports_alone.append({"port": j, "gain_dbi": m["rhcp_gain_dbi"] + gain_shift})
        rows[0]["each_port_alone_dbi"] = ports_alone
    return {
        "ports": rows,
        "budget": {
            "edges_unconverted_power_fraction_of_launched_wave": {str(j): v for j, v in edges.items()},
            "note": "x launch fraction for the accepted-power share; shared edges feed the neighbours' cells (pickup), the free edge diffracts",
            "source": "MODEL_ONLY",
        },
        "cells": cells,
        "feeds": [t.FEEDS[j] for j in t.FEEDS],
        "cut": cut,
        "eta_conv": float(np.mean([bd[j]["march"]["eta_conv"] for j in t.FEEDS])),
    }


def run_design(request, out_dir=None, progress=None, cancelled=None):
    progress = progress or (lambda frac, msg: None)
    cancelled = cancelled or (lambda: False)
    r = normalise_request(request)
    sub = r["substrate"]
    check = registry.check_request(
        FAMILY, substrate=(sub["eps_r"], sub["tan_d"], sub["h_mm"]), f_ghz=r["frequency_ghz"], period_mm=r["period_mm"]
    )
    fam = registry.get(FAMILY)
    base = {
        "request": r,
        "family": FAMILY,
        "validation_label": fam.validation_label,
        "known_biases": list(fam.known_biases),
        "registry_check": check.to_dict(),
        "labels": ["MODEL_ONLY", "SOLVER_SUPERCELL"],
    }
    if check.state != registry.AVAILABLE:
        return {**base, "state": registry.NOT_AVAILABLE, "reasons": list(check.reasons)}
    warnings = []
    for b in r["beams_deg"]:
        if abs(b[0] - HFSS_ENVELOPE_THETA_DEG) > 0.5:
            warnings.append(
                f"BEAM_OUTSIDE_HFSS_ENVELOPE: theta {b[0]:g} deg; the HFSS-checked design is theta 30 only (one run)"
            )
    if r["handedness"] == "lhcp":
        warnings.append(
            "LHCP_NOT_HFSS_CHECKED: the hand = -1 design is the mirror-sense of the checked RHCP one (symmetry test only)"
        )
    if r["quick"]:
        warnings.append(
            "QUICK_GRID: coarse pattern grid undersamples narrow beams (the four-hexagon beam loses up to 4 dB at the peak); not for reporting"
        )
    launch, launch_label = launch_fraction(r["pin_height_mm"])
    if launch_label != "HFSS_CALIBRATED":
        warnings.append(
            f"PIN_HEIGHT_{r['pin_height_mm']:g}_MM: launch fraction {launch:.3f} is MODEL_ONLY (calibrated only at 6.0 mm)"
        )
    progress(0.02, "loading supercell table")
    table = SupercellTable.from_files(t.TABLE_GLOB)
    res = _single(r, table, launch, progress) if r["board"] == "single" else _hex4(r, table, launch, progress)
    if cancelled():
        return {**base, "state": "CANCELLED"}
    files = {}
    if out_dir is not None:
        from pathlib import Path

        out_dir = Path(out_dir)
        files["cell_map_png"] = _cell_map_png(
            out_dir / "cell_map.png", res["cells"], res["feeds"], f"{r['board']}: {len(res['cells'])} slotted circles"
        )
        files["pattern_png"] = _pattern_png(
            out_dir / "pattern_cut.png", res["cut"], f"{r['board']} {r['handedness']} (anchored, MODEL_ONLY)"
        )
    return {
        **base,
        "state": "OK",
        "launch_fraction": launch,
        "launch_source": launch_label,
        "warnings": warnings,
        "n_cells": len(res["cells"]),
        "eta_conv": res["eta_conv"],
        "ports": res["ports"],
        "power_budget": res["budget"],
        "pattern_cut": res["cut"],
        "cells": res["cells"],
        "files": files,
    }
