"""End-to-end Part C design: phase synthesis -> psi0 -> direct per-cell radius solve -> predicted pattern."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import design as dz
from . import pattern as pt
from .cells import CellFamilyModel, wrap_deg
from .feed import field_amplitude, q_for_edge_taper
from .layout import LensLayout
from .synthesis import (
    ApertureFFT,
    beam_windows,
    binary_baseline_phase,
    gerchberg_saxton,
    peak_u_along_x,
    window_efficiency,
)

C0 = 299_792_458.0


@dataclass(frozen=True)
class DesignConfig:
    freq_hz: float
    layout: LensLayout
    beam_deg: float
    edge_taper_db: float
    model: CellFamilyModel
    min_feature_m: float
    d_seam_deg: float
    cap_fill: float
    guess_library_csv: str
    window_lambda_over_d: tuple = (0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9)  # candidate GS window radii; best pointing wins
    gs_iterations: int = 500
    workers: int = 4


def synthesize(cfg: DesignConfig, cell_amp=None):
    """Dual-beam aperture phase (GS) and the transmission phase base arg T = psi + k r (psi0 not yet applied)."""
    lay = cfg.layout
    lam = C0 / cfg.freq_hz
    k0 = 2 * math.pi / lam
    ix, iy, x, y = lay.cells()
    th, phi, r = lay.incidence(x, y)
    q = q_for_edge_taper(math.degrees(math.atan(lay.d_flats_m / 2 / lay.focal_m)), cfg.edge_taper_db)
    amp = field_amplitude(th, r, q) * (1.0 if cell_amp is None else cell_amp)
    amp = amp / np.sqrt(np.sum(amp**2))
    fft = ApertureFFT(ix, iy, lay.n_grid, lam, lay.pitch_m, 256)
    groups, _sign = quadrant_map(x, y, lay.pitch_m)
    ideal = np.sum(amp) ** 2
    u0 = math.sin(math.radians(cfg.beam_deg))
    best = None
    for w in cfg.window_lambda_over_d:
        cand, _eff = gerchberg_saxton(
            amp, x, fft, cfg.beam_deg, w * lam / lay.d_flats_m, cfg.gs_iterations, k0, shaped=True, mirror_groups=groups
        )
        up, f = peak_u_along_x(amp, cand, x, k0, u0, 0.06, 481)
        key = (round(abs(math.degrees(math.asin(up)) - cfg.beam_deg), 2), -f)  # pointing first, then level
        if best is None or key < best[0]:
            best = (key, cand, w)
    psi, window_used = best[1], best[2]
    bin_psi = binary_baseline_phase(x, k0, cfg.beam_deg)
    report_masks = beam_windows(fft, cfg.beam_deg, 1.2 * lam / lay.d_flats_m)

    def describe(ps):
        up, f = peak_u_along_x(amp, ps, x, k0, u0, 0.06, 481)
        _p, tot = window_efficiency(fft.forward(amp * np.exp(1j * ps)), report_masks)
        return dict(
            peak_angle_deg=math.degrees(math.asin(up)),
            peak_over_ideal_single=float(f / ideal),
            window_1p2_lambda_over_d_efficiency=tot,
        )

    gs = dict(
        window_lambda_over_d_used=window_used,
        synthesised=describe(psi),
        binary_baseline=describe(bin_psi),
        ideal_equal_split_peak_over_ideal_single=0.5,
    )
    return dict(
        x=x,
        y=y,
        theta=th,
        phi=phi,
        r=r,
        q=q,
        k0=k0,
        amp=amp,
        psi=psi,
        base_phase_deg=np.degrees(psi + k0 * r),
        weight=field_amplitude(th, r, q) ** 2 * np.cos(np.radians(th)),
        gs=gs,
    )


def quadrant_map(x, y, pitch):
    """Index of the first-quadrant representative of every cell and the mirror-sign of its cross-pol."""
    key = {}
    rep = np.empty(len(x), int)
    for i in range(len(x)):
        k = (round(abs(x[i]) / pitch * 2), round(abs(y[i]) / pitch * 2))
        if k not in key:
            key[k] = i
        rep[i] = key[k]
    sign = np.where((x < 0) ^ (y < 0), -1.0, 1.0)
    return rep, sign


def run_design(cfg: DesignConfig, cell_amp=None, psi0_deg=None):
    syn = synthesize(cfg, cell_amp)
    x, y, th, phi = syn["x"], syn["y"], syn["theta"], syn["phi"]
    lay = cfg.layout
    rep, sign = quadrant_map(x, y, lay.pitch_m)
    # the lens is mirror-symmetric; verify the synthesised target is too, then solve one quadrant
    base = syn["base_phase_deg"]
    dev = np.abs(wrap_deg(base - base[rep]))
    syn["symmetry_dev_deg"] = float(dev.max())
    base_sym = base[rep]
    limits = dz.FamilyLimits.build(lay.pitch_m, cfg.min_feature_m, cfg.cap_fill)
    rows = dz.load_guess_library(cfg.guess_library_csv)
    table = dz.guess_table(rows, cfg.cap_fill)
    guesses = [dz.cell_guess(table, t, p) for t, p in zip(th, phi, strict=True)]
    if psi0_deg is None:
        best = dz.choose_psi0(guesses, base_sym, syn["weight"], cfg.d_seam_deg)
        if best is None:
            raise RuntimeError("no global phase offset gives every cell a library guess within tolerance")
        psi0_deg = best[0]
    target = base_sym + psi0_deg
    uniq = np.unique(rep)
    sol = dz.assign_cells(
        cfg.model, limits, [guesses[i] for i in uniq], th[uniq], phi[uniq], target[uniq], cfg.d_seam_deg, cfg.workers
    )
    sol = dz.resolve_empty_cells(cfg.model, sol, th[uniq], phi[uniq], target[uniq], limits.pillar_min_m)
    by_rep = {int(i): s for i, s in zip(uniq, sol, strict=True)}
    cells = []
    for i in range(len(x)):
        s = by_rep[int(rep[i])]
        j = s["J"].copy()
        j[0, 1] *= sign[i]
        j[1, 0] *= sign[i]
        cells.append(dict(kind=s["kind"], R=s["R"], residual_deg=s["residual_deg"], J=j, n_solves=s["n_solves"]))
    return dict(syn=syn, limits=limits, psi0_deg=psi0_deg, target_deg=target, cells=cells, rep=rep, sign=sign)


def predict_pattern(cfg: DesignConfig, result, theta_deg=None, phi_deg=None):
    syn = result["syn"]
    lay = cfg.layout
    if theta_deg is None:
        theta_deg = np.arange(0.0, 90.01, 0.25)
    if phi_deg is None:
        phi_deg = np.arange(0.0, 360.0, 2.0)
    ex, ey, tau = pt.exit_field(
        lay, syn["q"], syn["x"], syn["y"], syn["theta"], syn["r"], syn["k0"], [c["J"] for c in result["cells"]]
    )
    u_inc = (syn["x"] / syn["r"], syn["y"] / syn["r"])
    ff = pt.far_field(syn["x"], syn["y"], ex, ey, syn["k0"], lay.pitch_m, theta_deg, phi_deg, u_inc=u_inc)
    p_in = pt.incident_power_per_cell(syn["theta"], syn["r"], syn["q"], lay.pitch_m)
    p_ap = float(np.sum(p_in * tau))
    metrics = pt.beam_metrics(ff, cfg.beam_deg)
    metrics.update(
        aperture_power_fraction_of_horn=p_ap,
        aperture_transmission_mean=float(np.sum(p_in * tau) / np.sum(p_in)),
        radiated_power_fraction_of_horn=ff.radiated_power(),
    )
    return ff, metrics, ex, ey, tau
