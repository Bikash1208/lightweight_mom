"""Phase 28 (literature step, model side): scalar 1-/2-bit cells with a SIDE (plane-like) feed against the centre monopole, and grouped control cells. IDEALISED_CELLS, MODEL_ONLY.

    python scripts/research/phase28_feed_grouping_model.py <out.json> [--workers N]

Same aperture as the item-2 reference (150 mm circumradius hexagon, 6.4 mm lattice), the 10.0 GHz supercell table, `coded_feasibility` pattern code unchanged.
* Centre monopole, scalar cell: exactly the item-2 `lp` variant (radial current direction, amplitude/phase from the Phase 26 march, eta_conv from the table).
* Side feed (ASSUMPTION, not from any paper): a plane-like TM0 surface wave enters at the hexagon edge on the -x side and travels along +x, uniform across y, phase exp(-j beta (x - x_in(y))) with
  beta = 1.053 k0 (the table mean), amplitude exp(-alpha (x - x_in(y))) with ONE uniform leakage alpha chosen so that the converted fraction equals the monopole case (0.504) for this hexagon; the cell current
  is x-directed (the TM0 field is along the propagation direction): J = A exp(-j phi_ref) m x_hat, m = cos(phi_ref - k_t . r) quantised to 2**bits equally spaced levels in [-1, 1].
* Grouping (ASSUMPTION: group sizes are not from [78]): g x g lattice cells share one state, sampled at the group centre; every cell keeps its own reference phase and amplitude.
Gain bookkeeping as item 2: D_ap + 10 log10(0.592 eta_conv) (launch fraction of the 10.54 GHz bare slab, held; for the side feed this is a placeholder, only the pattern terms are compared).
"""

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import json  # noqa: E402
import math  # noqa: E402
import sys  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scipy.optimize import brentq  # noqa: E402

from src.hologram import coded_feasibility as cf  # noqa: E402
from src.hologram.tensor_leaky_cp import SupercellTable  # noqa: E402

cf.set_frequency(10.0e9)
TABLE_GLOB = str(ROOT / "results/phase28/f10/supercell_table_*.json")
R_MM, PITCH = 150.0, 6.4e-3
BETA_OVER_K0 = 1.053
BEAMS = [(0.0, 0.0), (30.0, 0.0), (30.0, 90.0), (45.0, 0.0), (45.0, 90.0)]
GROUPS = [1, 2, 3, 4]


def hex_chord_x(poly, y):
    """x-interval of the convex polygon at height `y` (array)."""
    xs_lo = np.full(y.shape, np.inf)
    xs_hi = np.full(y.shape, -np.inf)
    for a, b in zip(poly, np.roll(poly, -1, axis=0), strict=True):
        lo, hi = sorted((a[1], b[1]))
        if hi - lo < 1e-12:
            continue
        ok = (y >= lo - 1e-12) & (y <= hi + 1e-12)
        x = a[0] + (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1])
        xs_lo = np.where(ok, np.minimum(xs_lo, x), xs_lo)
        xs_hi = np.where(ok, np.maximum(xs_hi, x), xs_hi)
    return xs_lo, xs_hi


def side_alpha(poly, eta_target):
    """Uniform alpha [Np/m] whose converted fraction, averaged over the y-uniform incident flux, equals `eta_target`."""
    y = np.linspace(-np.abs(poly[:, 1]).max() * 0.9999, np.abs(poly[:, 1]).max() * 0.9999, 4001)
    lo, hi = hex_chord_x(poly, y)
    L = hi - lo
    w = np.maximum(L, 0) > 0
    # flux per unit y is uniform over the entry edge at -x: the hexagon's left edge is slanted, the flux through it per unit y is the same
    f = lambda a: float(np.mean(1.0 - np.exp(-2.0 * a * L[w]))) - eta_target  # noqa: E731
    return brentq(f, 1e-3, 200.0)


def twin_lobe(jx, jy, x, y, green, beam, peak_tp):
    th, ph = np.meshgrid(np.arange(0.0, 90.0 + 1e-9, 2.0), np.arange(0.0, 360.0, 4.0), indexing="ij")
    et, ep = cf.far_field(jx, jy, x, y, green, th, ph, PITCH)
    pw = np.abs(et) ** 2 + np.abs(ep) ** 2
    d = cf._angle_deg(th, ph, peak_tp[0], peak_tp[1])
    far = np.where(d > 15.0, pw, -1.0)
    i = np.unravel_index(np.argmax(far), far.shape)
    return float(th[i]), float(ph[i]), float(far[i])


def case(job):
    feed, bits, beam, g = job
    table = SupercellTable.from_files(TABLE_GLOB)
    poly = cf.hexagon_m(R_MM * 1e-3)
    green = cf.green_function()
    th0, ph0 = math.radians(beam[0]), math.radians(beam[1])
    if feed == "monopole":
        x, y, d, march = cf.design_cells(table, poly, PITCH, beam[0], beam[1])
        eta = float(march["eta_conv"])
        m = cf.lp_modulation(d, x, y, *beam)
        grp_ids = None
        idx, mq = cf.lp_states(m, bits)
        jx, jy = cf.lp_currents(d, mq)
    else:
        x, y = cf.lattice(poly, PITCH, keepout_m=0.0)
        eta = 0.504  # same conversion as the monopole case (see module docstring)
        alpha = side_alpha(poly, eta)
        lo, hi = hex_chord_x(poly, y)
        beta = BETA_OVER_K0 * cf.K0
        s = x - lo
        amp = np.exp(-alpha * s)
        phi_ref = beta * s
        kt = cf.K0 * math.sin(th0)
        # sampling point of the state: the lattice cell itself (g = 1) or the centre of its g x g block
        i_idx, j_idx = np.round(x / PITCH).astype(int), np.round(y / PITCH).astype(int)
        gi, gj = i_idx // g, j_idx // g
        cx, cy = (gi * g + (g - 1) / 2.0) * PITCH, (gj * g + (g - 1) / 2.0) * PITCH
        lo_c, _ = hex_chord_x(poly, np.clip(cy, -np.abs(poly[:, 1]).max() * 0.9999, np.abs(poly[:, 1]).max() * 0.9999))
        phi_c = beta * (cx - lo_c)
        m = np.cos(phi_c - kt * (cx * math.cos(ph0) + cy * math.sin(ph0)))
        idx, mq = cf.lp_states(m, bits)
        j = amp * np.exp(-1j * phi_ref) * mq
        jx, jy = j, np.zeros_like(j)
        grp_ids = gi * 100000 + gj
    r = cf.analyse(jx, jy, x, y, green, PITCH, "lp", beam)
    r["gain_dbi"] = cf.launch_gain_dbi(r["d_ap_dbi"], eta)
    tl = twin_lobe(jx, jy, x, y, green, beam, (r["peak_theta_deg"], r["peak_phi_deg"]))
    n_groups = int(len(np.unique(grp_ids))) if grp_ids is not None else int(x.size)
    return {
        "feed": feed,
        "bits": bits,
        "beam": list(beam),
        "group": g,
        "n_cells": int(x.size),
        "n_groups": n_groups,
        "eta_conv": eta,
        "d_ap_dbi": r["d_ap_dbi"],
        "gain_dbi": r["gain_dbi"],
        "peak_theta_deg": r["peak_theta_deg"],
        "peak_phi_deg": r["peak_phi_deg"],
        "pointing_error_deg": r["pointing_error_deg"],
        "sll_db": r["sll_db"],
        "twin_theta_phi_deg": [tl[0], tl[1]],
        "twin_level_db_coarse": 10 * math.log10(max(tl[2], 1e-30) / max(r["radiated_power_arb"], 1e-30)) if False else None,
    }


def main(out, workers):
    jobs = [("monopole", b, beam, 1) for b in (1, 2) for beam in BEAMS]
    jobs += [("side", b, beam, g) for b in (1, 2) for beam in BEAMS for g in GROUPS]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        rows = list(ex.map(case, jobs))
    doc = {
        "labels": ["MODEL_ONLY", "IDEALISED_CELLS", "NOT_HFSS_VALIDATED"],
        "frequency_ghz": 10.0,
        "assumptions": "side feed, uniform leakage and group sizes are ASSUMPTIONS (see script docstring); launch fraction 0.592 held",
        "rows": rows,
    }
    Path(out).write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print("wrote", out, len(rows))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], int(a[a.index("--workers") + 1]) if "--workers" in a else 24)
