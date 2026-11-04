"""Cell search for a larger switched-cell swing (alpha_beam 5-10 Np/m target). Generic crossed two-diode disc cell on any stack: disc diameter, slit width g, slit half-length ls (rim = diam/2 - ls),
x diode at (0, +pos) bridging the vertical slit, y diode at (+pos, 0). MoM sheet tensor, measured MA4AGP907 at the stack's f0. MODEL_ONLY.
    python sidequest_cell_search.py ro3010|rt5880   -> results/sidequest_hybrid/cell_search_<stack>.json"""
import itertools
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")
import numpy as np

from src.composite_geometry import AxisRect
from src.dispersion import sheet_model as sm
from src.dispersion.leakage_model import leakage_alpha_first_order
from src.frequency import free_space_wavenumber
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import PolyComposite, regular_polygon
from src.ris import diodes_measured as dm

MM = 1e-3
P = 6.4e-3
STACKS = {"rt5880": (Substrate("rt5880", 2.2, 0.0009, 1.575e-3), 10.54e9), "ro3010": (Substrate("ro3010", 11.2, 0.0022, 0.635e-3), 10.05e9)}
ST = {"OFF": "OFF_0V", "ON": "ON_10mA"}


def pieces(ls, pos, wd, foot):
    if not foot:
        return [(-ls, ls)]
    segs = []
    if pos - wd / 2 > -ls:
        segs.append((-ls, pos - wd / 2))
    if pos + wd / 2 < ls:
        segs.append((pos + wd / 2, ls))
    return segs


def cell(stack, sx, sy, diam=6.0, g=0.4, ls=2.8, pos=1.6, wd=1.2, extra_slits=(), sub=None, f0=None):
    """sx, sy in OFF/ON/None (None = slit left open, no diode). extra_slits: static extra vertical+horizontal slits ((x_offset, half_len, width), ...) (loaded-slit variants)."""
    sb, ff = STACKS[stack]
    sb = sub or sb
    disc = regular_polygon(0.0, 0.0, diam / 2 * MM, 48, 0.0)
    ap, loads = [], []
    for axis, state in (("x", sx), ("y", sy)):
        for a0, a1 in pieces(ls, pos, wd, state is not None):
            if axis == "x":
                ap.append(AxisRect(0.0, (a0 + a1) / 2 * MM, g * MM, (a1 - a0) * MM))
            else:
                ap.append(AxisRect((a0 + a1) / 2 * MM, 0.0, (a1 - a0) * MM, g * MM))
        if state is not None:
            if axis == "x":
                loads.append(dm.load(ST[state], 0.0, pos * MM, g * MM, wd * MM, "x"))
            else:
                loads.append(dm.load(ST[state], pos * MM, 0.0, g * MM, wd * MM, "y"))
    for off, half, w in extra_slits:
        for s_ in (-1, 1):
            ap.append(AxisRect(s_ * off * MM, 0.0, w * MM, 2 * half * MM))
            ap.append(AxisRect(0.0, s_ * off * MM, 2 * half * MM, w * MM))
    return UnitCell(RectangularLattice(P, P), PolyComposite(metal=(disc,), apertures=tuple(ap)), sb, Metal(), backing="pec", loads=tuple(loads))


def bloch_loaded(unit_cell, f, mesh_n, kx0, ky0=0.0):
    """Bloch sheet tensor (siemens) of a grounded cell with lumped loads at transverse wavenumber (kx0, ky0): Phase 30 `tensor_sheet.bloch_sheet_tensor` plus the load matrix."""
    from src.basis_functions import MeshSettings, RooftopBasis
    from src.dispersion.tensor_sheet import LayeredGreenFunction
    from src.loads import load_matrix
    from src.mom_matrix import assemble_excitation_vector, assemble_impedance_matrix
    from src.solver import solve

    mesh = MeshSettings(mesh_n, mesh_n, mesh_n)
    green = LayeredGreenFunction(unit_cell.substrate, f, backing="pec")
    basis = RooftopBasis(unit_cell, mesh.cells_x, mesh.cells_y)
    z = assemble_impedance_matrix(basis, green, kx0, ky0, mesh.floquet_order)
    if unit_cell.loads:
        z = z + load_matrix(basis, list(unit_cell.loads), green.omega)
    v = np.stack([assemble_excitation_vector(basis, e, kx0, ky0) for e in np.eye(2)], axis=1)
    cur = solve(z, v)
    j = np.stack([basis.current_coefficients(cur[:, k], kx0, ky0) for k in range(2)], axis=1)
    e_tot = np.eye(2) + np.asarray(green.spectral_dyad(kx0, ky0), dtype=complex) @ j
    return j @ np.linalg.inv(e_tot)


def superpixel(stack, sx, sy, g_c=0.4, g_o=0.2, pos=None, wd=1.2, sub=None):
    """2 x 2 superpixel (6.4 mm): four rectangular sub-patches, central cross gap g_c (x diode bridges the vertical gap at y = +pos, y diode the horizontal gap at x = +pos), outer gap g_o to the neighbours.
    sx, sy in OFF/ON/None (open) / SHORT (metal in the footprint). Sub-patch edges: g_c/2 .. P/2 - g_o/2."""
    from src.composite_geometry import CompositeGeometry

    sb, _ = STACKS[stack]
    sb = sub or sb
    lo, hi = g_c / 2, 3.2 - g_o / 2
    size, ctr = hi - lo, 0.5 * (lo + hi)
    pos = ctr if pos is None else pos
    rects = [AxisRect(sxn * ctr * MM, syn * ctr * MM, size * MM, size * MM) for sxn in (-1, 1) for syn in (-1, 1)]
    loads = []
    for axis, state in (("x", sx), ("y", sy)):
        if state in ("OFF", "ON", "SHORT"):
            if axis == "x":
                rects.append(AxisRect(0.0, pos * MM, g_c * MM, wd * MM))
            else:
                rects.append(AxisRect(pos * MM, 0.0, wd * MM, g_c * MM))
        if state in ("OFF", "ON"):
            loads.append(dm.load(ST[state], 0.0, pos * MM, g_c * MM, wd * MM, "x") if axis == "x" else dm.load(ST[state], pos * MM, 0.0, g_c * MM, wd * MM, "y"))
    return UnitCell(RectangularLattice(P, P), CompositeGeometry(metal=tuple(rects)), sb, Metal(), backing="pec", loads=tuple(loads))


def superpixel_notch(stack, sx, sy, g=0.2, window=0.4, pos=1.6, wd=1.2, sub=None):
    """2 x 2 superpixel with a UNIFORM gap g between all sub-patches (so the sheet stays E1-like) and a local widening of the central gaps to `window` over wd, where the diode sits (x diode at y = +pos in the vertical
    gap, y diode at x = +pos in the horizontal gap). sx, sy as in `superpixel`."""
    from src.composite_geometry import CompositeGeometry

    sb, _ = STACKS[stack]
    sb = sub or sb
    lo, hi = g / 2, 3.2 - g / 2
    notch = (window - g) / 2
    rects = []
    for sxn in (-1, 1):
        for syn in (-1, 1):
            # main block: away from both central gaps by `notch`, then the strips along the gaps (broken by the notches)
            def R(x0, x1, y0, y1):
                rects.append(AxisRect(sxn * (x0 + x1) / 2 * MM, syn * (y0 + y1) / 2 * MM, (x1 - x0) * MM, (y1 - y0) * MM))

            R(lo + notch, hi, lo + notch, hi)
            # strip along the vertical central gap (x in lo..lo+notch): present except in the x-diode window (y in pos-wd/2..pos+wd/2 for the +y half)
            wins_x = [(pos - wd / 2, pos + wd / 2)] if syn == 1 else []
            wins_y = [(pos - wd / 2, pos + wd / 2)] if sxn == 1 else []
            for (a, b) in _complement(lo + notch, hi, wins_x):
                R(lo, lo + notch, a, b)
            for (a, b) in _complement(lo + notch, hi, wins_y):
                R(a, b, lo, lo + notch)
            R(lo, lo + notch, lo, lo + notch)  # corner
    loads = []
    for axis, state in (("x", sx), ("y", sy)):
        if state in ("OFF", "ON", "SHORT"):
            rects.append(AxisRect(0.0, pos * MM, window * MM, wd * MM) if axis == "x" else AxisRect(pos * MM, 0.0, wd * MM, window * MM))
        if state in ("OFF", "ON"):
            loads.append(dm.load(ST[state], 0.0, pos * MM, window * MM, wd * MM, "x") if axis == "x" else dm.load(ST[state], pos * MM, 0.0, wd * MM, window * MM, "y"))
    return UnitCell(RectangularLattice(P, P), CompositeGeometry(metal=tuple(rects)), sb, Metal(), backing="pec", loads=tuple(loads))


def _complement(a, b, wins):
    out, cur = [], a
    for w0, w1 in sorted(wins):
        if w0 > cur:
            out.append((cur, w0))
        cur = max(cur, w1)
    if cur < b:
        out.append((cur, b))
    return out


def alpha_pair(y_off, y_on, sub, f0, theta_deg=25.0):
    """beta/k0, alpha_beam, alpha_abs [Np/m] for a 50/50 pattern of the two sheet admittances (along the propagation axis)."""
    k0 = free_space_wavenumber(f0)
    y_mean = 0.5 * (y_on + y_off)
    beta = sm.find_beta_from_sheet_admittance(1j * y_mean.imag, sub, f0).beta_rad_per_m
    c1 = abs(y_on - y_off) / np.pi
    a_beam = leakage_alpha_first_order(beta, 1j * y_mean.imag, c1, np.radians(theta_deg), 0.0, sub, f0)
    db = 1e-5
    b2 = sm.find_beta_from_sheet_admittance(1j * (y_mean.imag + db), sub, f0).beta_rad_per_m
    return beta / k0, a_beam, (b2 - beta) / db * y_mean.real


def job(p):
    warnings.simplefilter("ignore")
    stack, kw, mesh, eps, label = p
    kw = dict(kw)
    kind = kw.pop("kind", "disc")
    build = superpixel_notch if kind == "notch" else superpixel if kind == "super" else (lambda stack_, sx, sy, sub=None, **k: cell(stack_, sx, sy, sub=sub, **k))
    sb, f0 = STACKS[stack]
    if eps is not None:
        sb = Substrate(sb.name, eps, sb.tan_delta, sb.thickness)
    out = {"stack": stack, "kw": {k: v for k, v in kw.items()}, "mesh": mesh, "eps": eps, "label": label, "Y_mS": {}}
    k0 = free_space_wavenumber(f0)
    beta = (1.3 if stack == "ro3010" else 1.07) * k0
    for it in range(2):  # Bloch tensor at the operating kt, one refinement of beta from the mean sheet
        y = bloch_loaded(build(stack, "OFF", "OFF", sub=sb, **kw), f0, mesh, beta)
        yb = 0.5 * (y[0, 0].imag + y[1, 1].imag)
        try:
            beta = sm.find_beta_from_sheet_admittance(1j * 0.5 * (yb + 0.0), sb, f0).beta_rad_per_m
        except Exception:  # noqa: BLE001
            break
    out["beta_used_over_k0"] = beta / k0
    for sx, sy in (("OFF", "OFF"), ("ON", "OFF"), ("OFF", "ON"), ("ON", "ON")) if label.startswith("full") else (("OFF", "OFF"), ("ON", "OFF")):
        c = build(stack, sx, sy, sub=sb, **kw)
        yx = bloch_loaded(c, f0, mesh, beta, 0.0)  # propagation along x
        yy = bloch_loaded(c, f0, mesh, 0.0, beta) if sy != sx or label.startswith("full") else yx
        y = np.array([[yx[0, 0], yx[0, 1]], [yy[1, 0], yy[1, 1]]])
        out["Y_mS"][f"{sx}/{sy}"] = [[complex(v).real * 1e3, complex(v).imag * 1e3] for v in y.ravel()]
    y = {k: np.array([complex(a, b) for a, b in v]).reshape(2, 2) * 1e-3 for k, v in out["Y_mS"].items()}
    d = y["ON/OFF"][0, 0] - y["OFF/OFF"][0, 0]
    iso = 0.5 * (0.5 * (y["ON/OFF"][0, 0] + y["OFF/OFF"][0, 0]) + 0.5 * (y["ON/OFF"][1, 1] + y["OFF/OFF"][1, 1]))
    try:
        bk, ab, aa = alpha_pair(iso - d / 2, iso + d / 2, sb, f0)
    except Exception as e:  # noqa: BLE001
        bk, ab, aa = None, None, None
    out["analysis"] = {"dY_mS": [d.real * 1e3, d.imag * 1e3], "Y_iso_mS": [iso.real * 1e3, iso.imag * 1e3], "beta_k0": bk, "alpha_beam": ab, "alpha_abs": aa}
    return out


if __name__ == "__main__":
    stack = sys.argv[1]
    grids = {
        "ro3010": [dict(kind="notch", g=0.2, window=0.4, wd=wd) for wd in (1.2, 0.8, 0.4)] + [dict(kind="notch", g=0.2, window=0.4, wd=1.2, pos=0.8)],
        "rt5880": [dict(diam=6.0, ls=2.8), dict(kind="super", g_c=0.4, g_o=0.4), dict(kind="super", g_c=0.4, g_o=0.2), dict(diam=6.0, ls=2.8, pos=0.8), dict(diam=6.0, ls=2.8, extra_slits=((1.0, 2.6, 0.2),)), dict(diam=6.0, ls=2.8, extra_slits=((2.0, 2.0, 0.2),))],
    }
    jobs = [(stack, kw, 64, None, f"v{i}") for i, kw in enumerate(grids[stack])]
    with ProcessPoolExecutor(3) as ex:
        res = list(ex.map(job, jobs))
    for r in res:
        a = r["analysis"]
        print(r["kw"], "dY %.2f%+.2fj Yiso %.2f%+.2fj beta %s alpha_beam %s abs %s" % (*a["dY_mS"], *a["Y_iso_mS"], a["beta_k0"] and round(a["beta_k0"], 3), a["alpha_beam"] and round(a["alpha_beam"], 3), a["alpha_abs"] and round(a["alpha_abs"], 3)), flush=True)
    Path(f"results/sidequest_hybrid/cell_search_{stack}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8", newline="\n")
