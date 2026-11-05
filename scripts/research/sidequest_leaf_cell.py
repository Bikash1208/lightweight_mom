"""Leaf-cell (Zhao, Pan, Shen, IEEE TAP 74(1) 2026, doi 10.1109/TAP.2025.3621516) in our MoM: Y_sheet at normal incidence -> TM0 beta -> X_s.  MODEL_ONLY.
The paper's Table I / Fig. 2 are not readable from the text: the leaf shape is RECONSTRUCTED (lens of two circles r1, r2, tip-to-tip length L, tips at the 0.5 mm gap, 0.6 mm pads). Not a calibration.
    python sidequest_leaf_cell.py <case.json-ish via argv: tag>  ->  results/sidequest_leaf/<tag>.json"""
import json, math, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent)); warnings.simplefilter("ignore")
from src.composite_geometry import AxisRect
from src.constants import ETA0
from src.dispersion import sheet_model as sm
from src.frequency import angular_frequency, free_space_wavenumber
from src.geometry import UnitCell
from src.hologram.coded_cell_eval import sheet_tensor
from src.loads import LumpedLoad
from src.materials import Metal, Substrate
from src.periodic_bc import RectangularLattice
from src.polygon_geometry import Poly
MM = 1e-3


def lens(r1, r2, L, x0, rot):
    """leaf: intersection of two discs (r1 below the axis, r2 above), tips at x0 and x0+L on the axis, rotated by rot [rad] about the origin; returns Poly (m)."""
    d1, d2 = math.sqrt(r1 ** 2 - L ** 2 / 4), math.sqrt(r2 ** 2 - L ** 2 / 4)
    xc = x0 + L / 2
    t1 = math.asin(L / 2 / r1); t2 = math.asin(L / 2 / r2)
    up = [(xc + r1 * math.sin(a), -d1 + r1 * math.cos(a)) for a in np.linspace(-t1, t1, 25)]       # upper arc of disc 1
    lo = [(xc + r2 * math.sin(a), d2 - r2 * math.cos(a)) for a in np.linspace(t2, -t2, 25)]       # lower arc of disc 2
    pts = up + lo[1:-1]
    c, s = math.cos(rot), math.sin(rot)
    return Poly(tuple(((c * x - s * y) * MM, (s * x + c * y) * MM) for x, y in pts))


def cell(a, L, rot_deg, state, sub, h_gap=0.5, pad=0.6, r1=6.4, r2=6.0, f=6.9e9, diode="smp"):
    rot = math.radians(rot_deg)
    metal = []
    for sgn in (+1, -1):
        rr = rot + (0 if sgn > 0 else math.pi)
        metal.append(lens(r1, r2, L, h_gap / 2, rr))
        pc = sgn * (h_gap / 2 + pad / 2)
        if rot_deg == 0:
            metal.append(AxisRect(pc * MM, 0.0, pad * MM, pad * MM))
    from src.polygon_geometry import oriented_rect
    if diode == "ma" and state == "OFF":  # HFSS-confirmed MA4AGP907 model: OFF = open gap, no load
        return UnitCell(RectangularLattice(a * MM, a * MM), __import__("src.polygon_geometry", fromlist=["PolyComposite"]).PolyComposite(metal=tuple(metal)), sub, Metal(), backing="pec", loads=())
    metal.append(AxisRect(0.0, 0.0, h_gap * MM, pad * MM) if rot_deg == 0 else oriented_rect(0.0, 0.0, h_gap * MM, pad * MM, rot))  # load footprint = conductor carrying the load
    w = pad * MM
    zf = (lambda ff: complex(2.44, 0.0)) if (diode == "ma") else (lambda ff: complex(2.0, 2 * math.pi * ff * 0.45e-9)) if state == "ON" else (lambda ff: complex(0.0, 2 * math.pi * ff * 0.45e-9 - 1 / (2 * math.pi * ff * 0.23e-12)))
    load = LumpedLoad(0.0, 0.0, h_gap * MM, w, rot_deg if rot_deg else "x", z_of_f=zf)
    return UnitCell(RectangularLattice(a * MM, a * MM), type("PC", (), {})() if False else __import__("src.polygon_geometry", fromlist=["PolyComposite"]).PolyComposite(metal=tuple(metal)), sub, Metal(), backing="pec", loads=(load,))


def xs_from_y(y, sub, f, H):
    """TM0 beta with the sheet Y and (a) open air above, (b) PEC lid H above; X_s = eta alpha/k0 [tanh(alpha H)]."""
    om, k0 = angular_frequency(f), free_space_wavenumber(f)
    eps_c = sub.eps_complex
    m0 = sm.find_beta_from_sheet_admittance(y, sub, f)
    def fl(b):
        al = math.sqrt(b * b - k0 * k0)
        yup = 1j * om * 8.8541878128e-12 / (al * math.tanh(al * H))
        return (y + yup + sm._y_slab(b, om, eps_c, k0, sub.thickness)).imag
    bs = np.linspace(k0 * 1.0005, m0.beta_rad_per_m * 1.6, 4000); v = np.array([fl(b) for b in bs]); i = np.where(np.diff(np.sign(v)) != 0)[0]
    bl = float(brentq(fl, bs[i[0]], bs[i[0] + 1])) if len(i) else float("nan")
    xl = ETA0 * math.sqrt((bl / k0) ** 2 - 1) * math.tanh(math.sqrt(bl * bl - k0 * k0) * H) if bl == bl else float("nan")
    db = 1e-3 * k0
    fi = lambda b: (y + 1j * om * 8.8541878128e-12 / math.sqrt(b * b - k0 * k0) + sm._y_slab(b, om, eps_c, k0, sub.thickness))
    bb = m0.beta_rad_per_m; dI = (fi(bb + db).imag - fi(bb - db).imag) / (2 * db); alpha = abs(fi(bb).real / dI)
    return {"alpha_np_per_m": alpha, "beta_over_k0_open": m0.beta_rad_per_m / k0, "X_open_ohm": m0.x_s_ohm, "beta_over_k0_lid": bl / k0, "X_lid_ohm": xl}


def job(c):
    warnings.simplefilter("ignore")
    f = c["f"]; sub = Substrate("s", c["er"], c["tand"], c["h"] * MM)
    uc = cell(c["a"], c["L"], c["rot"], c["state"], sub, c.get("gap", 0.5), c.get("pad", 0.6), c.get("r1", 6.4), c.get("r2", 6.0), f, c.get("diode", "smp"))
    y = sheet_tensor(uc, f, c["mesh"])
    yv = complex(y[0, 0]) if c["rot"] == 0 else complex(0.5 * (y[0, 0] + y[1, 1]))
    H = 3e8 / f / 4
    r = xs_from_y(yv, sub, f, H)
    return {**c, "Yxx_mS": [y[0, 0].real * 1e3, y[0, 0].imag * 1e3], "Yyy_mS": [y[1, 1].real * 1e3, y[1, 1].imag * 1e3], **r}


if __name__ == "__main__":
    cases = json.loads(sys.argv[2]); tag = sys.argv[1]
    with ProcessPoolExecutor(min(len(cases), 32)) as ex:
        res = list(ex.map(job, cases))
    for r in res:
        print(r["state"], "rot", r["rot"], "L", r["L"], "mesh", r["mesh"], "Y %.2f%+.2fj mS" % tuple(r["Yxx_mS"]), "X_open %.1f X_lid %.1f ohm alpha %.2f Np/m" % (r["X_open_ohm"], r["X_lid_ohm"], r["alpha_np_per_m"]), flush=True)
    Path(f"results/sidequest_leaf/{tag}.json").write_text(json.dumps(res, indent=1), encoding="utf-8", newline="\n")
