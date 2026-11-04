"""Pass-state matching of the 3-row gate: 1-2 transition rows of the hologram's average sheet on each side (single-mode line model). MODEL_ONLY."""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import numpy as np
from scipy.optimize import minimize

import sidequest_gate_cell as gc

P = gc.P


def abcd(y, n_rows):
    b = gc.beta_complex(y)
    if b is None:
        raise ValueError("no bound mode")
    z = gc.zm(b)
    th = b * n_rows * P
    return np.array([[np.cos(th), 1j * z * np.sin(th)], [1j * np.sin(th) / z, np.cos(th)]]), b


def cascade(yref, sections):
    """sections: list of (Y, n_rows) from the inner (feed) side outwards, embedded in the reference sheet yref on both sides.
    Returns (S11, S21) complex, referenced to the front and back planes of the structure (voltage waves, Zm of the reference)."""
    m = np.eye(2, dtype=complex)
    for y, n in sections:
        a, _ = abcd(y, n)
        m = m @ a
    br = gc.beta_complex(yref)
    zr = gc.zm(br)
    a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    den = a + b / zr + c * zr + d
    return (a + b / zr - c * zr - d) / den, 2.0 / den


def metrics(yref, sections):
    s11, s21 = cascade(yref, sections)
    return {"RL_db": float(-20 * np.log10(max(abs(s11), 1e-12))), "IL_db": float(-20 * np.log10(abs(s21))), "S11": complex(s11), "S21": complex(s21)}


def optimise(yref, y_gate, n_gate=3, n_t=1, b0=(4.0, 6.0), w_il=0.3):
    """Transition sheets (susceptance b_k in mS, k = 1..n_t rows each, outer to inner symmetric); minimise |S11| + w |1 - |S21||."""

    def sections(bs):
        tr = [(1j * b * 1e-3, 1) for b in bs]
        return tr[::-1] + [(y_gate, n_gate)] + tr

    def cost(bs):
        if np.any(np.asarray(bs) < 0.2):
            return 10.0
        try:
            s11, s21 = cascade(yref, sections(bs))
        except Exception:
            return 10.0
        return abs(s11) + w_il * abs(1 - abs(s21))

    best = None
    for start in ([b0[0]] * n_t, [3.0] * n_t, [5.5] * n_t, list(np.linspace(3.5, 6.5, n_t))):
        r = minimize(cost, start, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    bs = [float(b) for b in best.x]
    return bs, metrics(yref, sections(bs)), sections(bs)


if __name__ == "__main__":
    yref = 1j * 7.82e-3
    yg = 0.0009e-3 + 1j * 1.93e-3
    print("no transition", metrics(yref, [(yg, 3)]))
    for nt in (1, 2):
        bs, m, _ = optimise(yref, yg, 3, nt)
        print(nt, [round(b, 2) for b in bs], {k: round(v, 2) for k, v in m.items() if k.endswith("db")})
