"""Vertical-pin launcher in a grounded slab: spectral (Hankel) transmission-line model, axisymmetric TM.  MODEL_ONLY.
H(k,z) solves (H'/eps)' + (kz^2/eps) H = -k J_s/eps with H'(ground) = 0 and an outgoing wave above; J_s = I(z)/(2 pi); P = 1/2 Re int k dk int dz I (j k H / (w eps)).
Powers: k < k0 = space wave (air), k > k0 = surface wave (+ dielectric loss).  launch = P_sw / (P_sw + P_space).  R_in = 2 P / I0^2.  Thin pin, no probe-radius / reactance model (X not computed).
Check: lambda/4 monopole over PEC -> 36.5 ohm."""
import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.integrate import quad
EPS0, MU0, C0 = 8.8541878128e-12, 4e-7 * math.pi, 299792458.0


def make(f, er, tand, h, t, cur, dz=0.05e-3, ys=0.0, ys_sign=1.0):
    w = 2 * math.pi * f; k0 = w / C0
    ns, na = max(int(round(h / dz)), 4), max(int(round(t / dz)), 1) + 40   # 40 air cells above the pin top
    dzs, dza = h / ns, dz
    z = np.concatenate([-h + (np.arange(ns) + .5) * dzs, (np.arange(na) + .5) * dza])
    d = np.concatenate([np.full(ns, dzs), np.full(na, dza)])
    eps = np.concatenate([np.full(ns, er * (1 - 1j * tand)), np.ones(na)])
    I = np.array([cur(zz) for zz in z]) * (z <= t)
    return dict(w=w, k0=k0, z=z, d=d, eps=eps, I=I, er=er, h=h, ys=ys, isl=ns - 1, ys_sign=ys_sign)


def integrand(M, k):
    w, k0, z, d, eps, I = M["w"], M["k0"], M["z"], M["d"], M["eps"], M["I"]
    n = len(z); kz2 = w * w * MU0 * EPS0 * eps - k * k
    rc = 0.5 * (d[:-1] / eps[:-1] + d[1:] / eps[1:])          # series resistance of the face flux
    if M["ys"]:  # uniform sheet admittance j*ys [S] at z = 0 (between the last slab cell and the first air cell): H jumps by -+ Ys E_rho, E_rho = j flux/(w eps0)
        rc = rc.copy(); rc[M["isl"]] += M["ys_sign"] * (-1j * 1j * M["ys"]) / (w * EPS0)
    g = 1.0 / rc                                               # face conductance (H_{i+1}-H_i) -> flux
    kza = np.sqrt(k0 * k0 - k * k + 0j)
    if kza.imag > 0: kza = -kza
    ab = np.zeros((3, n), complex)
    diag = kz2 / eps * d
    diag[:-1] -= g; diag[1:] -= g
    diag[-1] += -1j * kza                                      # top face flux = -j kza H / 1
    ab[1] = diag; ab[0, 1:] = g; ab[2, :-1] = g
    rhs = -k * (I / (2 * math.pi)) / eps * d
    H = solve_banded((1, 1), ab, rhs)
    p = 0.5 * np.real(np.sum(I * (1j * k * H / (w * EPS0 * eps)) * d)) * k
    return p


def powers(M):
    k0 = M["k0"]; kmax = k0 * math.sqrt(M["er"]) * 1.02
    sp = quad(lambda k: integrand(M, k), 1e-6 * k0, k0 * (1 - 1e-9), limit=200)[0]
    # find the SW peak: scan
    ks = np.linspace(k0 * 1.0001, kmax, 4000); v = np.array([integrand(M, k) for k in ks]); kp = ks[int(np.argmax(v))]
    pts = [kp * (1 + s) for s in (-2e-3, -5e-4, -1e-4, 0, 1e-4, 5e-4, 2e-3)]
    pts = sorted(set([k0 * 1.0001] + [p for p in pts if k0 < p < kmax] + [kmax]))
    sw = sum(quad(lambda k: integrand(M, k), a, b, limit=400, epsabs=0, epsrel=1e-6)[0] for a, b in zip(pts[:-1], pts[1:]))
    return sp, sw, kp / k0


def run(f, er, tand, h, t, cur, i0, ys=0.0, ys_sign=1.0):
    M = make(f, er, tand, h, t, cur, ys=ys, ys_sign=ys_sign)
    sp, sw, kp = powers(M)
    tot = sp + sw
    return {"P_space": sp, "P_sw": sw, "launch": sw / tot, "direct": sp / tot, "R_in_ohm": 2 * tot / i0 ** 2, "beta_over_k0": kp}


if __name__ == "__main__":
    f = 10.5e9; lam = C0 / f
    # check: lambda/4 monopole over a PEC (tiny slab er 1)
    t = lam / 4
    r = run(f, 1.0, 0.0, 0.2e-3, t, lambda z: math.sin(2 * math.pi / lam * (t - z)) if z > 0 else math.sin(2 * math.pi / lam * t), math.sin(2 * math.pi / lam * t) if False else 1.0)
    print("lambda/4 monopole check: R_in %.1f ohm (36.5 expected)  launch %.3f" % (r["R_in_ohm"], r["launch"]))
