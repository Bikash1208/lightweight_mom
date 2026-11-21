# MODEL_ONLY (direction/ratio only). Edge-reflection injection, circle (inward + through-centre legs) vs square (4 first-order images), state 1 / state 2 (B cells x eps).
import sys, math
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_labB\results\sub6_cell_param\idc_cells\full_board_layout')
import fwd_beta as F
x, y, r, phi, isb, k0 = F.x, F.y, F.r, F.phi, F.isb, F.k0
BTA, BTB = 1.198, 1.298
bt = np.where(isb, BTB, BTA) * k0
psi = F.layout_psi(1.10, 1.10); c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
R0 = 150.0; a = 150.0; N = len(x)
RK = float(np.sort(r[isb])[149]); MODE = [1]; PD = [0.0]
print('B cells', int(isb.sum()), 'key radius (150th B cell) %.1f mm' % RK)
def alpha_r(rm):
    al = np.interp(rm * 1e3, [0, 75, 150], [1.11, 1.45, 1.83]); return np.where(rm * 1e3 > R0 + 1e-6, 0.0, al)
pB = lambda xmm: 0.5 * (1 - np.tanh(xmm / 30.0))
def line_int(x0, y0, x1, y1, eps2, n=48):   # integral of alpha_eff along a straight line (mm coordinates), Np
    t = np.linspace(0, 1, n); xs = x0 + (x1 - x0) * t; ys = y0 + (y1 - y0) * t; rm = np.hypot(xs, ys) * 1e-3
    f = pB(xs); fb = np.where(rm * 1e3 < RK, eps2, PD[0]); al = alpha_r(rm) * ((1 - f) + f * fb); L = math.hypot(x1 - x0, y1 - y0) * 1e-3
    return np.trapezoid(al, t) * L
TH = np.radians(np.arange(0, 90.1, 3)); PH = np.radians(np.arange(0, 360, 6))
TT, PP = np.meshgrid(TH, PH, indexing='ij'); DX = np.sin(TT) * np.cos(PP); DY = np.sin(TT) * np.sin(PP)
EF = np.exp(1j * k0 * (DX[..., None] * x + DY[..., None] * y))
def amp(w, g, mask=None):
    if mask is not None: w = w * mask
    Ox = c2 * np.cos(g) + s2 * np.sin(g); Oy = s2 * np.cos(g) - c2 * np.sin(g)
    Fx = (EF * (w * Ox)).sum(-1); Fy = (EF * (w * Oy)).sum(-1)
    return {h: np.conj(e[0]) * Fx + np.conj(e[1]) * Fy for h, e in (('R', F.eR), ('L', F.eL))}
def direct_w(eps2):
    return np.array([math.exp(-line_int(0, 0, x[k], y[k], eps2)) for k in range(N)]) * np.exp(-1j * bt * r)
def circle_waves(eps2):
    ex = R0 * np.cos(phi); ey = R0 * np.sin(phi)
    inward = np.array([math.exp(-line_int(0, 0, ex[k], ey[k], eps2) - line_int(ex[k], ey[k], x[k], y[k], eps2)) for k in range(N)]) * np.exp(-1j * bt * (2 * R0 - r))
    thru = np.array([math.exp(-line_int(0, 0, -ex[k], -ey[k], eps2) - line_int(-ex[k], -ey[k], x[k], y[k], eps2)) for k in range(N)]) * np.exp(-1j * bt * (2 * R0 + r) - 1j * math.pi / 2)
    return [(inward, phi + math.pi), (thru, phi)]
def square_waves(eps2):
    ws = []
    for sx, sy in ((2*a, 0), (-2*a, 0), (0, 2*a), (0, -2*a)):
        dx, dy = x - sx, y - sy; d = np.hypot(dx, dy)
        if sx: tt = (sx - np.sign(sx) * a) / (sx - x); ex = np.sign(sx) * a * np.ones(N); ey = y * tt
        else: tt = (sy - np.sign(sy) * a) / (sy - y); ey = np.sign(sy) * a * np.ones(N); ex = x * tt
        w = np.array([math.exp(-line_int(0, 0, ex[k], ey[k], eps2) - line_int(ex[k], ey[k], x[k], y[k], eps2)) for k in range(N)]) * np.exp(-1j * bt * d)
        ws.append((w, np.arctan2(dy, dx)))
    return ws
def total(waves, mask=None):
    out = {'R': 0, 'L': 0}
    for w, g in waves:
        A = amp(w, g, mask)
        for h in 'RL': out[h] = out[h] + A[h]
    return out
def win(A, h, th0, th1, ph0, ph1):
    P = np.abs(A[h]) ** 2; ti = [i for i, t in enumerate(np.degrees(TH)) if th0 <= t <= th1]
    pi = [j for j, p in enumerate(np.degrees(PH)) if (ph0 <= p <= ph1) or (ph0 > ph1 and (p >= ph0 or p <= ph1))]
    return P[np.ix_(ti, pi)].max()
db = lambda p: 10 * math.log10(p + 1e-30)
RW = (24, 39, 354, 6); LW = (24, 39, 174, 186)
keyc = isb & (r <= RK)
states = {'s1': (1.0, np.where(isb, np.where(keyc, 1.0, 0.0), 1.0), 1), 's2a (eps .13, eps2 .017)': (0.017, np.where(isb, np.where(keyc, 0.13, 0.0), 1.0), 2), 's2b (eps .28, eps2 .08)': (0.08, np.where(isb, np.where(keyc, 0.28, 0.0), 1.0), 2)}
def amp_to_(xe, ye, eps2):
    return math.exp(-line_int(0, 0, xe, ye, eps2))

PD[0] = 0.0
d1 = None; out = {}
for sk in ('s1', 's2a (eps .13, eps2 .017)'):
    eps2, mk, md = states[sk]; MODE[0] = md
    d = amp(direct_w(eps2), phi, mk); cw = circle_waves(eps2)
    print(sk, 'B-rim amp (-150,0) %.3f' % amp_to_(-150.0, 0.0, eps2))
    full = total(cw, mk)
    mkA0 = mk * np.where(isb, 1.0, 0.0)           # remove A cells from the INWARD leg only
    inw = total([cw[0]], mkA0); thr = total([cw[1]], mk)
    noA = {h: inw[h] + thr[h] for h in 'RL'}
    out[sk] = (d, full, noA)
d1, f1, n1 = out['s1']; d2, f2, n2 = out['s2a (eps .13, eps2 .017)']
supp0 = db(win(d1,'L',*LW)) - db(win(d2,'L',*LW))
for G in (0.02, 0.05, 0.1):
    for nm, (r1, r2) in (('full circle', (f1, f2)), ('A-inward leg removed', (n1, n2))):
        ss = []
        for ph in np.radians(np.arange(0, 360, 15)):
            a1 = {h: d1[h] + G*np.exp(1j*ph)*r1[h] for h in 'RL'}; a2 = {h: d2[h] + G*np.exp(1j*ph)*r2[h] for h in 'RL'}
            ss.append(db(win(a1,'L',*LW)) - db(win(a2,'L',*LW)))
        print('G %.2f %-22s suppression DROP worst %.2f dB (no-refl %.2f)' % (G, nm, supp0 - min(ss), supp0))
