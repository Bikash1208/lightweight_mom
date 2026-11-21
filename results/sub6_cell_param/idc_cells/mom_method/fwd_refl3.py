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
def alpha_r(rm):
    al = np.interp(rm * 1e3, [0, 75, 150], [1.11, 1.45, 1.83]); return np.where(rm * 1e3 > R0 + 1e-6, 0.0, al)
pB = lambda xmm: 0.5 * (1 - np.tanh(xmm / 30.0))
def line_int(x0, y0, x1, y1, eps2, n=48):   # integral of alpha_eff along a straight line (mm coordinates), Np
    t = np.linspace(0, 1, n); xs = x0 + (x1 - x0) * t; ys = y0 + (y1 - y0) * t; rm = np.hypot(xs, ys) * 1e-3
    f = pB(xs); al = alpha_r(rm) * ((1 - f) + f * eps2); L = math.hypot(x1 - x0, y1 - y0) * 1e-3
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
states = {'s1': (1.0, np.ones(N)), 's2a (eps .13, eps2 .017)': (0.017, np.where(isb, 0.13, 1.0)), 's2b (eps .28, eps2 .08)': (0.08, np.where(isb, 0.28, 1.0))}
# --- split check: reflected-only circle, A cells only vs B cells only (state 1)
cw = circle_waves(1.0)
for nm, mk in (('A only', (~isb).astype(float)), ('B only', isb.astype(float))):
    A = total(cw, mk); print('SPLIT circle reflected-only,', nm, ': R at (30,0) window %.2f dB, R at (30,180) window %.2f, L at (30,0) %.2f, L at (30,180) %.2f' % (db(win(A,'R',*RW)), db(win(A,'R',*LW)), db(win(A,'L',*RW)), db(win(A,'L',*LW))))
sw = square_waves(1.0)
for nm, mk in (('A only', (~isb).astype(float)), ('B only', isb.astype(float))):
    A = total(sw, mk); print('SPLIT square reflected-only,', nm, ': R(30,0) %.2f R(30,180) %.2f L(30,0) %.2f L(30,180) %.2f' % (db(win(A,'R',*RW)), db(win(A,'R',*LW)), db(win(A,'L',*RW)), db(win(A,'L',*LW))))
res = {}
for sk, (eps2, mk) in states.items():
    d = amp(direct_w(eps2), phi, mk)
    for shape, wf in (('circle', circle_waves), ('square', square_waves)):
        res[(sk, shape)] = (d, total(wf(eps2), mk))
    print(sk, 'direct only: R(30,0) %.2f dB, L(30,180) %.2f dB' % (db(win(d,'R',*RW)), db(win(d,'L',*LW))))
def tot(sk, shape, G, ph):
    d, rf = res[(sk, shape)]; return {h: d[h] + G * np.exp(1j * ph) * rf[h] for h in 'RL'}

cw = circle_waves(1.0)
for leg, nm in ((0, 'inward leg'), (1, 'through-centre leg')):
    for cl, mk in (('A cells', (~isb).astype(float)), ('B cells', isb.astype(float))):
        A = total([cw[leg]], mk)
        print('LEG %s, %s: R(30,0) %.2f  R(30,180) %.2f  L(30,0) %.2f  L(30,180) %.2f dB' % (nm, cl, db(win(A,'R',*RW)), db(win(A,'R',*LW)), db(win(A,'L',*RW)), db(win(A,'L',*LW))))
d1 = res[('s1', 'circle')][0]
GR = (0.02, 0.05, 0.08, 0.1, 0.2, 0.3)
for shape in ('circle', 'square'):
    for sk in ('s2a (eps .13, eps2 .017)', 's2b (eps .28, eps2 .08)'):
        d2 = res[(sk, shape)][0]
        supp0 = db(win(d1,'L',*LW)) - db(win(d2,'L',*LW)); dR0 = db(win(d2,'R',*RW)) - db(win(d1,'R',*RW))
        print('== %s, %s: no-reflection suppression %.2f dB (model, uncalibrated), RHCP change %.2f dB (NOT a prediction)' % (shape, sk, supp0, dR0))
        for G in GR:
            ss = []; dd = []
            for ph in np.radians(np.arange(0, 360, 15)):
                a1 = tot('s1', shape, G, ph); a2 = tot(sk, shape, G, ph)
                ss.append(db(win(a1,'L',*LW)) - db(win(a2,'L',*LW))); dd.append(db(win(a2,'R',*RW)) - db(win(a1,'R',*RW)) - dR0)
            print('   G %.2f: suppression DROP worst %.2f dB ; RHCP change caused by reflection worst |d| %.2f dB' % (G, supp0 - min(ss), max(abs(v) for v in dd)))
