# MODEL_ONLY (direction/ratio only): injected edge-reflected waves, circle vs square rim, state 1 and 2. Units: x,y,r mm; k0 rad/mm; attenuation uses metres.
import sys, math
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_labB\results\sub6_cell_param\idc_cells\full_board_layout')
import fwd_beta as F
x, y, r, phi, isb, k0 = F.x, F.y, F.r, F.phi, F.isb, F.k0
print('k0 rad/mm', k0)
BTA, BTB = 1.198, 1.298
bt = np.where(isb, BTB, BTA) * k0            # rad/mm
psi = F.layout_psi(1.10, 1.10)
c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
def alpha_r(rm):  return np.interp(rm * 1e3, [0, 75, 150, 400], [1.11, 1.45, 1.83, 1.83])   # Np/m at radius rm [m]
pB = lambda xmm: 0.5 * (1 - np.tanh(xmm / 30.0))
def aeff(xmm, rm, eps2):  # state 2: B-type share has alpha*eps2
    f = pB(xmm); return alpha_r(rm) * ((1 - f) + f * eps2)
def amp_to(xe, ye, eps2):   # amplitude of the outgoing wave at point (xe,ye) mm, integrated along the ray from the feed
    rr = math.hypot(xe, ye); n = 60; s = np.linspace(0, rr, n); xs = xe * s / max(rr, 1e-9); ys = ye * s / max(rr, 1e-9)
    al = aeff(xs, s * 1e-3, eps2); return math.exp(-np.trapezoid(al, s * 1e-3))
TH = np.radians(np.arange(0, 90.1, 3)); PH = np.radians(np.arange(0, 360, 6))
TT, PP = np.meshgrid(TH, PH, indexing='ij'); DX = np.sin(TT) * np.cos(PP); DY = np.sin(TT) * np.sin(PP)
EF = np.exp(1j * k0 * (DX[..., None] * x + DY[..., None] * y))      # th,ph,cells
def pattern(w, g):
    Ox = c2 * np.cos(g) + s2 * np.sin(g); Oy = s2 * np.cos(g) - c2 * np.sin(g)
    Fx = (EF * (w * Ox)).sum(-1); Fy = (EF * (w * Oy)).sum(-1)
    return {h: np.abs(np.conj(e[0]) * Fx + np.conj(e[1]) * Fy) ** 2 for h, e in (('R', F.eR), ('L', F.eL))}
def amp_dir(eps2): return np.array([amp_to(x[k], y[k], eps2) for k in range(len(x))])
def direct(eps2): return (amp_dir(eps2) * np.exp(-1j * bt * r), phi)
R0 = 150.0; a = 150.0
def circle(G, eps2):
    rim = np.array([amp_to(R0 * math.cos(p), R0 * math.sin(p), eps2) for p in phi])
    back = np.array([math.exp(-np.trapezoid(aeff(x[k] * np.linspace(1, R0 / max(r[k], 1e-9), 30), np.linspace(r[k], R0, 30) * 1e-3, eps2), np.linspace(r[k], R0, 30) * 1e-3)) for k in range(len(x))])
    d = 2 * R0 - r
    return (G * rim * back * np.exp(-1j * bt * d), phi + math.pi)
def square(G, eps2):
    waves = []
    for sx, sy in ((2*a, 0), (-2*a, 0), (0, 2*a), (0, -2*a)):
        dx, dy = x - sx, y - sy; d = np.hypot(dx, dy)
        if sx: tt = (sx - np.sign(sx) * a) / (sx - x); ex = np.sign(sx) * a * np.ones_like(x); ey = y * tt
        else: tt = (sy - np.sign(sy) * a) / (sy - y); ey = np.sign(sy) * a * np.ones_like(y); ex = x * tt
        rim = np.array([amp_to(ex[k], ey[k], eps2) for k in range(len(x))])
        dedge = np.hypot(x - ex, y - ey) * 1e-3
        back = np.exp(-alpha_r(0.15) * dedge)
        waves.append((G * rim * back * np.exp(-1j * bt * d), np.arctan2(dy, dx)))
    return waves
def lvl(P, th_deg, ph_deg):
    return {h: 10 * math.log10(P[h][int(round(th_deg / 3)), int(round(ph_deg / 6)) % len(PH)] + 1e-30) for h in P}
def run(tag, waves, ref):
    P = {h: sum(pattern(w, g)[h] for w, g in waves) for h in ('R', 'L')}
    print(tag, 'peaks R %.2f L %.2f' % tuple(10*math.log10(P[h].max()) for h in 'RL'), '| at (30,0) R %.2f L %.2f | at (30,180) R %.2f L %.2f' % (lvl(P,30,0)['R'], lvl(P,30,0)['L'], lvl(P,30,180)['R'], lvl(P,30,180)['L']))
    i = np.unravel_index(P['R'].argmax(), P['R'].shape); j = np.unravel_index(P['L'].argmax(), P['L'].shape)
    print('    peak R at th %.0f ph %.0f ; peak L at th %.0f ph %.0f' % (math.degrees(TH[i[0]]), math.degrees(PH[i[1]]), math.degrees(TH[j[0]]), math.degrees(PH[j[1]])))
    return P
# known-answer: duplicate the direct wave through the circle code path -> +6.02 dB
d1 = direct(1.0)
P1 = {h: pattern(*d1)[h] for h in 'RL'}
P2 = {h: pattern(d1[0] * 2, d1[1])[h] for h in 'RL'}
print('known answer (G=2 duplicate wave): R %.2f dB L %.2f dB (expect +6.02)' % (10*math.log10(P2['R'].max()/P1['R'].max()), 10*math.log10(P2['L'].max()/P1['L'].max())))
G = 0.3
for st, eps2 in (('state1 (all on-alpha)', 1.0), ('state2 (B off, eps2 0.017)', 0.017), ('state2 (B off, eps2 0.08)', 0.08)):
    print('==', st)
    d = direct(eps2); Pd = run('direct      ', [d], None)
    for nm, w in (('circle G.3', [circle(G, eps2)]), ('square G.3', square(G, eps2))):
        Pr = run('refl only ' + nm, w, None)
        Pt = run('direct+' + nm, [d] + w, None)
        for h in 'RL':
            pk = Pd[h].argmax(); rho = math.sqrt(Pr[h].ravel()[pk] / Pd[h].ravel()[pk]); print('    %s rho at direct peak (refl/direct field) %.3f -> worst case %.2f / %.2f dB' % (h, rho, 20*math.log10(1+rho), 20*math.log10(max(1-rho,1e-3))))
