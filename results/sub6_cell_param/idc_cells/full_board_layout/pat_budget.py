# MODEL_ONLY aperture-pattern budget of the v3 layout: directivity and beam share of the aperture-radiated power
import sys, math, numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_labB\results\sub6_cell_param\idc_cells\full_board_layout')
import make_layout_v3 as m
d = m.build(0, True); k0 = m.k0
x, y, r, phi, psi = d['x'], d['y'], d['r'], d['phi'], d['psi']
alpha_r = np.interp(r, [0, 75, 150], [1.11, 1.45, 1.83])
th = np.radians(np.arange(0, 91, 2)); ph = np.radians(np.arange(0, 360, 2))
c, s = np.cos(2 * psi), np.sin(2 * psi); Ex, Ey = np.cos(phi), np.sin(phi)
Ox = c * Ex + s * Ey; Oy = s * Ex - c * Ey
eR = np.array([1, -1j]) / math.sqrt(2); eL = np.array([1, 1j]) / math.sqrt(2)
def pattern(amp, beta):
    PR = np.zeros((len(th), len(ph))); PL = np.zeros_like(PR)
    for i, t in enumerate(th):
        for j, p in enumerate(ph):
            w = amp * np.exp(-1j * beta * r) * np.exp(1j * k0 * math.sin(t) * (math.cos(p) * x + math.sin(p) * y))
            ex, ey = np.sum(w * Ox), np.sum(w * Oy)
            PR[i, j] = abs(np.conj(eR) @ [ex, ey]) ** 2; PL[i, j] = abs(np.conj(eL) @ [ex, ey]) ** 2
    return PR, PL
def report(name, amp, beta):
    PR, PL = pattern(amp, beta)
    w = np.sin(th)[:, None] * np.radians(2) * np.radians(2)
    P = ((PR + PL) * w).sum()
    out = []
    for hn, Pp in (('R', PR), ('L', PL)):
        i, j = np.unravel_index(Pp.argmax(), Pp.shape)
        D = 4 * math.pi * Pp.max() / P
        tg, pg = np.meshgrid(th, ph, indexing='ij')
        cosd = np.cos(tg) * np.cos(th[i]) + np.sin(tg) * np.sin(th[i]) * np.cos(pg - ph[j])
        near = np.degrees(np.arccos(np.clip(cosd, -1, 1))) <= 8
        out.append('%s peak (%d,%d) D %.1f dB, share of aperture power within 8 deg %.3f' % (hn, round(math.degrees(th[i])), round(math.degrees(ph[j])), 10 * math.log10(D), (Pp * w * near).sum() / P))
    print(name, '|', ' | '.join(out))
amp_a = np.exp(-alpha_r * r * 1e-3)
amp_b = np.maximum(r * 1e-3, 0.0058) ** -0.5 * np.exp(-alpha_r * r * 1e-3)
for bn, beta in (('beta as built', m.beta), ('beta 1.16 k0', 1.16 * k0)):
    report('exp(-a r) taper, ' + bn, amp_a, beta)
    report('rho^-1/2 exp(-a r), ' + bn, amp_b, beta)
