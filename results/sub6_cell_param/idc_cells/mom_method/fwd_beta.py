# MODEL_ONLY per-hand effective-beta forward model on the v3 layout (af_jones physics, S4 model_core style): relative peak power vs the physical per-hand beta_t.
import sys, math
import numpy as np
sys.path.insert(0, r'C:\Users\Bikash\Desktop\lwmom_labB\results\sub6_cell_param\idc_cells\full_board_layout')
import make_layout_v3 as m
d = m.build(0, True); k0 = m.k0
x, y, r, phi, isb = d['x'], d['y'], d['r'], d['phi'], d['isb']
AMP = np.exp(-np.interp(r, [0, 75, 150], [1.11, 1.45, 1.83]) * r * 1e-3)
eR = np.array([1, -1j]) / math.sqrt(2); eL = np.array([1, 1j]) / math.sqrt(2)
TH = np.radians(np.arange(16, 50, 0.5)); 
def layout_psi(bdA, bdB, thA=30.0, thB=30.0):
    PhiA = -k0 * math.sin(math.radians(thA)) * x + bdA * k0 * r
    PhiB = +k0 * math.sin(math.radians(thB)) * x + bdB * k0 * r
    return 0.5 * (np.where(isb, -PhiB, PhiA) + phi)
def peaks(psi, btA, btB, amp=None):
    a = AMP if amp is None else amp
    bt = np.where(isb, btB, btA) * k0
    c, s = np.cos(2 * psi), np.sin(2 * psi)
    Ox = c * np.cos(phi) + s * np.sin(phi); Oy = s * np.cos(phi) - c * np.sin(phi); w = a * np.exp(-1j * bt * r)
    out = {}
    for hn, e, phis in (('R', eR, np.radians(np.arange(-12, 12.1, 2))), ('L', eL, np.radians(np.arange(168, 192.1, 2)))):
        TT, PP = np.meshgrid(TH, phis, indexing='ij'); DX = np.sin(TT) * np.cos(PP); DY = np.sin(TT) * np.sin(PP)
        Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
        for k in range(len(x)):
            ph_ = np.exp(1j * k0 * (DX * x[k] + DY * y[k])) * w[k]; Fx += ph_ * Ox[k]; Fy += ph_ * Oy[k]
        P = np.abs(np.conj(e[0]) * Fx + np.conj(e[1]) * Fy) ** 2
        i, j = np.unravel_index(P.argmax(), P.shape); out[hn] = (10 * math.log10(P.max()), math.degrees(TH[i]), math.degrees(phis[j]))
    return out
if __name__ == '__main__':
    psi0 = layout_psi(1.10 * k0 / k0 * 1.0, 1.10)
    ref = peaks(psi0, 1.198, 1.298)
    print('BOTH-like (beta_t A 1.198, B 1.298, layout 1.10): R %.2f dB at th %.1f | L %.2f dB at th %.1f' % (ref['R'][0], ref['R'][1], ref['L'][0], ref['L'][1]))
    on = peaks(psi0, 1.087, 0.88)
    print('ON-A-like (A 1.087, B 0.88): R %.2f dB (delta %.2f, measured -0.58) at th %.1f (HFSS 30.1) | L %.2f dB (delta %.2f, measured -2.64) at th %.1f (HFSS 25.1)' % (on['R'][0], on['R'][0]-ref['R'][0], on['R'][1], on['L'][0], on['L'][0]-ref['L'][0], on['L'][1]))
    off = peaks(psi0, 1.199, 1.298)
    print('OFF-like A 1.199: R %.2f dB (delta %.2f vs BOTH-like; measured OFF-D -0.37 / datasheet OFF -0.39)' % (off['R'][0], off['R'][0]-ref['R'][0]))
