"""Phase 2: generate fix layouts, assert identity to the v31 CSV (beam/gate/gate_for etc.), pre-solve af_jones-physics check on the NEW CSV."""
import csv, math, sys, os
sys.path.insert(0, '.'); sys.path.insert(0, r"C:\work\lwmom_s4\results\sub6_cell_param\idc_cells\full_board_layout")
import numpy as np
import make_layout_v4 as g
V31 = r"C:\work\lwmom_s4_ref\results\sub6_meander_model\full_board_gate\v31\layout_R150_th30_v31.csv"
ref = list(csv.DictReader(open(V31)))
os.makedirs('layouts_s4', exist_ok=True)
# 0. default reproduces v31 exactly (all columns)
d0 = g.build(0, 1.10, 1.10)
assert len(d0) == len(ref)
for a, b in zip(d0, ref):
    for k in a:
        assert str(a[k]) == str(b[k]) or abs(float(a[k])-float(b[k])) < 1e-9, (k, a[k], b[k])
print("default generator reproduces v31 CSV exactly (516 rows, all columns)")
def same_except_rot(new, tag, gate_off=False):
    for a, b in zip(new, ref):
        for k in ('i', 'j', 'x_mm', 'y_mm', 'r_mm', 'beam', 'cell', 'via_block'):
            assert str(a[k]) == str(b[k]), (tag, k)
        if not gate_off:
            for k in ('gate', 'gate_for'): assert str(a[k]) == str(b[k]), (tag, k)
    print(tag, "beam/gate/gate_for identical to v31" if not gate_off else "beam identical; gate=0 everywhere", "; rot_deg changed in", sum(a['rot_deg'] != b['rot_deg'] for a, b in zip(new, ref)), "cells")
n2 = g.build(0, 1.16, 1.29); same_except_rot(n2, "S2 beta 1.16/1.29"); g.write(n2, 'layouts_s4/layout_R150_th30_s4_beta1p16_1p29.csv')
n3 = g.build(0, 1.16, 1.29, no_gate=True); same_except_rot(n3, "S3 beta 1.16/1.29 + gate=0", True); g.write(n3, 'layouts_s4/layout_R150_th30_s4_beta1p16_1p29_nogate.csv')
assert [a['rot_deg'] for a in n2] == [a['rot_deg'] for a in n3]

# af_jones physics (radial launch rho_hat*exp(-j beta_t r)*taper, half-wave plates R(psi)diag(1,-1)R(-psi), projections on RHCP/LHCP) driven by the CSV rot_deg
from model_core import X, Y, R_, PHI, AMP, k0, est_peak, TH2, PH2, TT, DX, DY, eR, eL
from model_core import X as MX
import json
def far_csv(rows, bt_a, bt_b):
    key = {(int(r['i']), int(r['j'])): r for r in rows}
    from model_core import IX, IY
    psi = np.array([math.radians(float(key[(int(i), int(j))]['rot_deg'])) for i, j in zip(IX, IY)])
    isb = np.array([key[(int(i), int(j))]['beam'] == 'B' for i, j in zip(IX, IY)])
    bt = np.where(isb, bt_b, bt_a)*k0; c, s = np.cos(2*psi), np.sin(2*psi)
    Ox = c*np.cos(PHI)+s*np.sin(PHI); Oy = s*np.cos(PHI)-c*np.sin(PHI); w = AMP*np.exp(-1j*bt*R_)
    Fx = np.zeros(TT.shape, complex); Fy = np.zeros(TT.shape, complex)
    for k in range(len(X)):
        e = np.exp(1j*k0*(DX*X[k]+DY*Y[k]))*w[k]; Fx += e*Ox[k]; Fy += e*Oy[k]
    R = np.abs(np.conj(eR[0])*Fx+np.conj(eR[1])*Fy)**2; L = np.abs(np.conj(eL[0])*Fx+np.conj(eL[1])*Fy)**2
    return R, L
out = {}
for tag, rows, bts in (("v31 (as solved) beta_t=1.10 [unmodified af_jones physics]", ref, (1.10, 1.10)),
                       ("v31 (as solved) beta_t=fit 1.158/1.286", ref, (1.158, 1.286)),
                       ("NEW 1.16/1.29 beta_t=1.10 [unmodified af_jones physics]", n2, (1.10, 1.10)),
                       ("NEW 1.16/1.29 beta_t=fit 1.158/1.286", n2, (1.158, 1.286))):
    R, L = far_csv(rows, *bts)
    pa = est_peak(TH2, PH2, R, (30, 0)); pb = est_peak(TH2, PH2, L, (30, 180))
    # global peak per hand (hand check): grid argmax over theta>6
    gr = np.unravel_index(int(np.argmax(np.where(TH2[:, None] > 6, R, 0))), R.shape); gl = np.unravel_index(int(np.argmax(np.where(TH2[:, None] > 6, L, 0))), L.shape)
    print(f"{tag:62s} RHCP peak th {pa[0]:5.2f} ph {pa[1]:6.1f} (grid max at {TH2[gr[0]]:.0f},{PH2[gr[1]]:.0f}) | LHCP th {pb[0]:5.2f} ph {pb[1]:6.1f} (grid max at {TH2[gl[0]]:.0f},{PH2[gl[1]]:.0f}) | LHCP@(30,0) vs RHCP@(30,0): {10*math.log10(L[15,0]/R[15,0]):6.1f} dB")
    out[tag] = dict(rh=[pa[0], pa[1]], lh=[pb[0], pb[1]])
json.dump(out, open('p2_af_check.json', 'w'), indent=1)
