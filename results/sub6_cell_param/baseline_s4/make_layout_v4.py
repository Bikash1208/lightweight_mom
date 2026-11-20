"""v4 layout generator (S4 baseline track). Same rules as idc_cells/full_board_layout/make_layout_v3.py, but beta is a PER-HAND input.
CSV columns unchanged: i,j,x_mm,y_mm,r_mm,rot_deg,beam,cell,gate,gate_for,via_block.
Usage: python make_layout_v4.py OUT.csv --beta-a 1.16 --beta-b 1.29 [--seed 0] [--no-gate] [--theta-a 30 --theta-b 30]
Default (beta 1.10/1.10, seed 0, gates on) reproduces v3 / the v31 CSV exactly (asserted by p2_make_layouts.py).
Beta values are in units of k0. B-hand beta must not be frozen for a final layout until the diode-cell result is in."""
import argparse, csv, math, sys
import numpy as np
P = 11.6; R = 150.0; F = 5.68e9; K0 = 2*math.pi*F/299792458.0*1e-3   # rad/mm

def build(seed=0, beta_a=1.10, beta_b=1.10, theta_a=30.0, theta_b=30.0, az=True, no_gate=False, width=30.0):
    ix, iy = np.meshgrid(np.arange(-12, 13), np.arange(-12, 13)); ix = ix.ravel(); iy = iy.ravel()
    x, y = ix*P, iy*P; keep = (x*x+y*y <= R*R) & ((ix != 0) | (iy != 0)); ix, iy, x, y = ix[keep], iy[keep], x[keep], y[keep]
    r = np.hypot(x, y); phi = np.arctan2(y, x)
    pb = .5*(1-np.tanh(x/width)); isb = np.zeros(len(x), bool); rng = np.random.default_rng(seed)
    for c in np.unique(ix):
        ids = np.flatnonzero(ix == c); e = pb[ids].sum(); n = int(math.floor(e)+(rng.random() < e % 1))
        isb[rng.choice(ids, size=n, replace=False)] = True
    PhiA = -K0*math.sin(math.radians(theta_a))*x + beta_a*K0*r
    PhiB = +K0*math.sin(math.radians(theta_b))*x + beta_b*K0*r
    psi = 0.5*(np.where(isb, -PhiB, PhiA) + (phi if az else 0))
    rows = []
    for k in range(len(x)):
        rr = r[k]; b = bool(isb[k]); xx = x[k]
        gate = 0 if no_gate else (1 if (rr >= 5.8 and rr < 5.8+2*P) and ((xx >= 0 and not b) or (xx <= 0 and b)) else 0)
        rows.append(dict(i=int(ix[k]), j=int(iy[k]), x_mm=round(xx, 3), y_mm=round(y[k], 3), r_mm=round(rr, 3),
                         rot_deg=round(math.degrees(psi[k]) % 180, 3), beam='B' if b else 'A', cell='ST_R5p45_L8p72_W2p91',
                         gate=gate, gate_for='B' if (gate and b) else ('A' if gate else ''), via_block=0))
    return rows

def write(rows, path):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

if __name__ == '__main__':
    a = argparse.ArgumentParser(); a.add_argument('out'); a.add_argument('--beta-a', type=float, default=1.10); a.add_argument('--beta-b', type=float, default=1.10)
    a.add_argument('--theta-a', type=float, default=30.0); a.add_argument('--theta-b', type=float, default=30.0)
    a.add_argument('--seed', type=int, default=0); a.add_argument('--no-gate', action='store_true'); a.add_argument('--width', type=float, default=30.0)
    s = a.parse_args(); rows = build(s.seed, s.beta_a, s.beta_b, s.theta_a, s.theta_b, True, s.no_gate, s.width); write(rows, s.out)
    print(len(rows), 'cells; B', sum(r['beam'] == 'B' for r in rows), 'gate', sum(r['gate'] for r in rows))
