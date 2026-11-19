# MODEL_ONLY hologram layout: anchor ST_R5p45_L8p72_W2p91 everywhere, R<=150 mm (cell centre r<=150-5.8), pitch 11.6
# Rotation psi = 0.5*arg(gR + conj(gL)); RHCP uses e^{+j2psi}, LHCP e^{-j2psi} (geometric-phase convention).
# SIGN CONVENTION IS UNVERIFIED against HFSS far field: if RHCP/LHCP peaks come out swapped, negate psi.
import numpy as np, csv, json, os
P=11.6; R=150.0; TH=np.deg2rad(30); F=5.68e9; C=299792458.0
k0=2*np.pi*F/C*1e-3            # rad/mm
alpha=lambda r: 1.11+ (1.83-1.11)*r/150   # Np/m profile, informational
beta=1.10*k0                   # anchor beta/k0 ~1.10
n=int(R//P)+2; rows=[]
for i in range(-n,n+1):
  for j in range(-n,n+1):
    x,y=i*P,j*P; r=np.hypot(x,y)
    if r>R-5.8+1e-9: continue            # whole 5.8 mm half-pitch cell inside R
    ref=beta*r
    gR=np.exp(1j*(k0*np.sin(TH)*x - ref))      # beam at phi=0
    gL=np.exp(1j*(-k0*np.sin(TH)*x - ref))     # beam at phi=180
    psi=0.5*np.angle(gR+np.conj(gL))
    half='A' if x>0 else ('B' if x<0 else 'A')
    # gate: first two cell rows around the feed on each half (|x| centre <= 2 pitches from axis? ring of radius<=2 cell rings)
    gate = 1 if (r<=2.5*P and r>0.1) else 0   # study-1 2-row gate; centre cell (feed) excluded
    rows.append(dict(i=i,j=j,x_mm=round(x,3),y_mm=round(y,3),r_mm=round(r,3),rot_deg=round(np.degrees(psi)%180,3),
      cell='ST_R5p45_L8p72_W2p91' if not gate else 'GATE_DISC_R5p45',half=half,gate=gate,via_block=0))
d=os.path.dirname(__file__)
with open(d+'/layout_R150_th30.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
json.dump(dict(n_sites=len(rows),n_gate=sum(r['gate'] for r in rows),theta=30,f_GHz=5.68,pitch=P,R=R,
  note='MODEL_ONLY layout; rotation sign unverified; dither zone and alpha profile NOT yet applied per cell; via_block set per state by the run (A-only: half B gate via=1, B-only: half A gate via=1)'),open(d+'/manifest.json','w'),indent=1)
print(len(rows),sum(r['gate'] for r in rows))
