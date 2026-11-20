# v3 MODEL_ONLY layout. Each cell encodes ONE beam (A = RHCP at (30,0) , B = LHCP at (30,180)), assigned as in
# sub6_key_area_board_model.py: per column, p_B = 0.5*(1-tanh(x/30mm)), exact-count random pick (seed), centre cell (0,0) = feed, excluded.
# Cell set: lattice 25x25, r<=150 mm (same rule as the model). Rotation includes the radial-field azimuth term +phi/2 (sign fixed by the independent Jones check af_jones.py).
# Hand/sign conventions are model assumptions, UNVERIFIED against HFSS far field.
import numpy as np, csv, json, os, math, sys
P=11.6; R=150.0; TH=math.radians(30); F=5.68e9; k0=2*math.pi*F/299792458.0*1e-3; beta=1.10*k0
SEED=int(sys.argv[1]) if len(sys.argv)>1 else 0
AZ=(sys.argv[2]!='noaz') if len(sys.argv)>2 else True
def build(seed,az=True):
  ix,iy=np.meshgrid(np.arange(-12,13),np.arange(-12,13)); ix=ix.ravel(); iy=iy.ravel()
  x,y=ix*P,iy*P; keep=(x*x+y*y<=R*R)&((ix!=0)|(iy!=0)); ix,iy,x,y=ix[keep],iy[keep],x[keep],y[keep]
  r=np.hypot(x,y); phi=np.arctan2(y,x)
  pb=.5*(1-np.tanh(x/30.0)); isb=np.zeros(len(x),bool); rng=np.random.default_rng(seed)
  for c in np.unique(ix):
    ids=np.flatnonzero(ix==c); e=pb[ids].sum(); n=int(math.floor(e)+(rng.random()<e%1))
    isb[rng.choice(ids,size=n,replace=False)]=True
  PhiA=-k0*math.sin(TH)*x+beta*r      # cancel launch phase, aim at phi=0 (model: launch A = exp(-jk sin x))
  PhiB=+k0*math.sin(TH)*x+beta*r
  two_psi=np.where(isb,-PhiB,PhiA)+(phi if az else 0)
  psi=0.5*two_psi
  return dict(ix=ix,iy=iy,x=x,y=y,r=r,phi=phi,isb=isb,psi=psi)
if __name__=='__main__':
  d=build(SEED,AZ); d_=os.path.dirname(os.path.abspath(__file__))
  rows=[]
  for k in range(len(d['x'])):
    r=d['r'][k]; b=bool(d['isb'][k]); x=d['x'][k]
    gate=1 if (r>=5.8 and r<5.8+2*P) and ((x>=0 and not b) or (x<=0 and b)) else 0   # model rule: rho in [5.8mm, 5.8mm+2*P), half-only
    rows.append(dict(i=int(d['ix'][k]),j=int(d['iy'][k]),x_mm=round(x,3),y_mm=round(d['y'][k],3),r_mm=round(r,3),
      rot_deg=round(math.degrees(d['psi'][k])%180,3),beam='B' if b else 'A',cell='ST_R5p45_L8p72_W2p91',gate=gate,
      gate_for='B' if (gate and b) else ('A' if gate else ''),via_block=0))
  with open(d_+'/layout_R150_th30_v3.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
  json.dump(dict(version=3,n_sites=len(rows),n_A=sum(r['beam']=='A' for r in rows),n_B=sum(r['beam']=='B' for r in rows),
   n_gate_A=sum(r['gate_for']=='A' for r in rows),n_gate_B=sum(r['gate_for']=='B' for r in rows),seed=SEED,azimuth_term=AZ,
   rules='lattice 25x25 r<=150, centre (0,0)=feed excluded (model rule); beam per cell by 30 mm tanh dither as in sub6_key_area_board_model.py; gate = first 2 cell-rows (5.8<=rho<5.8+2P, model rule) of cells of that beam, half-only (A x>=0, B x<=0); gate cells KEEP the anchor slot geometry replaced by the gate cell (disc r5.45 plain + optional 0.3 mm via) at run time; axis cells x=0 are in the gate of whichever beam they encode',
   state_rule='A-only: BLOCK (via=1) gate_for=B; B-only: BLOCK gate_for=A; BOTH: all PASS',
   note='MODEL_ONLY; hand/sign unverified'),open(d_+'/manifest_v3.json','w'),indent=1)
  print(len(rows),sum(r['beam']=='B' for r in rows),sum(r['gate'] for r in rows))
