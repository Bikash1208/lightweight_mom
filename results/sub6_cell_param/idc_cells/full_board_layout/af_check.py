# MODEL_ONLY array-factor check: sites e^{-j beta r}-launched wave with hologram phases; RHCP uses beam-A cells? Field of cell:
# A-cell radiates RHCP e^{j(2psi+phi_term)} B-cell radiates LHCP e^{-j(2psi+phi)}. Gate removal = zero the cells beyond-gate? Implemented as in model:
# one-beam state: blocked half's cells (all cells of that beam outside gate region keep amplitude via transmission T) -> here T power.
import numpy as np, math, sys
sys.path.insert(0,'.'); import make_layout_v3 as m
k0=m.k0
def peak(d,amp,th,ph,hand,az):
  # field of RHCP from A cells, LHCP from B cells. phase = launch(-beta r) + hand*(2psi + phi_az)
  x,y,r,phi,psi=d['x'],d['y'],d['r'],d['phi'],d['psi']
  isb=d['isb']; sel=(~isb) if hand=='R' else isb
  s=1 if hand=='R' else -1
  ph_cell=-m.beta*r + s*(2*psi+(phi if az else 0))
  ux=math.sin(th)*math.cos(ph); uy=math.sin(th)*math.sin(ph)
  return abs(np.sum(amp*sel*np.exp(1j*(ph_cell+k0*(ux*x+uy*y)))))
for az in (True,False):
  d=m.build(0,az)
  amp=np.ones(len(d['x']))
  # amplitude taper by alpha
  amp=np.exp(-np.interp(d['r'],[0,75,150],[1.11,1.45,1.83])*d['r']*1e-3)
  def best(hand,amp):
    return max((peak(d,amp,math.radians(t),math.radians(p),hand,az),t,p) for t in range(20,41,2) for p in (0,180))
  bR=best('R',amp); bL=best('L',amp)
  print('az_in_layout',az,'BOTH  RHCP peak',round(bR[0],2),'at',bR[1:], ' LHCP',round(bL[0],2),'at',bL[1:])
  # one-beam: block half's gate (T=0.01 amplitude) -> everything beyond gate in that half loses power ~T
  for st in ('A-only','B-only'):
    blk='B' if st=='A-only' else 'A'
    a2=amp.copy()
    # gate blocks surface wave into half: cells of that half beyond the 2-row gate (x of blocked half) get amplitude T
    half=(d['x']<0) if blk=='B' else (d['x']>0)
    # perfect-gate: surface wave to cells of that half blocked beyond r>2P
    a2[half&(d['r']>2*m.P)]*=0.0316   # -30 dB
    bR2=best('R',a2); bL2=best('L',a2)
    print('  ',st,'RHCP',round(20*math.log10(bR2[0]/bR[0]),1),'dB LHCP',round(20*math.log10(bL2[0]/bL[0]),1),'dB vs BOTH')
