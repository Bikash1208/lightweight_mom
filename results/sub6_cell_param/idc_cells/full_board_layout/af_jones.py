# Independent MODEL_ONLY check: radial launch rho_hat*exp(-j beta r)*taper, each cell = half-wave plate R(psi)diag(1,-1)R(-psi),
# far field projected onto explicit RHCP e=(1,-j)/sqrt2 and LHCP e=(1,+j)/sqrt2 (Cartesian transverse, +z propagation).
import numpy as np, math, make_layout_v3 as m
d=m.build(0,True); k0=m.k0
x,y,r,phi,psi0=d['x'],d['y'],d['r'],d['phi'],d['psi']; isb=d['isb']
amp=np.exp(-np.interp(r,[0,75,150],[1.11,1.45,1.83])*r*1e-3)
def field(psi,th,ph):
  c,s=np.cos(2*psi),np.sin(2*psi)
  Ex,Ey=np.cos(phi),np.sin(phi)
  Ox=c*Ex+s*Ey; Oy=s*Ex-c*Ey
  w=amp*np.exp(-1j*m.beta*r)*np.exp(1j*k0*(math.sin(th)*(math.cos(ph)*x+math.sin(ph)*y)))
  return np.sum(w*Ox),np.sum(w*Oy)
eR=np.array([1,-1j])/math.sqrt(2); eL=np.array([1,1j])/math.sqrt(2)
def rep(name,psi):
  out=[]
  for (t,p) in ((30,0),(30,180)):
    ex,ey=field(psi,math.radians(t),math.radians(p))
    out.append((abs(np.conj(eR)@[ex,ey]),abs(np.conj(eL)@[ex,ey])))
  # also overall best over hemisphere for each hand
  best={}
  for hn,e in (('R',eR),('L',eL)):
    bb=(0,)
    for t in range(0,60,4):
      for p in range(0,360,20):
        ex,ey=field(psi,math.radians(t),math.radians(p)); v=abs(np.conj(e)@[ex,ey])
        if v>bb[0]: bb=(v,t,p)
    best[hn]=bb
  print(name,'at(30,0) R/L=%.1f/%.1f'%out[0],' at(30,180) R/L=%.1f/%.1f'%out[1],' global peak R',tuple(round(v,1) for v in best['R']),'L',tuple(round(v,1) for v in best['L']))
rep('LAYOUT (as written)',psi0)
rep('psi negated          ',-psi0)
PhiA=-k0*math.sin(m.TH)*x+m.beta*r; PhiB=k0*math.sin(m.TH)*x+m.beta*r
rep('az sign flipped      ',0.5*(np.where(isb,-PhiB,PhiA)+phi))
rep('az flipped + negated ',-0.5*(np.where(isb,-PhiB,PhiA)+phi))
rep('no azimuth term      ',0.5*np.where(isb,-PhiB,PhiA))
