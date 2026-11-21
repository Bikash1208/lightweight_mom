import numpy as np,math
exec(open('an.py').read().split("def rep")[0])
def ang(t1,p1,t2,p2):
    a1,a2=np.radians(t1),np.radians(t2);dp=np.radians(p1-p2)
    return np.degrees(np.arccos(np.clip(np.cos(a1)*np.cos(a2)+np.sin(a1)*np.sin(a2)*np.cos(dp),-1,1)))
def run(name,p,GR,GL):
    t,ph,et,ep=ffd(p)
    if abs(ph[-1]-ph[0]-360)<1e-6: et=et[:,:-1];ep=ep[:,:-1];ph=ph[:-1]
    dth=math.radians(t[1]-t[0]);dph=math.radians(ph[1]-ph[0]);w=np.sin(np.radians(t))[:,None]
    UR=abs(et+1j*ep)**2/2;UL=abs(et-1j*ep)**2/2
    P=((UR+UL)*w).sum()*dth*dph
    T,PH=np.meshgrid(t,ph,indexing='ij')
    print(name)
    for hand,U,Gk,(bt,bp) in (('R',UR,GR,(30,0)),('L',UL,GL,(30,180))):
        D=4*math.pi*U.max()/P; scale=10**(Gk/10)/(4*math.pi*U.max()/P)  # eff factor
        Gm=4*math.pi*U/P*scale   # gain, normalised to accepted power
        pk=np.unravel_index(U.argmax(),U.shape)
        a=ang(T,PH,t[pk[0]],ph[pk[1]])
        tot=(Gm*w).sum()*dth*dph/(4*math.pi)
        near=(Gm*w*(a<8)).sum()*dth*dph/(4*math.pi)
        bt2=ang(T,PH,30,0);bt3=ang(T,PH,30,180)
        out=(Gm*w*((bt2>=8)&(bt3>=8)&(T<=70))).sum()*dth*dph/(4*math.pi)
        low=(Gm*w*(T>70)).sum()*dth*dph/(4*math.pi)
        # HPBW area: pixels within 3 dB of peak
        hp=(U>=U.max()/2).sum()
        print(f" {hand}: peak th {t[pk[0]]:.0f} ph {ph[pk[1]]:.0f}  D {10*math.log10(D):.2f} dBi  G {Gk}  radiated frac of accepted {tot:.3f}  within8deg {near:.3f}  mid(outside both beams,th<=70) {out:.3f}  th>70 {low:.3f}  3dB-area pixels {hp}")
run('BOTH','BOTH.ffd',9.033,8.996)
run('ON-A','ONA.ffd',8.45,6.36)
run('OFF-D','OFFD.ffd',8.64,1.66)
