import numpy as np,sys
from rdffd import rd
def rh(et,ep): return (et+1j*ep)/np.sqrt(2)
base='C:/Users/Bikash/Desktop/sub6_board_hfss/board_v2/crop/'
F={'open0':base+'res_OPEN_IDEAL/OPEN_IDEAL_1.ffd','open90':base+'res_O90/O90_1.ffd','on0':base+'res_sw_on/sw_on_1.ffd','on90':base+'res_ON90/ON90_1.ffd','a0':base+'res_brA0/brA0_1.ffd','a90':base+'res_brA90/brA90_1.ffd'}
D={}
for k,p in F.items():
    th,ph,et,ep=rd(p); D[k]=(th,ph,et,ep); print(k,ph[0],ph[-1],len(ph))
def reg(k):
    th,ph,et,ep=D[k]; return th,ph,et,ep
# common phi grid: use 0..360 mapping
def tograd(k,phg):
    th,ph,et,ep=D[k]
    out=[]
    for a in (et,ep):
        o=np.empty((len(th),len(phg)),complex)
        for j,p in enumerate(phg):
            q=((p-ph[0])%360)+ph[0]
            idx=np.argmin(abs(ph-q)); o[:,j]=a[:,idx]
        out.append(o)
    return out
phg=np.arange(0,360,2.0); th=D['open0'][0]
G={k:tograd(k,phg) for k in D}
TH,PH=np.meshgrid(th,phg,indexing='ij'); w=np.sin(np.radians(TH))
def orient(a,b):
    et=(G[a][0]-G[b][0])/2; ep=(G[a][1]-G[b][1])/2; return et,ep
for st,(a,b) in {'open':('open0','open90'),'ON':('on0','on90'),'A':('a0','a90')}.items():
    et,ep=orient(a,b); D[st]=(et,ep)
    P=(abs(et)**2+abs(ep)**2)
    print(st,'oriented power total (rel)',(P*w).sum()*2*2*np.pi**2/180**2/1)

eta=376.73; dd=np.radians(2.0)**2
def Wt(x,m): return (abs(x)**2*m*w).sum()*dd/(2*eta)
def cone2(rad,tb,pb):
    c=np.cos(np.radians(TH))*np.cos(np.radians(tb))+np.sin(np.radians(TH))*np.sin(np.radians(tb))*np.cos(np.radians(PH-pb))
    return (np.degrees(np.arccos(np.clip(c,-1,1)))<=rad)*1.0
for (tb,pb,lab) in ((25,0,'design_B(25,0)'),(20,22,'argmax(20,22)')):
    for rad in (10,15):
        m=cone2(rad,tb,pb)
        o=rh(*D['open']);n=rh(*D['ON']);q=rh(*D['A'])
        Po,Pn,Pq=Wt(o,m),Wt(n,m),Wt(q,m)
        print(f"{lab} r{rad}: open {Po:.5f} W  ON {Pn:.5f}  shortA {Pq:.5f} W  shortA/open field {np.sqrt(Pq/Po):.3f}  ON/open {np.sqrt(Pn/Po):.3f}")
m=cone2(10,25,0)
for k in ('open0','open90','a0','a90'):
    print(k,'RHCP W design_B 10deg',round(Wt(rh(G[k][0],G[k][1]),m),5))
