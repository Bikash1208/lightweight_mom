import pandas as pd, numpy as np
c=299792458.0
sub=pd.read_csv('SUB.csv'); rlc=pd.read_csv('RLC.csv'); ref=pd.read_csv('bare_grounded_slab_reference.csv'); m=pd.read_csv('mom_reference.csv')
cols=lambda d,i,j: d[f're(S(FloquetPort1:{i},FloquetPort1:{j})) []']+1j*d[f'im(S(FloquetPort1:{i},FloquetPort1:{j})) []']
for d in (sub,rlc): d.columns=[x.strip() for x in d.columns]
# fit plane offset from slab state 1
s=sub[sub['state []']==1]
rows=[]
for r in s.itertuples(index=False):
    th,f=r[0],r[2]
    for mode,lab in ((1,'1 (TE)'),(2,'2 (TM)')):
        h=getattr(r,'_%d'%0) if False else None
    pass
k=lambda f: 2*np.pi*f*1e9/c
A=[];B=[]
for _,r in s.iterrows():
    th,f=r['scan_theta [deg]'],r['Freq [GHz]']
    for mode,lab in ((1,'1 (TE)'),(2,'2 (TM)')):
        h=r[f're(S(FloquetPort1:{mode},FloquetPort1:{mode})) []']+1j*r[f'im(S(FloquetPort1:{mode},FloquetPort1:{mode})) []']
        a=ref[(ref.f_ghz==f)&(ref.theta_deg==th)&(ref.floquet_mode==lab)].iloc[0]
        dphi=np.angle(np.exp(1j*np.radians(a.S11_phase_deg))/h)  # rad
        A.append(2*k(f)*np.cos(np.radians(th))); B.append(dphi)
A=np.array(A);B=np.array(B); d=np.sum(A*B)/np.sum(A*A)
res=np.degrees(B-A*d); print('plane offset d = %.4f mm, residual max %.2f deg'%(d*1e3,np.abs(res).max()))
def corr(df,f,th): return np.exp(1j*2*k(f)*np.cos(np.radians(th))*d)
def hfss_S(df):
    out={}
    for _,r in df.iterrows():
        th,st,f=r['scan_theta [deg]'],int(r['state []']),r['Freq [GHz]']
        g=lambda i,j: (r[f're(S(FloquetPort1:{i},FloquetPort1:{j})) []']+1j*r[f'im(S(FloquetPort1:{i},FloquetPort1:{j})) []'])*corr(df,f,th)
        out[(st,f,th)]={'S_TE_TE':g(1,1),'S_TM_TM':g(2,2),'S_TE_TM':g(1,2),'S_TM_TE':g(2,1)}
    return out
def circ(S,sx=1):
    TT,EE,ME,EM=S['S_TM_TM'],S['S_TE_TE'],sx*S['S_TM_TE'],sx*S['S_TE_TM']
    return {'r_RR':0.5*(TT-EE+1j*ME+1j*EM),'r_LL':0.5*(TT-EE-1j*ME-1j*EM),'r_LR':0.5*(TT+EE+1j*ME-1j*EM),'r_RL':0.5*(TT+EE-1j*ME+1j*EM)}
H=hfss_S(rlc)
mm={}
for _,r in m.iterrows(): mm.setdefault((int(r.state),r.f_ghz,r.theta_deg),{})[r.quantity]=10**(r.mom_mag_dB/20)*np.exp(1j*np.radians(r.mom_phase_deg))
db=lambda x:20*np.log10(abs(x)); ph=lambda x:np.degrees(np.angle(x))
for sx in (1,-1):
    e=[]
    for key in H:
        if key in mm:
            ch=circ(H[key],sx); e.append(abs(ch['r_RR']-mm[key]['r_RR']))
    print('cross sign',sx,'mean |dr_RR|',np.mean(e))

def table(f,th):
    print(f'\n== f={f} GHz theta={th} ==')
    print('state | HFSS |rRR| dB  phase | MoM |rRR| dB phase | HFSS opp dB | MoM opp dB')
    hp=[];mp=[]
    for st in (1,2,3,4):
        h=circ(H[(st,f,th)]); o=mm[(st,f,th)]
        hp.append(ph(h['r_RR'])); mp.append(ph(o['r_RR']))
        print(st,'| %6.2f %8.1f | %6.2f %8.1f | %6.1f | %6.1f'%(db(h['r_RR']),ph(h['r_RR']),db(o['r_RR']),ph(o['r_RR']),max(db(h['r_LR']),db(h['r_RL'])),max(db(o['r_LR']),db(o['r_RL']))))
    st=lambda p:[((p[i+1]-p[i])%360) for i in range(3)]+[((p[0]-p[3])%360)]
    print('steps HFSS',np.round(st(hp),1),' MoM',np.round(st(mp),1))
for th in (0.0,30.0):
    for f in (9.5,10.0,10.5,11.0): table(f,th)
# bandwidth summary: freq where all |rRR|>=-1.5, opp<=-15, steps 90+-15
print()
for th in (0.0,30.0):
  for src in ('HFSS','MoM'):
    ok=[]
    for f in sorted({k[1] for k in H}):
        rr=[(circ(H[(s,f,th)]) if src=='HFSS' else mm[(s,f,th)]) for s in (1,2,3,4)]
        p=[ph(x['r_RR']) for x in rr]; steps=[((p[(i+1)%4]-p[i])%360) for i in range(4)]
        dev=max(abs(x-90) for x in steps); amp=min(db(x['r_RR']) for x in rr); opp=max(max(db(x['r_LR']),db(x['r_RL'])) for x in rr)
        ok.append((f,round(dev,1),round(amp,2),round(opp,1)))
    print(th,src,ok)
