import json, math, numpy as np
exec(open(r'C:/Users/Bikash/AppData/Local/Temp/3/claude/032c4363339a0143/e8d958c9-c3e5-43ca-aee8-490255dbb25c/scratchpad/an.py').read().split("def rep")[0].replace("import numpy as np,sys",""))
scan = json.load(open(r'C:\Users\Bikash\Desktop\lwmom_labB\bar90_tune_work\s4\p1_c_scan.json'))
bt = np.array([b for b in sorted(float(k) for k in scan) if b <= 1.375]); thR = np.array([scan[str(round(b,3)) if str(round(b,3)) in scan else str(b)][0][0] for b in bt]); thL = np.array([scan[str(round(b,3)) if str(round(b,3)) in scan else str(b)][1][0] for b in bt])
def est_peak(th, ph, pw, near):
    i0 = int(np.argmin(abs(th-near[0]))); j0 = int(np.argmin(abs(ph-near[1]))); n = len(ph)-1 if ph[-1] >= 360 else len(ph)
    best = None
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di; j = (j0+dj) % n
            if 0 <= i < len(th) and (best is None or pw[i, j] > pw[best]): best = (i, j)
    i0, j0 = best; U=[];V=[];Y=[]; dph = ph[1]-ph[0]
    for di in range(-2, 3):
        for dj in range(-2, 3):
            i = i0+di
            if not (0 <= i < len(th)): continue
            j = (j0+dj) % n; t = math.radians(th[i]); p = math.radians(ph[j0]+dj*dph)
            U.append(math.sin(t)*math.cos(p)); V.append(math.sin(t)*math.sin(p)); Y.append(10*math.log10(max(pw[i, j], 1e-300)))
    U=np.array(U);V=np.array(V);Y=np.array(Y); A=np.column_stack([np.ones_like(U),U,V,U*U,U*V,V*V]); c=np.linalg.lstsq(A,Y,rcond=None)[0]
    u,v=np.linalg.solve(np.array([[2*c[3],c[4]],[c[4],2*c[5]]]),-np.array([c[1],c[2]]))
    return math.degrees(math.asin(min(1,math.hypot(u,v)))), math.degrees(math.atan2(v,u))%360
def inv(theta, arr):
    return float(np.interp(theta, arr, bt)) if arr.min()<=theta<=arr.max() else float('nan')
print('scan theta range R %.1f-%.1f L %.1f-%.1f for beta_t %.3f-%.3f' % (thR.min(), thR.max(), thL.min(), thL.max(), bt.min(), bt.max()))
for nm, f in (('BOTH','BOTH.ffd'),('V2 ON-A','ONA.ffd'),('V2 OFF-D','OFFD.ffd'),('B-ALL-SHORT','SHORT.ffd'),('A-only-IDEAL','AONLYI.ffd'),('B-ALL-OFF-IDEAL','BOFFI.ffd')):
    t,ph,et,ep = ffd(r'C:/Users/Bikash/AppData/Local/Temp/3/claude/032c4363339a0143/e8d958c9-c3e5-43ca-aee8-490255dbb25c/scratchpad/'+f)
    pl=abs(et+1j*ep)**2; mi=abs(et-1j*ep)**2
    pR,pL=(pl,mi)
    a = est_peak(t, ph, pR, (31,0)); b = est_peak(t, ph, pL, (30,180))
    print('%-16s RHCP peak (%.2f, %.0f) -> beta_t,A %.3f | LHCP peak (%.2f, %.0f) -> beta_t,B %.3f' % (nm, a[0], a[1], inv(a[0], thR), b[0], b[1], inv(b[0], thL)))
