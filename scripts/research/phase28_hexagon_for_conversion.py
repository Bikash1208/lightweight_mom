"""Phase 28: smallest hexagon (6.4 mm pitch, 30/225 beam) whose march converts >= 50 % and >= 70 % of the surface wave, for the 10.0 GHz table, the table capped at the item-2 alpha_cap values, and the table capped at G's alpha."""
import sys,json
sys.path.insert(0,'scripts/research')
import phase28_small_board as sb
from src.hologram import coded_feasibility as cf
from src.hologram.tensor_leaky_cp import SupercellTable
T=SupercellTable.from_files(sb.TABLE_GLOB)
res={}
for n,c in {'table':None,'cap3.66':3.661,'cap2.09':2.089,'G_0.45':0.4457}.items():
    tab=T if c is None else cf.capped_table(T,c)
    def eta(R):
        x,y,d,m=cf.design_cells(tab,cf.hexagon_m(R*1e-3),6.4e-3,30.0,225.0); return len(x),float(m['eta_conv'])
    for tgt in (0.5,0.7):
        lo,hi=100.0,100.0
        while eta(hi)[1]<tgt and hi<1500: lo=hi; hi*=1.4
        if eta(hi)[1]<tgt: res[(n,tgt)]=None; print(n,tgt,'not reached by',hi,eta(hi),flush=True); continue
        for _ in range(9):
            mid=(lo+hi)/2
            if eta(mid)[1]<tgt: lo=mid
            else: hi=mid
        nc,e=eta(hi); res[(n,tgt)]=(round(hi,1),nc,8*nc,round(e,3)); print(n,tgt,res[(n,tgt)],flush=True)
json.dump({f'{k[0]}|{k[1]}':v for k,v in res.items()},open('results/phase28/f10/g_refine/hexagon_for_conversion_50_70.json','w'))
