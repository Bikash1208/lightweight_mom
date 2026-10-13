import sys, json, numpy as np
sys.path.insert(0,'.'); import fixcells
from pathlib import Path
from src.app.routers.hmsa import _resolve_two_state_choice, TwoStateSourceRequest
from src.hologram.two_state_cells import mixed_surface_mode
from src.hologram.coded_synthesis import synthesize_coded_one_bit
from src.hologram.multi_port_synthesis import synthesize_ports
from src.hologram.ports import custom_layout
from src.dispersion.analytic_tm0 import tm0_beta
from src.materials import Substrate
from src.phase8b.object_wave import Beam
OUT=Path(sys.argv[1]); OUT.mkdir(parents=True,exist_ok=True); GAP=float(sys.argv[2])
F=10.54e9; P=6.4; NC=47; S=NC*P; L=5.85; TH=30.0; POL='rhcp'
sub=Substrate('user',2.2,0.0009,1.575e-3)
# patch-plane reflection coefficients: empty (exact, grounded slab) and the 5.85 mm square (mesh-converged estimate)
g_empty=np.exp(1j*np.radians(138.14))*10**(-0.0005/20); g_sq=np.exp(1j*np.radians(-42.0))*10**(-0.033/20)
ts=TwoStateSourceRequest(source='uploaded',uploaded_label_0='empty (bare substrate)',uploaded_value_0_real=g_empty.real,uploaded_value_0_imag=g_empty.imag,
   uploaded_label_1=f'square patch {L} mm',uploaded_value_1_real=g_sq.real,uploaded_value_1_imag=g_sq.imag,
   uploaded_geometry_0={'kind':'empty'},uploaded_geometry_1={'kind':'rectangle','L_mm':L,'W_mm':L}).model_dump()
choice,_=_resolve_two_state_choice(ts,None,F,sub,P*1e-3,beta_spread_limit=float('inf'))
beta=mixed_surface_mode(choice,sub,F,P*1e-3).beta_rad_per_m
g=GAP/2
rects=[(g,g+S,g,g+S),(-g-S,-g,g,g+S),(-g-S,-g,-g-S,-g),(g,g+S,-g-S,-g)]
regions=custom_layout([tuple(v*1e-3 for v in r) for r in rects])
BEAMS=[(TH,45.0),(TH,135.0),(TH,225.0),(TH,315.0)]
pb={r.port_id:[Beam(theta_deg=BEAMS[r.port_id][0],phi_deg=BEAMS[r.port_id][1],weight=1+0j,desired_polarization=POL)] for r in regions}
per_port=lambda ref,beams,feed: synthesize_coded_one_bit(ref,beams,F,choice,P*1e-3,feed_xy_m=feed,surface_wave_betas_rad_per_m=(beta,beta))
multi=synthesize_ports(regions,pb,per_port,P*1e-3,k_swg=tm0_beta(F,2.2,1.575e-3))
cells=[]; ports=[]
for port,(x0,x1,y0,y1) in zip(multi.ports,rects):
    d=port.design; lev=np.asarray(d.level_map); ny,nx=lev.shape
    assert (nx,ny)==(NC,NC), (nx,ny)
    # exact cell centres from integer indices; check they agree with the synthesis grid
    xc=np.array([x0+(j+0.5)*P for j in range(NC)]); yc=np.array([y0+(i+0.5)*P for i in range(NC)])
    dev=max(np.max(np.abs(d.x_m[0,:]*1e3-xc)),np.max(np.abs(d.y_m[:,0]*1e3-yc)))
    fx,fy=(x0+x1)/2,(y0+y1)/2
    keep=0
    for i in range(NC):
        for j in range(NC):
            st=int(lev[i,j])
            if abs(xc[j]-fx)<P*0.6 and abs(yc[i]-fy)<P*0.6:   # 4 cells touching the feed corner -> empty
                if st: keep+=1
                st=0
            cells.append(dict(port=port.port_id+1,i=i,j=j,x=round(xc[j],4),y=round(yc[i],4),state=st))
    ports.append(dict(port=port.port_id+1,rect_mm=[round(v,4) for v in (x0,x1,y0,y1)],feed_mm=[round(fx,4),round(fy,4)],grid_dev_mm=float(dev),
                      ones=int(sum(c['state'] for c in cells if c['port']==port.port_id+1)),keepout_cleared=keep,beam=BEAMS[port.port_id]))
    np.save(OUT/f'port{port.port_id+1}_aperture.npy',np.stack([d.x_m,d.y_m,d.aperture_field.real,d.aperture_field.imag]))
json.dump(dict(frequency_GHz=10.54,cell_mm=P,cells_per_port=NC,port_side_mm=S,gap_mm=GAP,patch_L_mm=L,polarization=POL,beta_rad_per_m=beta,
     states={'0':'empty','1':f'square {L} mm'},ports=ports,cells=cells),open(OUT/'design.json','w'),indent=0)
print(json.dumps(ports,indent=0)); print('beta',beta,'cells',len(cells),'patches',sum(c['state'] for c in cells))
