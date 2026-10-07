import sys, json, numpy as np, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
spec=importlib.util.spec_from_file_location('pb','examples/make_paper_board_26ghz.py'); pb=importlib.util.module_from_spec(spec); spec.loader.exec_module(pb)
from src.hologram.two_state_cells import LIBRARY_REFERENCE_OFFSET_M, patch_plane_gamma
from src.app.routers.hmsa import _annular_library_base_config, _run_synthesize_ports_job, HmsaSynthesizePortsRequest, TwoStateSourceRequest
from src.library_sweep import build_annular_slot_library_cells
from src.export.array_export import export_multi_port_array
OUT=Path(sys.argv[1]); OUT.mkdir(parents=True,exist_ok=True); GAP=float(sys.argv[2]); POL=sys.argv[3]; TH=float(sys.argv[4]); BMODE=sys.argv[5] if len(sys.argv)>5 else 'model'
F=10.54e9; P=6.4; b=5.56; R=1.922; r=1.179; LAM=299792458/F*1e3; S=round(5*LAM,1)
pb.F_HZ=F; pb.PERIOD_MM=P; pb.PERIOD_M=P*1e-3
BEAMS=[(TH,45.0),(TH,135.0),(TH,225.0),(TH,315.0)]; pb.BEAMS=BEAMS
ge=patch_plane_gamma(pb._gamma_top_empty(F),F,LIBRARY_REFERENCE_OFFSET_M)
body=dict(period_mm=P,eps_r=2.2,tan_delta=0.0009,thickness_mm=1.575,frequency_hz=F,accuracy='accurate',a_mm=b,b_mm=b,r_mm=r)
c=build_annular_slot_library_cells(_annular_library_base_config(body,placeholder_R_mm=R),b,b,r,[R],F)[0]
gs=patch_plane_gamma(complex(c['S11'][0]['real'],c['S11'][0]['imag']),F,LIBRARY_REFERENCE_OFFSET_M)
ge5=ge
ts=TwoStateSourceRequest(source='uploaded',uploaded_label_0='empty (bare substrate)',uploaded_value_0_real=ge5.real,uploaded_value_0_imag=ge5.imag,
  uploaded_label_1=f'{b} mm square, annular slot r {r} / R {R} mm',uploaded_value_1_real=gs.real,uploaded_value_1_imag=gs.imag,
  uploaded_geometry_0={'kind':'empty'},uploaded_geometry_1={'kind':'annular_slot','a_mm':b,'b_mm':b,'R_mm':R,'r_mm':r}).model_dump()
g=GAP/2
regions=[[g,g+S,g,g+S],[-g-S,-g,g,g+S],[-g-S,-g,-g-S,-g],[g,g+S,-g-S,-g]]
beam=lambda t,p:{'theta_deg':t,'phi_deg':p,'weight_real':1.0,'weight_imag':0.0,'initial_phase_deg':0.0,'desired_polarization':POL,'oam_charge':0}
req=HmsaSynthesizePortsRequest(formulation='coded',frequency_hz=F,period_mm=P,layout={'preset':'custom','regions_mm':regions,'gap_mm':GAP},
  ports=[{'beams':[beam(*BEAMS[i])]} for i in range(4)],port_eps_r=2.2,port_thickness_mm=1.575,port_tan_delta=0.0009,two_state=ts,beta_mode=BMODE).model_dump()
res=_run_synthesize_ports_job(pb._Job(f'coded4_gap{int(GAP)}'),req)
print('synth done',flush=True)
pats=[(p,pb._port_pattern(p,F)[0]) for p in res['ports']]
mets=[pb._port_metrics(p,pat) for p,pat in pats]
pb.figure_patterns(pats, OUT/'port_patterns.png')
pb.figure_layout(res, OUT/'board_layout.png')
feeds=[(p['feed_x_mm'],p['feed_y_mm']) for p in res['ports']]
rep={'frequency_GHz':10.54,'cell_mm':P,'port_side_mm':S,'gap_mm':GAP,'polarization':POL,'beta_mode':BMODE,'beams':BEAMS,
 'cell_pair':{'empty_patch_plane_deg':float(np.degrees(np.angle(ge))),'slotted_patch_plane_deg':float(np.degrees(np.angle(gs))),'slotted':{'b_mm':b,'R_mm':R,'r_mm':r}},
 'two_state':res.get('two_state'),'isolation':res['isolation'],'board_mm':[res['board_width_mm'],res['board_height_mm']],'cells':[len(res['lw_grid_m'][0]),len(res['lw_grid_m'])],
 'feeds_mm':feeds,'ports':mets}
json.dump(rep,open(OUT/'design_report.json','w'),indent=1,default=float); print(json.dumps(rep,indent=1,default=float)[:5000],flush=True)
exp=export_multi_port_array(res['lw_grid_m'],P*1e-3,P*1e-3,pb.SUBSTRATE,pb.METAL,OUT,include_step=True,feeds_mm=feeds,ports_info=[(k+1,x,y) for k,(x,y) in enumerate(feeds)])
print('step',exp.file_size_bytes,flush=True)
