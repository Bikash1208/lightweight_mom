import sys, json, re, shutil, numpy as np
from pathlib import Path
from decimal import Decimal
GAP=int(sys.argv[1]); SRC=Path(f'out/sq300_gap{GAP}'); RUN=Path(f'out/sq4_300_gap{GAP}')
OUT=Path(f'out/final/four_port_1bit_square_300mm_gap{GAP}'); OUT.mkdir(parents=True,exist_ok=True)
d=json.load(open(SRC/'design.json')); rep=json.load(open(RUN/'design_report.json'))
P=Decimal('6.4'); L=Decimal('5.85'); H=Decimal('1.575'); S=47*P; G=Decimal(GAP)/2
W=2*S+2*G; RPIN=Decimal('0.635'); ROUT=Decimal('2.05'); MARK=Decimal('0.3')
# exact decimal geometry (no float arithmetic)
rects=[(G,G+S,G,G+S),(-G-S,-G,G,G+S),(-G-S,-G,-G-S,-G),(G,G+S,-G-S,-G)]
feeds=[((a+b)/2,(c+e)/2) for a,b,c,e in rects]
patches=[]
for c in d['cells']:
    if not c['state']: continue
    x0,_,y0,_=rects[c['port']-1]
    cx=x0+(c['j']*2+1)*P/2; cy=y0+(c['i']*2+1)*P/2
    assert abs(float(cx)-c['x'])<1e-6 and abs(float(cy)-c['y'])<1e-6
    patches.append((c['port'],c['i'],c['j'],cx-L/2,cy-L/2))      # lower-left corner, exact
def s(v):
    t=format(v.normalize(),'f'); return '0' if t in('-0','') else t
# ---------------- layout.csv
with open(OUT/'layout.csv','w') as f:
    f.write('port,row,col,x_min_mm,y_min_mm,x_max_mm,y_max_mm\n')
    for p,i,j,x,y in patches: f.write(','.join([str(p),str(i),str(j),s(x),s(y),s(x+L),s(y+L)])+'\n')
# ---------------- STEP (every vertex placed directly at its exact coordinate; no transforms)
import cadquery as cq
V=lambda x,y,z: cq.Vector(float(x),float(y),float(z))
def rect(x0,y0,x1,y1,z): return cq.Face.makeFromWires(cq.Wire.makePolygon([V(x0,y0,z),V(x1,y0,z),V(x1,y1,z),V(x0,y1,z)],close=True))
hw=W/2
sub=cq.Solid.makeBox(float(W),float(W),float(H),V(-hw,-hw,0))
gnd=cq.Face.makePlane(float(W),float(W),cq.Vector(0,0,0),cq.Vector(0,0,1))
pat=cq.Compound.makeCompound([rect(x,y,x+L,y+L,H) for _,_,_,x,y in patches])
mk=cq.Compound.makeCompound([cq.Face.makeFromWires(cq.Wire.makeCircle(float(MARK),V(fx,fy,0),cq.Vector(0,0,1))) for fx,fy in feeds])
asm=cq.Assembly(name='HMSA_4port')
asm.add(sub,name='Substrate'); asm.add(gnd,name='Ground'); asm.add(pat,name='Patches'); asm.save(str(OUT/'board.step'))
# verify: every coordinate in the STEP is an exact 0.001 mm multiple and is printed without float residue
txt=open(OUT/'board.step').read()
nums=[]
for m in re.finditer(r'CARTESIAN_POINT\([^,]*,\(([^)]*)\)\)',txt): nums+= [t.strip() for t in m.group(1).split(',')]
bad=[t for t in nums if abs(float(t)*1000-round(float(t)*1000))>1e-9]
ugly=[t for t in nums if re.search(r'(9999|0000\d)',t)]
allnum=re.findall(r'[-\d.]+E-\d+',txt); print('tiny exp values',allnum[:5]); print('STEP points',len(nums),'off-grid',len(bad),'ugly-print',len(ugly),ugly[:5])
# ---------------- AEDT IronPython script
pl=',\n'.join('(%s,%s)'%(s(x),s(y)) for _,_,_,x,y in patches)
fl=', '.join('(%s,%s)'%(s(x),s(y)) for x,y in feeds)
tpl=open('aedt_template.py').read()
tpl=tpl.replace('@@PATCHES@@',pl).replace('@@FEEDS@@',fl).replace('@@W@@',s(W)).replace('@@L@@',s(L)).replace('@@H@@',s(H))\
       .replace('@@RPIN@@',s(RPIN)).replace('@@ROUT@@',s(ROUT)).replace('@@GAP@@',str(GAP)).replace('@@N@@',str(len(patches)))
open(OUT/'hfss_build_board.py','w').write(tpl)
for fn in ('port_patterns.png',): shutil.copy(RUN/fn,OUT/fn)
out=dict(frequency_GHz=10.54,substrate='Rogers RT/duroid 5880, er 2.2, tand 0.0009, 1.575 mm, PEC ground',cell_mm=6.4,
  states={'0':'empty (bare substrate)','1':'square copper patch 5.85 x 5.85 mm'},
  pair_patch_plane_phase_deg={'empty':138.14,'square_5.85':-42.0,'separation':180.1},
  board_mm=[float(W),float(W)],port_side_mm=float(S),cells_per_port='47 x 47',gap_mm=GAP,feeds_mm=[[float(a),float(b)] for a,b in feeds],
  feed_keepout='the cell directly under each feed (feed at a cell centre) is forced empty',patches=len(patches),
  beta_rad_per_m=rep['two_state']['beta_used_rad_per_m'],ports=rep['ports'],isolation=rep['isolation'])
json.dump(out,open(OUT/'design_report.json','w'),indent=1,default=float)
print('patches',len(patches),'W',s(W),'feeds',fl)
