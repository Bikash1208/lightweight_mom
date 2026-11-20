import os,sys
os.environ["ANSYSEM_ROOT261"]=r"C:\Program Files\ANSYS Inc\v261\AnsysEM"
from pathlib import Path
t=Path(r"C:\work\s4_work\audit"); os.environ.update(TEMP=str(t),TMP=str(t))
from ansys.aedt.core import Hfss
h=Hfss(project=r"C:\work\s4_work\s1_ALL_PLAIN\BOARD_R150_GATE_V31_ALL_PLAIN.aedt",design="GATE_ALL_PLAIN",solution_type="Modal",version="2026.1",non_graphical=True,new_desktop=True,close_on_exit=False)
try:
    p=h.setups[0].props
    print("FREQ repr:",repr(p["Frequency"]), type(p["Frequency"]))
    print("setups",[s.name for s in h.setups],"solver",p.get("DrivenSolverType"), "keys",[k for k in p if "Freq" in k or "Solver" in k])
    port=h.boundaries_by_type if hasattr(h,"boundaries_by_type") else None
    for b in h.boundaries:
        print("BND",b.name,b.type, getattr(b,"props",{}).get("Objects") if hasattr(b,"props") else None, {k:str(v)[:60] for k,v in b.props.items() if k in("Faces","Objects","DoDeembed","RenormalizeAllTerminals","Impedance","CharImp")})
    print("PORT obj",h.modeler["PORT"].bounding_box, "faces",len(h.modeler["PORT"].faces))
    pin=[n for n in h.modeler.object_names if n.startswith("METAL_PEC")]
    # pin: smallest-radius tall object; ground: widest
    bbs={n:h.modeler[n].bounding_box for n in pin}
    pins=[n for n,b in bbs.items() if b[5]-b[2]>10]; gr=[n for n,b in bbs.items() if b[3]-b[0]>300]
    print("pin",pins,[bbs[n] for n in pins]); print("ground",gr,[bbs[n] for n in gr])
    print("mesh ops",len(h.mesh.meshoperations))
finally:
    h.release_desktop(close_projects=True,close_desktop=True)
