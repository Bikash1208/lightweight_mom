import sys, time, numpy as np
sys.path.insert(0,'.')
from src.geometry import UnitCell, Geometry
from src.periodic_bc import RectangularLattice
from src.materials import Substrate, Metal
from src.basis_functions import MeshSettings
from src.excitation import ReferencePlanes
from src.sparameters import compute_s_parameters, scattering_from_fields
from src.green_function import LayeredGreenFunction
F=float(sys.argv[1])*1e9; P=6.4e-3; H=1.575e-3; n=int(sys.argv[2]); Ls=[float(x) for x in sys.argv[3].split(',')]
sub=Substrate('RT5880',2.2,0.0009,H); planes=ReferencePlanes(top_air_thickness=0.0,bottom_air_thickness=None,substrate_thickness=H)
ge=complex(scattering_from_fields(LayeredGreenFunction(sub,F,backing='pec'),planes,np.zeros((1,2,2),complex)).s11[0,0])
print(f'f={F/1e9} n={n} empty {np.degrees(np.angle(ge)):.2f} deg',flush=True)
for L in Ls:
    t=time.time()
    uc=UnitCell(RectangularLattice(P,P),Geometry('rectangle',L*1e-3,L*1e-3),sub,Metal(),backing='pec')
    s=compute_s_parameters(uc,F,MeshSettings(n,n,n),planes,warn_grating=False)
    g=complex(s.s11[0,0]); sep=(np.degrees(np.angle(g))-np.degrees(np.angle(ge))+180)%360-180
    print(f' L={L:.2f} mm: {np.degrees(np.angle(g)):8.2f} deg |G| {20*np.log10(abs(g)):.3f} dB  sep {abs(sep):6.1f}  [{time.time()-t:.0f}s]',flush=True)
