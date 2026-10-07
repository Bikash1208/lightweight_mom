"""Finite N x N square-loop panel on an infinite RT5880 slab, plane-wave (x-pol, normal) incidence.
Cell currents = the periodic MoM solution (local periodicity). Plane-defined S11/S21:
  S11 = <E_refl . x>_A(z=+d) * exp(+j k0 d) / E0,   E_refl = slab reflection + panel scattered field
  S21 = <E_tot  . x>_A(z=-h-d) * exp(+j k0 d) / E0, E_tot = slab transmission + panel scattered field
<.>_A = average over a centred L x L plane.  Spectral integral on a contour lifted above the real kt axis
(exp(+jwt): branch point and lossy slab poles sit just below the axis)."""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))
from fssgeom import unit_cell, H, P
from src.green_function import LayeredGreenFunction

def kz(k, kt):
    v=np.sqrt(k*k-kt*kt+0j); return np.where(v.imag>0,-v,v)

class Panel:
    def __init__(self, npz, N):
        d=np.load(npz); self.f=float(d['f']); self.I=d['coeffs']; self.cx=d['cx']; self.cy=d['cy']; self.dirn=d['dirn']
        self.dx=float(d['dx']); self.dy=float(d['dy']); self.J00=d['J00']; self.s11=d['s11']; self.s21=d['s21']; self.N=N
        self.uc=unit_cell(); self.g=LayeredGreenFunction(self.uc.substrate,self.f,backing='air'); self.k0=self.g.k0; self.ks=self.g.k_substrate
    def slab(self, z):
        """unpatterned-stack reflected (z>0) or transmitted (z<-h) x-field for unit incident E at z=0"""
        e,_=self.g.stack_fields(0.0,0.0,np.array([1.0,0.0]),np.array(z,float),port='top')
        e=np.asarray(e)[...,0]
        return e-np.exp(1j*self.k0*z) if z>0 else e
    def infinite(self, d):
        """Floquet (0,0) check: must equal the solver's S11/S21 (planes at the substrate faces)"""
        out=[]
        for z,sg in ((d,1),(-H-d,-1)):
            es,_=self.g.sheet_fields(0.0,0.0,self.J00,np.array(z,float)); es=np.asarray(es)[...,0]
            out.append((self.slab(z)+es)*np.exp(1j*self.k0*d))
        return out
    def spectral(self, t_max=2500.0, dt=1.5, n_alpha=720, delta=20.0, du=2e-5):
        k0=self.k0; T=2.0*np.real(self.ks)
        t=np.arange(dt/2,t_max,dt); bump=np.where(t<T,np.sin(np.pi*t/T),0.0); dbump=np.where(t<T,np.pi/T*np.cos(np.pi*t/T),0.0)
        kt=t+1j*delta*bump; dk=(1+1j*delta*dbump)*dt
        al=(np.arange(n_alpha)+0.5)*2*np.pi/n_alpha; da=2*np.pi/n_alpha
        a=kz(k0,kt); b=kz(self.ks,kt)
        self.kt,self.w=kt,kt*dk*da/(2*np.pi)**2
        self.al=al; self.a=a
        # per-mode sheet voltages at z=0+ and at z=-h- (unit source), complex-kt safe internals
        mt_up=self.g._sheet_mode_transfer(a,b,np.full(kt.shape,1e-9),'above')
        mt_dn=self.g._sheet_mode_transfer(a,b,np.full(kt.shape,-H-1e-9),'above')
        self.V={'up':{m:mt_up[m][0]*np.exp(1j*a*1e-9) for m in ('te','tm')},
                'dn':{m:mt_dn[m][0]*np.exp(1j*a*1e-9) for m in ('te','tm')}}
        # cell current spectrum per angle via projection binning
        xs=self.dirn==0
        umax=np.max(np.hypot(self.cx,self.cy))+du; edges=np.arange(-umax,umax+du,du); uc=0.5*(edges[1:]+edges[:-1])
        E=np.exp(1j*np.outer(kt,uc))  # (nt, nu)
        Px=np.empty((n_alpha,kt.size),complex); Py=np.empty_like(Px)
        Nn=self.N; idx=np.arange(Nn)-(Nn-1)/2
        AFx=np.empty_like(Px); AFy=np.empty_like(Px)
        for i,alp in enumerate(al):
            ca,sa=np.cos(alp),np.sin(alp); u=self.cx*ca+self.cy*sa
            hx,_=np.histogram(u[xs],edges,weights=self.I[xs].real); hx=hx+1j*np.histogram(u[xs],edges,weights=self.I[xs].imag)[0]
            hy,_=np.histogram(u[~xs],edges,weights=self.I[~xs].real); hy=hy+1j*np.histogram(u[~xs],edges,weights=self.I[~xs].imag)[0]
            kx,ky=kt*ca,kt*sa; sx=np.sinc(kx*self.dx/(2*np.pi)); sy=np.sinc(ky*self.dy/(2*np.pi))
            Px[i]=self.dx*self.dy*sx**2*sy*(E@hx); Py[i]=self.dx*self.dy*sx*sy**2*(E@hy)
            AFx[i]=np.exp(1j*np.outer(kx,idx*P)).sum(1); AFy[i]=np.exp(1j*np.outer(ky,idx*P)).sum(1)
        self.Jx,self.Jy=Px*AFx*AFy,Py*AFx*AFy
        ca,sa=np.cos(al)[:,None],np.sin(al)[:,None]
        self.Jp=self.Jx*ca+self.Jy*sa; self.Jq=-self.Jx*sa+self.Jy*ca; self.ca,self.sa=ca,sa
    def plane(self, L, d):
        kx=self.kt[None,:]*self.ca; ky=self.kt[None,:]*self.sa
        S=np.sinc(kx*L/(2*np.pi))*np.sinc(ky*L/(2*np.pi))
        res=[]
        for side,z in (('up',d),('dn',-H-d)):
            prop=np.exp(-1j*self.a*d)[None,:]
            Ex=(self.V[side]['tm'][None,:]*self.Jp*self.ca - self.V[side]['te'][None,:]*self.Jq*self.sa)*prop
            Ey=(self.V[side]['tm'][None,:]*self.Jp*self.sa + self.V[side]['te'][None,:]*self.Jq*self.ca)*prop
            ex=np.sum(Ex*S*self.w[None,:]); ey=np.sum(Ey*S*self.w[None,:])
            tot=(self.slab(z)+ex)*np.exp(1j*self.k0*d)
            res.append((tot, ey*np.exp(1j*self.k0*d)))
        return res
if __name__=='__main__':
    p=Panel(sys.argv[1],int(sys.argv[2])); lam=2*np.pi/p.k0
    i11,i21=p.infinite(lam/2)
    print('f',p.f/1e9,'solver S11',p.s11[0,0],'inf-check',i11,'| solver S21',p.s21[0,0],'inf-check',i21)
