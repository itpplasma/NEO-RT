#!/usr/bin/env python3
"""Write an opt-in POTATO native-spectrum file from axisymmetric BC geometry.

Geometry: not-a-knot cubics in sqrt(s), exact angular Fourier sums. Perturbation:
piecewise linear on every original radial node, without mode filtering. Output
uses cm, G and the single-n cylindrical complex amplitude, with no conjugation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.interpolate import CubicSpline


def read_bc(path):
    surfaces,blocks=[],[]
    for line in Path(path).read_text().splitlines()[6:]:
        try:
            row=[float(v.replace('D','E')) for v in line.split()]
        except ValueError:
            continue
        if len(row)==6:
            surfaces.append(row);blocks.append([])
        elif len(row)==10 and blocks:
            blocks[-1].append(row)
    return np.asarray(surfaces),np.asarray(blocks)


def support_disk(cs,gm,surface,inner):
    """Conservative support disk from a Fourier derivative bound and winding.

    |d(R,Z)/dtheta| <= sum |m| (|C_m|_2+|S_m|_2). A nearest
    angular sample is at most h/2 away. The inner disk is used only when
    its center has winding number one and the bound excludes the origin.
    """
    c=cs(np.sqrt(surface))[:,:2]
    center=c[gm==0].real.sum(axis=0)
    th=np.arange(1024)*2*np.pi/1024
    xy=np.real(np.exp(1j*np.outer(th,gm))@c)-center
    radii=np.linalg.norm(xy,axis=1)
    lipschitz=np.sum(abs(gm)*(np.linalg.norm(c.real,axis=1)+np.linalg.norm(c.imag,axis=1)))
    gap=lipschitz*np.pi/len(th)+1e-8
    if inner:
        lower=max(0.,float(min(radii)-gap))
        angles=np.angle(xy[:,0]+1j*xy[:,1])
        winding=np.sum(np.angle(np.exp(1j*(np.roll(angles,-1)-angles))))/(2*np.pi)
        if abs(abs(winding)-1)>1e-8:
            lower=0.
        radius=lower
    else:
        radius=float(max(radii)+gap)
    return np.r_[center,radius**2]


def write_native(equilibrium,components,output,n_tor):
    es,ec=read_bc(equilibrium)
    if np.any(ec[:,:,1]!=0) or n_tor==0:
        raise ValueError('need axisymmetric geometry and nonzero perturbation n')
    rho=np.sqrt(es[:,0]);gm=ec[0,:,0].astype('i4')
    coeff=(ec[:,:,2:8:2]-1j*ec[:,:,3:8:2])*[100,100,2*np.pi]
    cs=CubicSpline(rho,coeff,axis=0)
    with np.load(components) as d:
        s,m,a=d['boozer_s'],d['boozer_m'].astype('i4'),d['boozer_total']
    if np.any(np.diff(s)<=0) or a.shape!=(len(s),len(m)):
        raise ValueError('invalid native spectrum')
    if s[0]<rho[0]**2 or s[-1]>rho[-1]**2:
        raise ValueError('spectrum support extends beyond geometry')
    disks=np.r_[support_disk(cs,gm,s[0],True),support_disk(cs,gm,s[-1],False)]
    rh,th=np.meshgrid(np.linspace(rho[0],rho[-1],33),np.arange(128)*2*np.pi/128,indexing='ij')
    geometry=np.real(np.einsum('nmk,nm->nk',cs(rh.ravel()),np.exp(1j*np.outer(th.ravel(),gm))))
    seed=np.array([geometry[:,0],geometry[:,1],rh.ravel(),th.ravel()])
    with Path(output).open('wb') as f:
        f.write(b'POTNAT02')
        f.write(struct.pack('=6i',len(rho),len(gm),len(s),len(m),n_tor,seed.shape[1]))
        for x,dtype in ((rho,'f8'),(gm,'i4'),(cs.c.real,'f8'),(-cs.c.imag,'f8'),
                        (s,'f8'),(m,'i4'),(a*1e4,'c16'),(seed,'f8')):
            f.write(np.asarray(x,dtype=dtype).tobytes(order='F'))
        f.write(disks.astype('f8').tobytes())
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    meta=dict(format='POTNAT02',n_tor=n_tor,native_radial_nodes=len(s),modes=m.tolist(),
        radial_interpolation='all original nodes, piecewise linear; no smoothing',
        geometry='original BC not-a-knot cubic in rho=sqrt(s); exact Fourier sum',
        phase='exp(-i*n*(phi_geom-phi_B)); no conjugation or factor two',
        units='cm, G; geometry toroidal shift radians',outside_support='exact zero',
        zero_support_disks_cm_squared=disks.tolist(),
        sources={str(p):sha(p) for p in (equilibrium,components)},output_sha256=sha(output))
    Path(str(output)+'.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('equilibrium','components','output'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--n-tor',type=int,required=True)
    a=p.parse_args()
    print(json.dumps(write_native(a.equilibrium,a.components,a.output,a.n_tor)))
