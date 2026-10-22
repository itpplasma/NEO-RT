"""Independent real sine/cosine oracle across a synthetic sharp response layer."""
import importlib.util
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

import numpy as np

tool=Path(__file__).parents[2]/'tools/write_native_perturbation.py'
spec=importlib.util.spec_from_file_location('write_native',tool)
writer=importlib.util.module_from_spec(spec);spec.loader.exec_module(writer)


def main(binary):
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);eq=root/'synthetic.bc'
        lines=['synthetic']*6
        for rho in [.05,.25,.5,.75,1.]:
            lines.append(f'{rho*rho} -.5 0 0 0 0')
            lines.append('0 0 6 0 0 0 0 0 5 0')
            lines.append(f'1 0 {2*rho} 0 0 {3*rho} 0 {.02*rho/(2*np.pi)} 0 0')
        eq.write_text('\n'.join(lines)+'\n')
        sn=np.array([.05,.4,.400001,.8,.95])
        coeff=np.array([.002+.001j,.004+.002j,-.003+.004j,.003-.002j,.005+.001j])
        np.savez(root/'spectrum.npz',boozer_s=sn,boozer_m=np.array([4]),boozer_total=coeff[:,None])
        writer.write_native(eq,root/'spectrum.npz',root/'bmod_native.dat',3)
        rng=np.random.default_rng(412)
        s=np.r_[rng.uniform(.06,.94,300),.4000001,.4000005,.4000009,.000001,.01,.99,1.1]
        theta=rng.uniform(-np.pi,np.pi,len(s));rho=np.sqrt(s)
        r=600+200*rho*np.cos(theta);z=300*rho*np.sin(theta)
        amp=np.interp(s,sn,coeff.real)+1j*np.interp(s,sn,coeff.imag)
        arg=4*theta-3*.02*rho*np.sin(theta)
        expected=(amp.real*np.cos(arg)-amp.imag*np.sin(arg)) + 1j*(amp.real*np.sin(arg)+amp.imag*np.cos(arg))
        expected[(s<sn[0])|(s>sn[-1])]=0
        with (root/'query.bin').open('wb') as f:
            f.write(struct.pack('=i',len(s)));f.write(r.astype('f8').tobytes());f.write(z.astype('f8').tobytes())
        for threads in (1,8):
            env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
            subprocess.run([binary],cwd=root,env=env,check=True,timeout=60)
            b=np.fromfile(root/'readback.bin',dtype='f8');actual=b[:len(s)]+1j*b[len(s):]
            np.testing.assert_allclose(actual,expected,rtol=2e-7,atol=2e-11)
            np.testing.assert_array_equal(actual[-4:],0)


if __name__=='__main__':
    main(str(Path(sys.argv[1]).resolve()))
