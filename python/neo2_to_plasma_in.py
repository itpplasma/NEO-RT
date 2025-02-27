# %%
import sys

sys.path.append("/var/tmp/ert/src/NEO-2/PythonScripts")

# %%
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# %%

# basedir = "/temp/ert/CONDOR/Neo2/DEMO/2023-07-04_mars_perturbation_100kAt_n1_positive/"
# basedir = "/temp/buchholz/Scans/vmec/DEMO/mars_perturbation_100kAt_n_minus1/neo-2_phi/"

# filename = "neo2_multispecies_out.h5"

# with h5py.File(os.path.join(basedir, filename)) as f:
#     print(f.keys())

# %%
infile = "/temp/grassl_g/GPEC_NEO2_AUG_30835/RUN_stor_lag6/aug_30835.in"

with h5py.File(infile) as f:
    s_neo2 = np.array(f['boozer_s'])
    n_neo2 = np.array(f['n_prof'])
    T_neo2 = np.array(f['T_prof'])

ks0 = 36

n_interp = interp1d(s_neo2[ks0:], n_neo2[:,ks0:],
    kind='cubic', fill_value='extrapolate')
T_interp = interp1d(s_neo2[ks0:], T_neo2[:,ks0:],
    kind='cubic', fill_value='extrapolate')
ns = 1000

s_out = np.arange(0.0, 1.0, 1.0/ns) + 0.5/ns
n_out = n_interp(s_out)
T_out = T_interp(s_out)

plt.figure()
plt.plot(s_out, n_out.T)

plt.figure()
plt.plot(s_out, T_out.T)

# %%
eV = 6.2415091e+11
with open('plasma.in', 'w') as fout:
    fout.write(' % N am1 am2 Z1 Z2' + '\n')
    fout.write('{ns}         0.2D+01         0.2D+01         '
               '0.1D+01         0.1D+01' + '\n')
    fout.write(' % s ni_1[cm^-3] ni_2[cm^-3] '
               ' Ti_1[eV] Ti_2[eV] Te[eV]' + '\n')
    for ks in range(ns):
        fout.write(
            f'{s_out[ks]:.8e} {n_out[1,ks]:.8e} {0.0:.8e} '
            f'{eV*T_out[1,ks]:.8e} {0.0:.8e} {eV*T_out[0,ks]:.8e}' + '\n')

# %%

vth =

with open('profile.in', 'w') as fout:
