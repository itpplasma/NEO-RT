#%%
from ctypes import cdll, CDLL, POINTER, RTLD_GLOBAL, c_double, c_int, byref
from _ctypes import dlclose
from numpy import pi, sin, cos, array, zeros, ones, sum, linspace
from matplotlib.pyplot import figure, plot, show, xlabel, ylabel, legend

import sys
sys.path.insert(0, '/home/calbert/code/pytchfort')
from pytchfort import fortran_module # pylint: disable=import-error

# close library if open
try:
    dlclose(libdriftorbit._handle) # pylint: disable=used-before-assignment
except:
    pass
    
libdriftorbit = cdll.LoadLibrary('/home/calbert/build/NEO-RT/libdriftorbit.so')
magfie = fortran_module(libdriftorbit,'do_magfie_mod')
magfie_pert = fortran_module(libdriftorbit,'do_magfie_pert_mod')

#%%
magfie.do_magfie_init()
print('equilibrium field initialised')

magfie_pert.do_magfie_pert_init()
print('perturbation field initialised')

print(magfie.nflux(c_int))
print(magfie_pert.nflux(c_int))


