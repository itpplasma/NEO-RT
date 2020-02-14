# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 09:41:35 2016

@author: Christopher Albert
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

n = 101
s = np.arange(n)*1.0/(n-1)

infile = "../RUN/Mtprofile_Hamiltonian.dat"
outfile = "../RUN/Mtprofile.dat"

data = np.loadtxt(infile)

Mt  = InterpolatedUnivariateSpline(data[:,0], data[:,1])
vth = InterpolatedUnivariateSpline(data[:,0], data[:,2])

out = np.zeros([n,3])
out[:,0] = s
out[:,1] = Mt(s)
out[:,2] = vth(s)

plt.plot(s, out[:,1])

np.savetxt(outfile, out)
