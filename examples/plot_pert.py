#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 18:07:06 2017

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi

data = np.loadtxt('orbit/bnoverb0_max_TFripple_RMP90_pgf.dat')

sint = np.append(data[:,2]**2,1.0)
epsint = np.append(data[:,3]/2.0e-3, data[-1,3]/2.0e-3)

relpert = spi.InterpolatedUnivariateSpline(sint, epsint)

stest = np.linspace(0.0,1.0,100)


fin = open('orbit/tok-synch2-n3.bc','r')
fout = open('orbit/in_file_pert','w')

ns = 63
nm = 37

for k in range(4):
    fout.write(fin.readline())

fout.write(fin.readline())

inline = fin.readline()
globdata = np.fromstring(inline,sep=' ')
print(globdata)
fout.write(inline)

fsdata = np.zeros([ns,6])
for ks in range(ns):
    fout.write(fin.readline())
    fout.write(fin.readline())
    inline = fin.readline()
    fsdata[ks,:] = np.fromstring(inline,sep=' ')
    fout.write(inline)
        
    fout.write(fin.readline())
    for km in range(nm):
        inline = fin.readline()
        dat = np.fromstring(inline,sep=' ')
        dat[-1] = dat[-1]*relpert(fsdata[ks,0])
        fout.write('{:6d}  {:3d}  {:.8e}  {:.8e}  {:.8e}  {:.8e}\n'.format(
                int(dat[0]), int(dat[1]), dat[2], dat[3], dat[4],dat[5]))

plt.plot(stest,relpert(stest))
plt.plot(fsdata[:,0], 1.0/fsdata[:,1])
#plt.plot([0,1],[5.0/4.0,5.0/4.0])
#plt.plot([0,1],[4.0/3.0,4.0/3.0])
#plt.plot([0,1],[3.0/2.0,3.0/2.0])
#plt.plot([0,1],[2,2])


#CC Boozer-coordinate data file Version 0.1 Author J.Geiger Created: 26.07.2010
#CC based on calculation: tok_examples/boozer.tok02
#CC <beta> = 0.006432
#CC Configurationname = tok02, aspect ratio ca 3.8, q_0 > 1, q_a < 3
# m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]
# 36    0    63   1   -1.330000e+00   0.46000   1.64377
#       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)
#                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper
#  2.3438E-02  8.9997E-01 -1.7754E+07 -2.3937E+04 -1.7231E+04 -7.7240E+00
#    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     bmn/[T]
#   -18    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.81554063e-18
#   -17    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -7.88256975e-19
#   -16    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.51987273e-18
#   -15    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -6.72614485e-19
#   -14    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.07698780e-19
#   -13    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.07644868e-17
#   -12    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.44805638e-16
#   -11    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.38982449e-15
#   -10    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.11444928e-14
#    -9    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.11606197e-13
#    -8    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.22676736e-12
#    -7    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  7.75435510e-12
#    -6    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -6.67967800e-12
#    -5    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.86756216e-10
#    -4    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.86813817e-09
#    -3    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.07059940e-09
#    -2    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.62253553e-07
#    -1    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -3.77039425e-05
#     0    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.96335464e-03
#     1    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -3.77039425e-05
#     2    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.62253553e-07
#     3    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.07059940e-09
#     4    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.86813817e-09
#     5    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.86756216e-10
#     6    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -6.67967800e-12
#     7    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  7.75435510e-12
#     8    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.22676736e-12
#     9    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.11606197e-13
#    10    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.11444928e-14
#    11    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.38982449e-15
#    12    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.44805638e-16
#    13    3  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.07644868e-17
#    14    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -4.07698780e-19
#    15    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -6.72614485e-19
#    16    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -1.51987273e-18
#    17    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -7.88256975e-19
#    18    3  0.00000000e+00  0.00000000e+00  0.00000000e+00 -2.81554063e-18
