#!/usr/bin/python

import numpy as np
import scipy.interpolate as spi
import scipy.integrate as spint
import matplotlib.pyplot as plt
import os
import re
#from exportfig import exportfig
from noexportfig import exportfig

plt.close('all')

datadirs = {}
s = {}
data = {}
sbdata = {}
D11do = {}
D12do = {}
D12D11do = {}
D11sb = {}
D12sb = {}
D11i = {}
D12i = {}
D12D11i = {}

neodata = np.loadtxt('ntv_out_asdex30835_rmp_n2_HatFun14.dat')
datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170118_n1_noshear_ql'
datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170118_n1_noshear_nonlin'

#neodata = np.loadtxt('ntv_out_nonlocal.dat')
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/2016-08-24_n1_shear_newmagfie'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/2016-08-24_n1_shear_nonlin'
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170118_n1_ql'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170118_n1_nonlin'

#neodata = np.loadtxt('ntv_out_asdex30835_rmp_n2_HatFun14_ErOnly.dat', skiprows=1)
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/2016-08-24_n1_nomag'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/161019_n1_nomag_nonlin'

neos = neodata[:,0]
neoorder = np.argsort(neos)
avnabs = neodata[:,14]
avnabs_interp = spi.interpolate.interp1d(neos, avnabs, bounds_error=False)
neos = neos[neoorder]
neodata = neodata[neoorder,:]
D11neo = neodata[:,7]
D12neo = neodata[:,8]
D12D11neo = D12neo/D11neo


profname1 = os.path.join('/proj/plasma/RMP/DOC/Diss_Albert.orig/code/driftorbit/Mtprofile_Hamiltonian.dat')
profname2 = os.path.join('/proj/plasma/RMP/DOC/Diss_Albert.orig/code/driftorbit/ASDEX/profiles/Mtprofile190.dat')
profdata1 = np.loadtxt(profname1)
profdata2 = np.loadtxt(profname2)

for k in datadirs.keys():
    datadir = datadirs[k]
    pattern1 = re.compile(r'driftorbit([0-9]+)\.out')
    pattern2 = re.compile(r'driftorbit_*([0-9]+\.[0-9]+)\.out')
    files = os.listdir(datadir)
    files1 = [f for f in files if pattern1.match(f)] # new format (> 01/2016)
    files2 = [f for f in files if pattern2.match(f)] # old format

    if len(files1) > 0:
        files = files1
        infiles = [f.replace('out','in') for f in files]
        s[k] = []
        for f in infiles:
            fp = open(os.path.join(datadirs[k], f))
            lines = fp.readlines()
            s[k].append(float(lines[3].split()[0]))
            fp.close()
        s[k] = np.array(s[k])
    else:
        files = files2
        s[k] = np.array([float(pattern2.match(f).groups(1)[0]) for f in files])

    sbdata[k] = []
    data[k] = []
    for f in files:
        #print(f)
        data[k].append(np.loadtxt(os.path.join(datadir,f)))
        sbdata[k].append(np.loadtxt(os.path.join(datadir,f.replace('.out','_integral.out'))))

    data[k] = np.array(data[k])
    sbdata[k] = np.concatenate(sbdata[k])
    sbdata[k] = sbdata[k][np.abs(sbdata[k][:,1])<0.5]

    order = np.argsort(s[k])
    s[k] = s[k][order]
    data[k] = data[k][order,:]
    sbdata[k] = sbdata[k][order,:]

    # correction of avnabs = dsdreff = <|nabla s|>
    a = 49.41
    avnabs_old = 2.0/a*np.sqrt(s[k])
    for l in range(1, data[k][1,:].size):
        data[k][:,l] = data[k][:,l]*(avnabs_old/avnabs_interp(s[k]))**2
    for l in range(2,9):
        sbdata[k][:,l] = sbdata[k][:,l]*(avnabs_old/avnabs_interp(s[k]))**2

    if len(files1) > 0: # new file format
        D11do[k] = data[k][:,4]
        D12do[k] = data[k][:,8]
        D12D11do[k] = data[k][:,8]/data[k][:,4]
        D11sb[k] = sbdata[k][:,4]
        D12sb[k] = sbdata[k][:,8]
    else:
        D11do[k] = data[k][:,3]
        D12do[k] = data[k][:,6]
        D12D11do[k] = data[k][:,6]/data[k][:,3]

    D11i[k] = spi.interpolate.interp1d(s[k], D11do[k], bounds_error=False)
    D12i[k] = spi.interpolate.interp1d(s[k], D12do[k], bounds_error=False)
    D12D11i[k] = spi.interpolate.interp1d(s[k], D12D11do[k], bounds_error=False)
    
sMt1  = profdata1[:,0]
Mt1  = profdata1[:,1]
sMt2  = profdata2[:,0]
Mt2  = profdata2[:,1]

#facnew = 46./49.41
#facnew = 6816./6845.*46./49.41
#facnew = (46./49.41)**2
#facnew = 1.0
D11do[0] = D11do[0]
D11do[1] = D11do[1]
D12do[0] = D12do[0]
D12do[1] = D12do[1]

#%%%
qe = 4.8032e-10 # elementary charge
c=2.9979e10 # speed of light
Zi=1; # ion charge number (deuterium)
mi=3.3436e-24 # ion mass (deuterium)
me=9.1094e-28 # electron mass

iota = spi.interp1d(neos, neodata[:,10])
R0 = neodata[0,11]
def q(s) : return 1.0/iota(s)
Bref = spi.interp1d(neos, neodata[:,12])
avnabpsi = spi.interp1d(neos, neodata[:,14])
sqrtgBth = spi.interp1d(neos, neodata[:,15])
Mt = spi.interp1d(sMt2, Mt2)

profdata = np.loadtxt('kappa_profile_asdex30835_axi.dat', skiprows=1)
sprof = profdata[:,1]
spol = spi.interp1d(sprof, profdata[:,0])
ne = spi.interp1d(sprof, profdata[:,9])
Ti = spi.interp1d(sprof, profdata[:,11]*(qe))
dlogneds = spi.interp1d(sprof, profdata[:,12])
dlogTids = spi.interp1d(sprof, profdata[:,14])
def vTi(s) : return np.sqrt(2.0*Ti(s)/mi)
def rholi(s) : return (mi*c*vTi(s))/(qe*Bref(s))

def Dp(s) : return np.pi*q(s)*vTi(s)*(rholi(s)**2)/(16.0*R0)
def Er(s) : return Mt(s)*vTi(s)*sqrtgBth(s)/(R0*c)
def A1(s) : return avnabpsi(s)*dlogneds(s)-qe*Er(s)/Ti(s)-1.5*avnabpsi(s)*dlogTids(s)
def A2(s) : return avnabpsi(s)*dlogTids(s)

Tphi0=-(sqrtgBth(s[0])/c)*qe*(-ne(s[0])*(D11do[0]*Dp(s[0])*A1(s[0])+D12do[0]*Dp(s[0])*A2(s[0])))
Tphi1=-(sqrtgBth(s[1])/c)*qe*(-ne(s[1])*(D11do[1]*Dp(s[1])*A1(s[1])+D12do[1]*Dp(s[1])*A2(s[1])))
Tphineo=-(sqrtgBth(neos)/c)*qe*(-ne(neos)*(D11neo*Dp(neos)*A1(neos)+D12neo*Dp(neos)*A2(neos)))


#%%


plt.figure(2)
plt.semilogy(np.sqrt(neos), np.abs(Tphineo)/10.0, '-', linewidth=2, mew=1, ms=5)
plt.semilogy(np.sqrt(s[0]), np.abs(Tphi0)/10.0, '-.', linewidth=2, mew=1, ms=5) 
plt.semilogy(np.sqrt(s[1]), np.abs(Tphi1)/10.0, '--', linewidth=2, mew=1, ms=5)
plt.ylim([1e-3,1e0])
plt.grid(True, which="both")
plt.xlabel(r'Normalized effective radius')
plt.ylabel(r'Torque density / N/m^2')
plt.legend(['NEO-2 quasilinear','NEO-RT quasilinear','NEO-RT non-linear'],loc='upper left')
ax2 = plt.gca().twinx()
ax2.plot(np.sqrt(sMt1), q(sMt1), 'b')
#ax2.plot(np.sqrt(sMt2), Mt2, 'r--')
ax2.set_ylabel('safety factor', color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
#plt.title('Electric drift only')

s[1] = s[1][np.invert(np.isnan(Tphi1))]
Tphi1 = Tphi1[np.invert(np.isnan(Tphi1))]

Tphi0i = spi.InterpolatedUnivariateSpline(s[0], np.abs(Tphi0)*np.sqrt(s[0]))
Tphi1i = spi.InterpolatedUnivariateSpline(s[1], np.abs(Tphi1)*np.sqrt(s[1]))
Tphineoi = spi.InterpolatedUnivariateSpline(neos, np.abs(Tphineo)*np.sqrt(neos))


#smin = np.min(s[0])
#smax = np.max(s[0])
smin = 0.025
smax = 0.78
si = np.linspace(smin,smax,200)


#plt.figure(3)
#plt.semilogy(neos, np.abs(Tphineo)/10.0, '-', linewidth=2, mew=1, ms=5)
#plt.semilogy(s[0], np.abs(Tphi0)/10.0, 'o--', linewidth=2, mew=1, ms=5) 
#plt.semilogy(s[1], np.abs(Tphi1)/10.0, 'd--', linewidth=2, mew=1, ms=5)
#plt.ylim([1e-3,1e0])
#plt.grid(True, which="both")
#plt.xlabel(r's')
#plt.ylabel(r'Torque density / N/m^2')
#plt.legend(['NEO-2 quasilinear','NEO-RT quasilinear','NEO-RT non-linear'],loc='upper left')
#ax2 = plt.gca().twinx()
#ax2.plot(sMt1, q(sMt1), 'b')
#ax2.plot(np.sqrt(sMt1), Mt1, 'b')
##ax2.plot(np.sqrt(sMt2), Mt2, 'r--')
#ax2.set_ylabel('safety factor', color='b')
#for tl in ax2.get_yticklabels():
#    tl.set_color('b')
#plt.title('Full magnetic drift')
#ax2.set_ylim([-0.3, 0.3])

#plt.figure(3)
#plt.plot(si, Tphineoi(si))
#plt.plot(si, Tphi0i(si))
#plt.plot(si, Tphi1i(si))

Tneo = spint.quad(Tphineoi, smin, smax)
T0 = spint.quad(Tphi0i, smin, smax)
T1 = spint.quad(Tphi1i, smin, smax)

print(Tneo)
print(T0)
print(T1)