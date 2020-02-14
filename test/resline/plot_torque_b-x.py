#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import scipy.interpolate as spi
import scipy.integrate as spint
#from exportfig import exportfig
#from noexportfig import exportfig

plt.close('all')

data = {}
datab = {}
s = {}
q = {}
rhopol = {}
dVds = {}
sbox = {}
ds = {}
Mt = {}
dTphi_int_ds = {}
dTphi_int_box = {}
Tphi = {}
Tphi_box = {}

datadirs = {}
datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170613_n1_ql_torque_full_190'
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_CIRC/170613_s_ql_torque_full'
pattern = re.compile(r'driftorbit([0-9]+)_torque\.out')

smin = 0.0
smax = 0.77

profdata = np.loadtxt(os.path.join(datadirs[0],'profile.in'))
s2 = profdata[:,0]
sbox2 = np.append(0.0,np.append(s2[:-1]+np.diff(s2)/2.0,1.0))
ds2 = np.diff(sbox2)

profdata3 = np.loadtxt('/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170613_n1_ql_torque_finite_190/driftorbit1_profile.out')
s3 = profdata3[:,0]
dVds3 = profdata3[:,2]
q3 = profdata3[:,3]

iotaspl = spi.splrep(s3,1.0/q3)
iotaint = spint.quad(lambda x: spi.splev(x, iotaspl), 0.0, 1.0)[0]
dVdsspl = spi.splrep(s3,dVds3*1e-6)

for k in range(len(datadirs)):
    datadir = datadirs[k]
    files = os.listdir(datadir)
    files = [f for f in files if pattern.match(f)]
    files2 = [f.replace('torque','magfie_param') for f in files]
    filesco = [f.replace('torque','torque_box_co') for f in files]
    filesctr = [f.replace('torque','torque_box_ctr') for f in files]
    filest = [f.replace('torque','torque_box_t') for f in files]
    data[k] = []
    datab[k] = []
    
    q[k] = []
    kf = -1
    for f in files:
        kf = kf+1
        dat = np.loadtxt(os.path.join(datadir,f))
        if dat.size == 0:
            continue
        
        with open(os.path.join(datadir,files2[kf])) as f2:
            for kl in range(12):
                f2.readline()
            dat2 = f2.readline()
            q[k].append(float(dat2.split()[-1]))
            
        data[k].append(dat) 
        
        dat = np.loadtxt(os.path.join(datadir,filesco[kf]))
        if dat.size > 0:
            datab[k].append(np.sum(dat,0))
        dat = np.loadtxt(os.path.join(datadir,filesctr[kf]))
        if dat.size > 0:
            datab[k][-1] = datab[k][-1] + np.sum(dat,0)
        dat = np.loadtxt(os.path.join(datadir,filest[kf]))
        if dat.size > 0:
            datab[k][-1] = datab[k][-1] + np.sum(dat,0)
    
    data[k] = np.array(data[k])
    datab[k] = np.array(datab[k])
    #datab[k] = datab[k]*du # BUG    
    
    q[k] = np.array(q[k])
    order = np.argsort(data[k][:,0])
    data[k] = data[k][order,:]
    q[k] = q[k][order]     
    
    condi = (data[k][:,0]<smax)*(data[k][:,0]>smin)   
    data[k] = data[k][condi]
    datab[k] = datab[k][order][condi]
    q[k] = q[k][condi]
    
    
    s[k] = data[k][:,0]
    rhopol[k] = np.sqrt([spi.splint(0.0, sk, iotaspl)/iotaint for sk in s[k]])   
    
    dVds[k] = data[k][:,1]*1e-6
    sbox[k] = np.append(smin,np.append(s[k][:-1]+np.diff(s[k])/2.0,smax))
    ds[k] = np.diff(sbox[k])
    Mt[k] = data[k][:,2]
    
    dTphi_int_ds[k] = np.sum(data[k][:,3:6],1)*1e-7
    dTphi_int_box[k] = np.dot(ds[k],datab[k][:,1:])*1e-7
    Tphi[k] = dTphi_int_ds[k]/dVds[k]
    Tphi_box[k] = dTphi_int_box[k][1:-1]/(ds2*spi.splev(s2,dVdsspl))
    
rhopol2 = np.sqrt([spi.splint(0.0, sk, iotaspl)/iotaint for sk in s2])  

plt.figure(1)
plt.semilogy(rhopol[0], Tphi[0])
plt.semilogy(rhopol2, Tphi_box[0])

plt.figure(2)
plt.plot(rhopol[0], np.cumsum(dTphi_int_ds[0]*ds[0]))

plt.plot(rhopol2, np.cumsum(dTphi_int_box[0][1:-1]))

plt.figure(3)
plt.semilogy(rhopol[0], dTphi_int_ds[0]*ds[0])
plt.semilogy(rhopol2, dTphi_int_box[0][1:-1])
plt.show()

print('Integral torque        : {} Nm'.format(np.sum(dTphi_int_ds[0]*ds[0])))
print('Integral torque box    : {} Nm'.format(np.sum(dTphi_int_box[0][1:-1])))
print('Integral torque box inf: {} Nm'.format(np.sum(dTphi_int_box[0])))

