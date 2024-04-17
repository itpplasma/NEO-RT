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
s = {}
q = {}
rhopol = {}
dVds = {}
sbox = {}
ds = {}
Mt = {}
dTphi_int_co_ds = {}
dTphi_int_ctr_ds = {}
dTphi_int_t_ds = {}
dTphi_int_ds = {}
Tphi = {}

datadirs = {}
datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170613_n1_ql_torque_190'
datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170613_n1_ql_torque_full_190'
pattern = re.compile(r'driftorbit([0-9]+)_torque\.out')

smin = 0.0
smax = 0.77

for k in range(len(datadirs)):
    datadir = datadirs[k]
    files = os.listdir(datadir)
    files = [f for f in files if pattern.match(f)]
    files2 = [f.replace('torque','magfie_param') for f in files]
    data[k] = []
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
    
    data[k] = np.array(data[k])
    q[k] = np.array(q[k])
    order = np.argsort(data[k][:,0])
    data[k] = data[k][order,:]
    q[k] = q[k][order]       
    
    condi = (data[k][:,0]<smax)*(data[k][:,0]>smin)   
    data[k] = data[k][condi]
    q[k] = q[k][condi]
    
    
    s[k] = data[k][:,0]
    rhopol[k] = np.sqrt([spi.splint(0.0, sk, iotaspl)/iotaint for sk in s[k]])   
    
    dVds[k] = data[k][:,1]*1e-6
    sbox[k] = np.append(smin,np.append(s[k][:-1]+np.diff(s[k])/2.0,smax))
    ds[k] = np.diff(sbox[k])
    Mt[k] = data[k][:,2]
    dTphi_int_co_ds[k] = data[k][:,3]*1e-7
    dTphi_int_ctr_ds[k] = data[k][:,4]*1e-7
    dTphi_int_t_ds[k] = data[k][:,5]*1e-7
    dTphi_int_ds[k] = np.sum(data[k][:,3:6],1)*1e-7
    Tphi[k] = dTphi_int_ds[k]/dVds[k]
    
plt.figure(1)
plt.semilogy(rhopol[0], Tphi[0])
plt.semilogy(rhopol[1], Tphi[1])
plt.semilogy(rhopol2, Tphi_box[0])

plt.figure(2)
plt.plot(rhopol[0], np.cumsum(dTphi_int_ds[0]*ds[0]))
plt.plot(rhopol[1], np.cumsum(dTphi_int_ds[1]*ds[1]))
rhopol2 = np.sqrt([spi.splint(0.0, sk, iotaspl)/iotaint for sk in s2])  
plt.plot(rhopol2, np.cumsum(dTphi_int_box[0][1:-1]))

plt.figure(3)
plt.plot(rhopol[0], np.cumsum(dTphi_int_ds[0]*ds[0]),'k--')
plt.plot(rhopol2, np.cumsum(dTphi_int_box[0][1:-1]),'k-')

plt.show()

print('Integral torque 0   : {} Nm'.format(np.sum(dTphi_int_ds[0]*ds[0])))
print('Integral torque 1   : {} Nm'.format(np.sum(dTphi_int_ds[1]*ds[1])))
print('Integral torque box : {} Nm'.format(np.sum(dTphi_int_box[0][1:-1])))

