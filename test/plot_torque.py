#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
import re

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
datadirs[0] = '~/mnt/marconi_scratch/NEO-RT/AUG/181113_test_igochine_kink/1e5'
datadirs[0] = os.path.expanduser(datadirs[0])
#datadirs[0] = r'M:\run\NEO-RT\AUG\181113_test_igochine_kink'
#datadirs[0] = r'M:\run\NEO-RT\AUG\181113_n1_ql_torque'
#datadirs[0] = r'M:\run\NEO-RT\CIRC\181109_s_ql_torque'

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
    
    dVds[k] = data[k][:,1]*1e-6
    sbox[k] = np.append(smin,np.append(s[k][:-1]+np.diff(s[k])/2.0,smax))
    ds[k] = np.diff(sbox[k])
    Mt[k] = data[k][:,2]
    dTphi_int_co_ds[k] = data[k][:,3]*1e-7
    dTphi_int_ctr_ds[k] = data[k][:,4]*1e-7
    dTphi_int_t_ds[k] = data[k][:,5]*1e-7
    dTphi_int_ds[k] = np.sum(data[k][:,3:6],1)*1e-7
    Tphi[k] = dTphi_int_ds[k]/dVds[k]
#%%
plt.figure(1)
plt.title('Shift: {}, torque: {:1.2f} Nm'.format(
            datadirs[0].split('/')[-1],
            np.sum(dTphi_int_ds[0]*ds[0])))
plt.plot(var['boozer_s'], var['TphiNA_spec'][:,1]*0.1)
plt.plot(s[0], -Tphi[0])
plt.xlabel('s')
plt.ylabel('Tphi [Nm]')
plt.ylim([-5,5])
plt.grid()
plt.legend(['NEO-2', 'NEO-RT (inverted)'])

#plt.figure(2)
#plt.plot(s[0], np.cumsum(dTphi_int_ds[0]*ds[0]))
#
#plt.figure(3)
#plt.plot(s[0], np.cumsum(dTphi_int_ds[0]*ds[0]),'k--')

plt.show()

print('Integral torque   : {} Nm'.format(np.sum(dTphi_int_ds[0]*ds[0])))


