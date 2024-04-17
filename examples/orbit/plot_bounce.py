#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 17:31:12 2017

@author: ert
"""

import numpy as np
import matplotlib.pyplot as plt
import os

folder = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170613_n1_ql_torque_full_190'
prefix = os.path.join(folder, 'driftorbit28_')

#prefix = 'orbit_test_'
plt.close()


plt.figure()
data = np.loadtxt(prefix+'bounce_zeroorder.out')
plt.plot(np.sqrt(data[:,1])*np.cos(data[:,3]), np.sqrt(data[:,1])*np.sin(data[:,3]),'k:')
plt.plot(np.sqrt(data[:,8])*np.cos(data[:,3]), np.sqrt(data[:,8])*np.sin(data[:,3]),'r-')
data = np.loadtxt(prefix+'bounce_full.out')
plt.plot(np.sqrt(data[:,1])*np.cos(data[:,3]), np.sqrt(data[:,1])*np.sin(data[:,3]),'k--')
data = np.loadtxt(prefix+'bounce_rela.out')
plt.plot(np.sqrt(data[:,1])*np.cos(data[:,3]), np.sqrt(data[:,1])*np.sin(data[:,3]),'k-.')

plt.xlim(-1,1)
plt.ylim(-1,1)
plt.axis('equal')
#%%

plt.figure()
data = np.loadtxt(prefix+'bounce_full.out')
q = data[:,5]
plt.plot(data[:,3]/(2*np.pi), (data[:,2]-q*data[:,3])/(2*np.pi),'k')
data = np.loadtxt(prefix+'bounce_zeroorder.out')
plt.plot(data[:,3]/(2*np.pi), (data[:,2]-q*data[:,3]+data[:,7]*data[:,0]*data[:,6])/(2*np.pi),'k--')
plt.plot(data[:,3]/(2*np.pi), (data[:,2]-data[:,5]*data[:,3]+data[:,7]*data[:,0]*data[:,6])/(2*np.pi),'k:')

#plt.figure()
#plt.plot(data[:,0], data[:,4])
#plt.plot(data[:,0], data[:,8])
#plt.plot(data[:,0], data[:,12])
#%%

def initvars(datafile):
    global data,sz,taub,tau,thorb,ph,q
    data = np.loadtxt(datafile)
    sz = np.shape(data[:,0])
    taub = data[0,6]
    tau = data[:,0]*taub
    thorb = data[:,3]-data[0,3]
    ph = data[:,2]
    q = data[:,5]

plt.figure()
initvars(prefix+'bounce_full.out')
Omcan = np.mean(np.diff(ph)/np.diff(tau))
deltaphi = (ph-Omcan*tau)
phfull = ph
sfullmean = np.mean(data[:,1])
initvars(prefix+'bounce_zeroorder.out')
Omcan2 = data[0,7] + q*2*np.pi/taub
plt.plot(thorb, Omcan2*tau - Omcan*tau,'k-')
szeromean = np.mean(data[:,1])

plt.figure()

plt.plot(np.sqrt(data[:,8]),'r-')