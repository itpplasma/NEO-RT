#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
import re
#from exportfig import exportfig
#from noexportfig import exportfig

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
smin = 0.0
smax = 1.0

#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_CIRC/170608_s_ql_quadpack'
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_CIRC/170608_s_ql'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_CIRC/170608_s_ql_finite'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_CIRC/170608_s_ql_full'

datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170118_n1_ql'
#datadirs[0] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170608_n1_ql_finite'
#datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170608_n1_ql_full' 
datadirs[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170608_n1_ql_transp' 

for k in datadirs.keys():
    datadir = datadirs[k]
    pattern = re.compile(r'driftorbit([0-9]+)\.out')
    files = os.listdir(datadir)
    files = [f for f in files if pattern.match(f)]

    infiles = [f.replace('out','in') for f in files]
    s[k] = []
    sbdata[k] = []
    data[k] = []
    
    kf = -1
    for f in files:
        kf = kf+1
        dat = np.loadtxt(os.path.join(datadir,f))
        if dat.size == 0:
            continue
        dat2 = np.loadtxt(os.path.join(datadir,f.replace('.out','_integral.out')))
        if dat2.size == 0:
            continue
            
        data[k].append(dat)
        sbdata[k].append(dat2)
                
        fp = open(os.path.join(datadirs[k], infiles[kf]))
        lines = fp.readlines()
        s[k].append(float(lines[3].split()[0]))
        fp.close()
    
    s[k] = np.array(s[k])
    data[k] = np.array(data[k])
    #sbdata[k] = np.concatenate(sbdata[k])
    #sbdata[k] = sbdata[k][np.abs(sbdata[k][:,1])<0.5]

    order = np.argsort(s[k])
    s[k] = s[k][order]
    data[k] = data[k][order,:]
    #sbdata[k] = sbdata[k][order,:]
    condi = (s[k]<smax)*(s[k]>smin)

    data[k] = data[k][condi,:]
    #sbdata[k] = sbdata[k][condi,:]
    s[k] = s[k][condi]

    D11do[k] = data[k][:,4]
    D12do[k] = data[k][:,8]
    D12D11do[k] = data[k][:,8]/data[k][:,4]
    #D11sb[k] = sbdata[k][:,4]
    #D12sb[k] = sbdata[k][:,8]


plt.figure(1)
plt.semilogy(s[1], D11do[1], 'rx-')
plt.semilogy(s[0], D11do[0], '.--')
plt.ylim([2e-5,2e-2])



