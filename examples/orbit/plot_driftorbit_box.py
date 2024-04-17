#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os
import re
#from exportfig import exportfig
#from noexportfig import exportfig

def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

plt.close('all')

datadir = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/170608_n1_ql_full' 


pattern = re.compile(r'driftorbit([0-9]+)\.out')
files = os.listdir(datadir)
files = [f for f in files if pattern.match(f)]
infiles = [f.replace('out','in') for f in files]
boxfilest = [f.replace('.out','_box_t.out') for f in files]
boxfilesco = [f.replace('.out','_box_cop.out') for f in files]
boxfilesctr = [f.replace('.out','_box_ctr.out') for f in files]

infiles = [f.replace('out','in') for f in files]
s = []
data = []

kf = -1
for f in files:
    kf = kf+1
    dat = np.loadtxt(os.path.join(datadir,f))
    if dat.size == 0:
        continue
    dat2 = np.loadtxt(os.path.join(datadir,f.replace('.out','_integral.out')))
    if dat2.size == 0:
        continue
        
    data.append(dat)
            
    fp = open(os.path.join(datadir, infiles[kf]))
    lines = fp.readlines()
    s.append(float(lines[3].split()[0]))
    fp.close()
