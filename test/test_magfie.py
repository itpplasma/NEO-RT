import numpy as np
import matplotlib.pyplot as plt

files = ['../driftorbit_test_magfie.out', '../driftorbit_test_magfie_neo.out']
#files[1] = '/temp/ert/CONDOR/driftorbit/RUNS_ASDEX/2016-08-24_n1_shear_newmagfie/'\
#    + 'driftorbit0_magfie.out'

plt.figure()
f = files[0]
data = np.loadtxt(f)
for k in range(data.shape[1]):
    plt.subplot(4,4,k+1)
    plt.plot(data[:,0],data[:,k])
    plt.title(k)
        
f = files[1]
data = np.loadtxt(f)
for k in range(data.shape[1]):
    plt.subplot(4,4,k+1)
    plt.plot(data[:,0],data[:,k],'--')
    plt.title(k)