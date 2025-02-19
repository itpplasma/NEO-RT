import os
import numpy as np

os.system('rm -f *.out')
os.system('../../build/neo_rt.x driftorbit_test')

#%%
epsmn = 1e-3

data = np.loadtxt('driftorbit_test.out')
D11 = data[4]
D12 = data[8]

with open('driftorbit_test_magfie_param.out', 'r') as f:
    for line in f:
        if 'Drp' in line:
            Drp = float(line.split()[-1])
            break

print(f'NEO-RT:    D11/eps**2   = {D11/epsmn**2:.3e}')
print(f'Reference: D11rp/eps**2 = {Drp:.3e}')
print()
print(f'NEO-RT:    D12/D11     = {D12/D11:.3e}')
print(f'Reference: D12rp/D11rp = {3.0:.3e}')
print()

assert(abs(D12/D11 - 3.0)/3.0 < 0.08)

print('Ripple plateau test: OK')
