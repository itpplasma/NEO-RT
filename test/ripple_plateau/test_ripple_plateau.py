import os
import numpy as np

os.system('rm *.out')
os.system('../../BUILD/neo_rt.x driftorbit_test')

#%%
epsmn = 1e-3

data = np.loadtxt('driftorbit_test.out')
D11 = data[4]
D12 = data[8]

print(f'D11/epsmn**2 = {D11/epsmn**2}')
print(f'D12/epsmn**2 = {D12/epsmn**2}')
print(f'D12/D11 = {D12/D11}')

with open('driftorbit_test_magfie_param.out', 'r') as f:
    for line in f:
        if 'Drp' in line:
            Drp = float(line.split()[-1])
            break

print(f'Drp = {Drp}')

assert(abs(D12/D11 - 3.0) < 0.15)
assert(abs((D11/epsmn**2 - Drp)/Drp) < 0.05)

print('Ripple plateau test: OK')
