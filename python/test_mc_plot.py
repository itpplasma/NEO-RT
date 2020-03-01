# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fort.3000')
tau = data[:, 0]
zs = data[:, 1:]

plt.figure()
plt.plot(zs[:, 0]*np.cos(zs[:, 2]), zs[:, 0]*np.sin(zs[:, 2]))
plt.show()
