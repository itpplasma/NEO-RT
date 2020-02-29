# %%
import numpy as np
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt_mc

print(libneo_rt_mc.loaded)

x = np.zeros(3)
bmod = libneo_rt_mc.new('double', 0.0)
sqrtg = libneo_rt_mc.new('double', 0.0)
bder = np.zeros(3)
hcovar = np.zeros(3)
hctrvr = np.zeros(3)
hcurl = np.zeros(3)

x[0] = 1e-8      # s
x[1] = 0.0       # theta
x[2] = 0.0       # phi

libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

print(bmod[0])
print(sqrtg[0])
print(bder)
print(hcovar)
print(hctrvr)
print(hcurl)

# %%
nth = 100
th = np.linspace(-np.pi, np.pi, nth)
B = np.empty_like(th)

x[0] = 0.15      # s
x[1] = 0.0       # varphi
for kth in np.arange(nth):
    x[2] = th[kth]
    libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    B[kth] = bmod[0]

plt.figure()
plt.plot(th/np.pi, B)
plt.grid()
plt.xlabel(r'$\vartheta/\pi$')
plt.ylabel(r'$B / \mathrm{G}$')
plt.title(r'$s = {}$'.format(x[0]))
plt.show()
