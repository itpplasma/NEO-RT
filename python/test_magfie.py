# %%
import numpy as np
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt, magfie  # , magfie_pert

# %%
magfie.inp_swi = 8
magfie.s = 0.15
magfie.do_magfie_init()
psipr = magfie.psi_pr

# magfie_pert.do_magfie_pert_init()

# %%
x = np.zeros(3)
bmod = libneo_rt.new('double', 0.0)
sqrtg = libneo_rt.new('double', 0.0)
bder = np.zeros(3)
hcovar = np.zeros(3)
hctrvr = np.zeros(3)
hcurl = np.zeros(3)

x[0] = magfie.s  # s
x[1] = 0.0       # theta
x[2] = 0.0       # phi

magfie.do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

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

x[0] = magfie.s  # s
x[1] = 0.0       # varphi
for kth in np.arange(nth):
    x[2] = th[kth]
    magfie.do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    B[kth] = bmod[0]

plt.figure()
plt.plot(th/np.pi, B)
plt.grid()
plt.xlabel(r'$\vartheta/\pi$')
plt.ylabel(r'$B / \mathrm{G}$')
plt.title(r'$s = {}$'.format(magfie.s))

# %% For 2D plots
#
# nth = 30
# nph = 32

# th = np.linspace(0, 2*np.pi, nth)
# ph = np.linspace(0, 2*np.pi, nph)
# TH, PH = np.meshgrid(th, ph)
# B = np.empty_like(TH)
# Bpert = np.empty_like(TH)

# x[0] = magfie.s
# for kth in np.arange(nth):
#     for kph in np.arange(nph):
#         x[1] = th[kth]
#         x[2] = ph[kph]
#         magfie.do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl, psipr)
#         B[kph, kth] = bmod[0]
#         Bpert[kph, kth] = bmod[0]
# # %%

# plt.contour(TH, PH, B)
