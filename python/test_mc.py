# %%
import numpy as np
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt_mc, parmot_mod
from scipy.integrate import solve_ivp

s0 = 0.4  # Starting radial position

# %% Variables for magnetic field
x = libneo_rt_mc._ffi.new('double[3]')
bmod = libneo_rt_mc.new('double', 0.0)
sqrtg = libneo_rt_mc.new('double', 0.0)
bder = libneo_rt_mc._ffi.new('double[3]')
hcovar = libneo_rt_mc._ffi.new('double[3]')
hctrvr = libneo_rt_mc._ffi.new('double[3]')
hcurl = libneo_rt_mc._ffi.new('double[3]')


# %% Initialize magnetic field
x[0] = 1e-8      # s
x[1] = 0.0       # theta
x[2] = 0.0       # phi

libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

print(bmod[0])
print(sqrtg[0])
print(bder[0:3])
print(hcovar[0:3])
print(hctrvr[0:3])
print(hcurl[0:3])

# %% Plot magnetic field over flux surface poloidal cut
nth = 100
th = np.linspace(-np.pi, np.pi, nth)
B = np.empty_like(th)

x[0] = s0        # s
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

# %% Test orbit integration
c = 2.9979e10
e_charge = 4.8032e-10
e_mass = 9.1094e-28
p_mass = 1.6726e-24
ev = 1.6022e-12

parmot_mod.rmu = 1e5  # Inverse relativistic temperature

# This will translate to dimensionless units via parmot_mod.ro0
bmod_ref = 1e4            # 1 Tesla in Gauss
bmod00 = 1.0              # 1 Tesla in Tesla
tempi1 = 0.17549561306e4  # ion temperature
am1 = 2                   # Atomic mass 2 of deuterium ions
Zb = 1                    # Charge 1 of deuterium ions
amb = am1                 # Same ions test (beam) particles and main plasma

v0 = np.sqrt(2.0*tempi1*ev/(am1*p_mass))       # Reference (thermal) velocity
# Reference Larmor radius of thermal particles
rlarm = v0*amb*p_mass*c/(Zb*e_charge*bmod_ref)
parmot_mod.ro0 = rlarm*bmod00                  # Rescaled Larmor radius

print('ro0: {}'.format(parmot_mod.ro0))

# Initial conditions
z = libneo_rt_mc._ffi.new('double[5]')
z[0] = s0
z[1] = 0.0
z[2] = 0.7*np.pi
z[3] = 1.38  # normalized velocity module  v / v_0
z[4] = 0.00  # pitch v_\parallel / v:

x[0:3] = z[0:3]
libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

nt = 10000
dtau = 1.0

z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)

tau = libneo_rt_mc.new('double', 0.0)

# Horrible global variable, just for testing fast.
vz = libneo_rt_mc._ffi.new('double[5]')
def velo(t, y):
    z[0:5] = y[:]
    libneo_rt_mc.velo(tau, z, vz)
    return np.frombuffer(libneo_rt_mc._ffi.buffer(vz), dtype=np.float64)

times = np.linspace(0, nt*dtau, nt)

sol = solve_ivp(velo, [times[0], times[-1]], z0, t_eval=times, method='LSODA')
zs = sol.y
#%%
plt.figure()
plt.plot(zs[0,:]*np.cos(zs[2,:]), zs[0,:]*np.sin(zs[2,:]))
plt.show()
