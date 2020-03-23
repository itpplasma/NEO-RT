# %%
import numpy as np
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt_mc, parmot_mod
from scipy.integrate import solve_ivp

#
# Normalized time variable tau of integrator is given in units of length !
# 
# Physical time t = tau/v0
# where v0 is the thermal velocity (or some reference velocity)
#
# This means that after tau=1 a particle with velocity v0 
# has passed one unit length (1 cm if v0 is in CGS units)
#
# Reason: Can easily estimate integration steps for "worst case" strongly 
# passing particles that have the fastest timescale for doing one turn.
# A strongly passing particle with v0 requires time t = 2*pi*R/v0 to run in 
# a circle of cylinder radius R (e.g. the major radius of the device).
# In terms of tau it takes tau = 2*pi*R, independently from thermal velocity v0.
# Then we can set e.g. dtau = tau/N to be sure we have at least N
# integation points per turn (in-between we can also integrate adaptively) 
#

s0 = 0.4  # Starting radial position

# %% Variables for magnetic field
x = libneo_rt_mc._ffi.new('double[3]')
bmod = libneo_rt_mc.new('double', 0.0)
sqrtg = libneo_rt_mc.new('double', 0.0)
bder = libneo_rt_mc._ffi.new('double[3]')
hcovar = libneo_rt_mc._ffi.new('double[3]')
hctrvr = libneo_rt_mc._ffi.new('double[3]')
hcurl = libneo_rt_mc._ffi.new('double[3]')


# %% Initialize axisymmetric magnetic field from in_file
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

# Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
parmot_mod.rmu = 1e5  

# This will translate to dimensionless units via parmot_mod.ro0
# Here we translate an input file in Tesla to computation in Gauss/CGS
bmod_ref = 1e4            # 1 Tesla in Gauss
bmod00 = 1.0              # 1 Tesla in Tesla
tempi1 = 1.7549561306e3   # ion temperature
am = 2                    # Atomic mass 2 of deuterium ions
Zb = 1                    # Atomic charge 1 of deuterium ions

v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))        # Reference (thermal) velocity
# Reference Larmor radius of thermal particles
rlarm = v0*am*p_mass*c/(Zb*e_charge*bmod_ref)  # Larmor radius in bmod_ref
parmot_mod.ro0 = rlarm*bmod00                  # Rescaled Larmor radius

print('ro0: {}'.format(parmot_mod.ro0))

# %% Integrate with a given number of timesteps
# Initial conditions
z = libneo_rt_mc._ffi.new('double[5]')
z[0] = s0         # s = psi_tor/psi_tor_a
z[1] = 0.0        # varphi
z[2] = 0.9*np.pi  # vartheta
z[3] = 1.0        # normalized velocity module  v / v_0
z[4] = 0.00       # pitch v_\parallel / v:

x[0:3] = z[0:3]   # Initial conditions in space
# Testing if magnetic field works for initial conditions
libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

nt = 10000  # Number of timesteps
dtau = 1.0  # Maximum timestep, just setting empirically here
# usually, use estimate with major radius, as described above

# create numpy array with data at memory location of vz
z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)

# Horrible global variables, just for testing fast.
tau = libneo_rt_mc.new('double', 0.0)
vz = libneo_rt_mc._ffi.new('double[5]')  # future: vz = np.zeros(5)

def velo(t, y):
    """ To translate Fortran velo for use with numpy arrays"""
    z[0:5] = y[:]
    libneo_rt_mc.velo(tau, z, vz)
    # create numpy array with data at memory location of vz
    return np.frombuffer(libneo_rt_mc._ffi.buffer(vz), dtype=np.float64)
    # return vz  ... future usage when fffi will be better

times = np.linspace(0, nt*dtau, nt)

sol = solve_ivp(velo, [times[0], times[-1]], z0, t_eval=times, method='LSODA')

# %%  Plot orbit in the poloidal and toroidal plane
zs = sol.y
# Plot with "pseudo-toroidal" radius r_tor = sqrt(s), as the toroidal
# flux roughly scales with r^2. # Pretend major R0=2 for visualization
R0 = 2.0
rho_tor = np.sqrt(zs[0,:])
varphi = zs[1,:]
vartheta = zs[2,:]
R = R0 + rho_tor*np.cos(vartheta)
Z = rho_tor*np.sin(vartheta)
X = R*np.cos(varphi)
Y = R*np.sin(varphi)

plt.figure()  # poloidal projection
plt.plot(R, Z)

plt.figure()  # toroidal projection
plt.plot(X, Y)

# %% Integrate with events: banana tips and circular crossings of varphi=2*pi

def event_banantip(t, z):
    """ Trapped orbit - banana tip at vpar=0 """
    return z[4]

event_banantip.direction = 1.0
event_banantip.terminal = True

def event_circ(t, z):
    """ Passing orbit - full turn at varphi=2pi """
    return z[2] - 2*np.pi

event_circ.terminal = True

def run_bananas(s):
    # Initial conditions
    z[0] = s
    z[1] = 0.0
    z[2] = -0.9*np.pi
    z[3] = 1.38   # normalized velocity module  v / v_0
    z[4] = 1e-13  # pitch v_\parallel / v, be careful with events !!
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)

    integ = solve_ivp(
        velo, (0, 5e5*dtau), z0, events=(event_banantip,event_circ), 
        max_step=5*dtau, method='LSODA', rtol=1e-6, atol=1e-8)

    print(integ.t_events)

    zs = integ.y
    plt.plot(np.sqrt(zs[0,:])*np.cos(zs[2,:]), np.sqrt(zs[0,:])*np.sin(zs[2,:]))
 
    # Om^th = om_b = 2*pi/t_bounce = 2*pi*v0/tau_bounce
    omb = 2*np.pi*v0/integ.t_events[0][0]
    
    # Toroidal precession frequency
    Omphi = omb*integ.y_events[0][0][1]/(2*np.pi)
    # Om^ph = om_b*varphi(t_bounce)/theta_canonical(t_bounce)
    # Reason: During bounce, theta_canonical changes by 2*pi.
    # Ratio of precession and bounce frequency is ratio 
    # between change in varphi to change in theta_canonical.
    # see Diss_Albert for detailed explanation
    
    return [omb, Omphi] 

# Factor of larmor radius to check thin bananas
# TODO: compare to other thin orbit formula from NEO-RT driftorbit.f90
# TODO: passing orbits

fac1 = 1.0
fac2 = 1e-2

parmot_mod.ro0 = rlarm*bmod00*fac1 # Bigger Larmor radius
plt.figure(figsize=(4.0, 4.0))
freqs = []
srange = np.linspace(0.001, 0.7, 200)**2
#srange = np.array([0.002, 0.005, 0.01, 0.02, 0.03, 0.12, 0.22, 0.32, 0.52, 0.72])*0.7
for sk in srange:
    freqs.append(run_bananas(sk))
plt.xlabel(r'$\rho_\mathrm{tor}\,\cos\vartheta_\mathrm{B}$')
plt.ylabel(r'$\rho_\mathrm{tor}\,\sin\vartheta_\mathrm{B}$')
plt.title("Drift orbits")
plt.xlim([-1.0, 1.0])
plt.ylim([-1.0, 1.0])

parmot_mod.ro0 = fac2*rlarm*bmod00  # Smaller Larmor radius
plt.figure(figsize=(4.0, 4.0))
freqs_thin = []
for sk in srange:
    freqs_thin.append(run_bananas(sk))
plt.xlabel(r'$\rho_\mathrm{tor}\,\cos\vartheta_\mathrm{B}$')
plt.ylabel(r'$\rho_\mathrm{tor}\,\sin\vartheta_\mathrm{B}$')
plt.title("Thin orbit approximation")
plt.xlim([-1.0, 1.0])
plt.ylim([-1.0, 1.0])
# %%
freqs = np.array(freqs)
freqs_thin = np.array(freqs_thin)


plt.figure()
plt.plot(np.sqrt(srange), freqs_thin[:,0] , 'k--', zorder=10)
plt.plot(np.sqrt(srange), freqs[:,0], color='tab:red')
plt.xlabel(r'$\rho_\mathrm{tor}$')
plt.ylabel(r'$\omega_b / \mathrm{s}^{-1}$')
plt.title("Bounce frequency")
plt.xlim([0, 0.9])
plt.legend(['thin', 'full'])

plt.figure()
# Rescale by factor 1e-3 due to other Larmor radius
plt.plot(np.sqrt(srange), freqs_thin[:,1]*fac1/fac2, 'k--', zorder=10)
plt.plot(np.sqrt(srange), freqs[:,1], color='tab:red')
plt.xlabel(r'$\rho_\mathrm{tor}$')
plt.ylabel(r'$\Omega_{\mathrm{tB}} / \mathrm{s}^{-1}$')
plt.title("Toroidal precession frequency (magnetic)")
plt.xlim([0, 0.9])
plt.legend(['thin', 'full'])

# %%

plt.show()
