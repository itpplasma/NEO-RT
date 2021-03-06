# %%
import numpy as np
from random import random
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt_mc, parmot_mod
from scipy.interpolate import UnivariateSpline
#from interpolation import interpol, suppvec2
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import Rbf, splrep, splev, sproot, InterpolatedUnivariateSpline
from scipy.integrate import solve_ivp
from scipy.interpolate import UnivariateSpline
#from interpolation import interpol, suppvec2
from IPython import get_ipython
import os

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
x_cyl = libneo_rt_mc._ffi.new('double[3]')
bmod = libneo_rt_mc.new('double', 0.0)
sqrtg = libneo_rt_mc.new('double', 0.0)
bder = libneo_rt_mc._ffi.new('double[3]')
hcovar = libneo_rt_mc._ffi.new('double[3]')
hctrvr = libneo_rt_mc._ffi.new('double[3]')
hcurl = libneo_rt_mc._ffi.new('double[3]')


# Store s, iota from in_file
data_s_iota = np.empty((0,2),float)     # [s, iota]
line_nr = -1
line_nr_flux = -1
idx = 0
with open(os.path.expanduser('~/src/NEO-RT/python/in_file'), 'r') as file:
    for k in range(4): next(file)
    for line, data in enumerate(file):
        data = data.split()
        if line_nr == line:
            data_s_iota = np.append(data_s_iota, np.array([[float(data[0]), float(data[1])]]), axis=0)
        if line_nr_flux == line:
            flux = float(data[idx])
        if data[0] == 's':
            line_nr = line + 2
        if 'flux/[Tm^2]' in data:
            idx = data.index('flux/[Tm^2]')
            line_nr_flux = line + 1

psi_pr = abs(flux)/(2*np.pi)*1e8        # *1e8 bec. conversion from T*m² to Gauss*cm²
spl = UnivariateSpline(data_s_iota[:,0], psi_pr*data_s_iota[:,1],s=0)
psi_pol = spl.antiderivative()

# %% Initialize magnetic field
x[0] = 1e-8      # s
x[1] = 0.0       # theta
x[2] = 0.0       # phi

libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

#print(bmod[0])
#print(sqrtg[0])
#print(bder[0:3])
#print(hcovar[0:3])
#print(hctrvr[0:3])
#print(hcurl[0:3])

# %% Plot magnetic field over flux surface poloidal cut
nth = 100
th = np.linspace(-np.pi, np.pi, nth)
B = np.empty_like(th)
R = np.empty_like(th)
P = np.empty_like(th)
Z = np.empty_like(th)

x[0] = s0        # s
x[1] = 0.0       # varphi
for kth in np.arange(nth):
    x[2] = th[kth]
    libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    B[kth] = bmod[0]
    libneo_rt_mc.cyl_coord(x, x_cyl)
    R[kth] = x_cyl[0]
    P[kth] = x_cyl[1]
    Z[kth] = x_cyl[2]

plt.figure()
plt.plot(th/np.pi, B)
plt.grid()
plt.xlabel(r'$\vartheta/\pi$')
plt.ylabel(r'$B / \mathrm{T}$')
plt.title(r'$s = {}$'.format(x[0]))
plt.show()

#%%
plt.figure()
plt.plot(th/np.pi, R-150)
plt.plot(th/np.pi, 1e3*P)
plt.plot(th/np.pi, Z)
plt.grid()
plt.xlabel(r'$\vartheta/\pi$')
plt.title(r'$s = {}$'.format(x[0]))
plt.legend(['R - 150cm', '1e3*PHI', 'Z'])
plt.show()

# %% Test orbit integration
c = 2.9979e10           # cm/s
e_charge = 4.8032e-10   # franklin ( = 3.336e-10C)
e_mass = 9.1094e-28     # g
p_mass = 1.6726e-24     # g
ev = 1.6022e-12         # erg ( = 1e-7J)

# Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
parmot_mod.rmu = 1e5

# This will translate to dimensionless units via parmot_mod.ro0
# Here we translate an input file in Tesla to computation in Gauss/CGS
bmod_ref = 1e4            # 1 Tesla in Gauss
bmod00 = 1.0              # 1 Tesla in Tesla
tempi1 = 0.17549561306e4  # ion temperature
am = 2                    # Atomic mass 2 of deuterium ions
Zb = 1                    # Atomic charge 1 of deuterium ions

v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))        # Reference (thermal) velocity
# Reference Larmor radius of thermal particles
rlarm = v0*am*p_mass*c/(Zb*e_charge*bmod_ref)  # Larmor radius in bmod_ref
parmot_mod.ro0 = rlarm*bmod00                  # Rescaled Larmor radius

#print('ro0: {}'.format(parmot_mod.ro0))

# %% Integrate with a given number of timesteps
# Initial conditions
z = libneo_rt_mc._ffi.new('double[5]')
z[0] = s0         # s = psi_tor/psi_tor_a
z[1] = 0.0        # varphi
z[2] = 0.7*np.pi  # vartheta
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
#    print(np.frombuffer(libneo_rt_mc._ffi.buffer(vz), dtype=np.float64))
    return np.frombuffer(libneo_rt_mc._ffi.buffer(vz), dtype=np.float64)
    # return vz  ... future usage when fffi will be better

def velo_thin(t,y,v,eta):
    ret = np.zeros(5)
    libneo_rt_mc.magfie([y[0],y[1],y[2]], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    ret[0] = 0
    ret[1] = y[4]*v*hctrvr[1]               # varphi
    ret[2] = y[4]*v*hctrvr[2]               # vartheta
    ret[3] = 0
    ret[4] = -eta*v*hctrvr[2]*bder[2]       # vpar/v
    return ret

times = np.linspace(0, nt*dtau, nt)
sol = solve_ivp(velo, [times[0], times[-1]], z0, t_eval=times, method='LSODA')
# %%  Plot orbit in the poloidal and toroidal plane
zs = sol.y
# Plot with "pseudo-toroidal" radius r_tor = sqrt(s), as the toroidal
# flux roughly scales with r^2. # Pretend major R0=2 for visualization
R0 = 2.0
rho_tor = np.sqrt(zs[0,:])
#varphi = zs[1,:]
#vartheta = zs[2,:]
#R = R0 + rho_tor*np.cos(vartheta)
#Z = rho_tor*np.sin(vartheta)
#X = R*np.cos(varphi)
#Y = R*np.sin(varphi)
#
#plt.figure()  # poloidal projection
#plt.plot(R, Z)
#
#plt.figure()  # toroidal projection
#plt.plot(X, Y)

# %% Integrate with events: banana tips and circular crossings of varphi=2*pi

def event_bananatip(t, z):
    """ Trapped orbit - banana tip at vpar=0 """
    return z[4]

def event_circ(t, z):
    """ Passing orbit - full turn at theta=2pi """
    return z[2] - 2*np.pi


def orbit(r0=0.3, th0=-0.5*np.pi, v=1.38, vpar=0., plotting=False, method='thin'):
    """
    Calculate [omb, Omphi, pphi, H, mu, orbit_type]
    Methods: {'full', 'thin'}
    """

    # Initial values
    z = libneo_rt_mc._ffi.new('double[5]')
    z[0] = r0           # s
    z[1] = 0.           # varphi
    z[2] = th0          # vartheta
    z[3] = v            # normalized velocity module (v/v0)
    z[4] = vpar         # normalized parallel velocity (vpar/v)
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)


    # Calculate constants of motion
    libneo_rt_mc.magfie(z[0:3], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    v = v*v0; vpar = z[4]*v; vperp = np.sqrt(v**2 - vpar**2)    # cm/s
    mu = am*p_mass*vperp**2/(2*bmod[0]*ev)                      # eV/T
    H = am*p_mass*v**2/(2*ev)                                   # eV

    # Integrate equations of motion for thin or full orbits:
    if method == 'thin':
        eta = vperp**2/(v**2*bmod[0])     # T^(-1)
        zdot = lambda t, y: velo_thin(t,y,v,eta)
        integ = solve_ivp(zdot, (0, 1e5*dtau/v0), z0, events=(event_bananatip,event_circ),
                          max_step=10*dtau/v0, method='LSODA')
    else:
        integ = solve_ivp(velo, (0, 5e5*dtau), z0, events=(event_bananatip,event_circ),
                          max_step=10*dtau, method='LSODA')


    zs = integ.y
    time = integ.t
    y_events = integ.y_events
    time_events = integ.t_events
    rphi = zs[0,-1]                         # tip radius
    pphi = - Zb*e_charge*psi_pol(rphi)/c    # eVs/cm


    # Integrate again from tip or theta=0 for a complete orbit:
    if len(time_events[1])==0:
        zs2 = zs[:,-1]
        zs2[2] = zs2[2]%(2*np.pi)
        if zs2[2] > np.pi: zs2[2] = -(2*np.pi - zs2[2])     # project angle theta back in (-pi,pi) (as to not wrongly trigger event_circ)
        zs2[4] = 1e-15                                      # set vpar to numerical zero at tip
        if method == 'thin':
            integ2 = solve_ivp(zdot, (0, 1e5*dtau/v0), zs2, events=(event_bananatip,event_circ),
                               max_step=10*dtau/v0, method='LSODA')
        else:
            integ2 = solve_ivp(velo, (0, 5e5*dtau), zs2, events=(event_bananatip,event_circ),
                               max_step=10*dtau, method='LSODA')
        time_events = integ2.t_events
        time = integ2.t
        y_events = integ2.y_events
        zs = integ2.y

    # Plot the orbit:
    if plotting == True:
        plt.plot(100*np.sqrt(zs[0,:])*np.cos(zs[2,:]), 100*np.sqrt(zs[0,:])*np.sin(zs[2,:]))
        plt.title('Poloidal Projection (flux coord. field)')
        plt.xlabel('r / cm')
        plt.ylabel('z / cm')
        plt.axis('equal')
#        plt.show()

    # Calculate bounce/drift frequency
    if len(time_events[0])>0:
        omb = 2*np.pi/time_events[0][0]             # bounce frequency
        Omphi = omb*y_events[0][0][1]/(2*np.pi)     # toroidal precession frequency
    else:
        omb = 2*np.pi/time_events[1][0]
        Omphi = omb*y_events[1][0][1]/(2*np.pi)

    lambda_0 = mu*bmod[0]/H
    r_av = np.mean(abs(zs[0,:]))

    if method != 'thin':
        omb = omb*v0
        Omphi = Omphi*v0

        # If orbit width != 0: classify orbit type according to paper:
        # "Lagrangian neoclassical transport theory applied to the region near
        # the magnetic axis" by Satake et al. (2002)

        thdot = np.zeros(len(time))
        for k in range(len(time)):
            thdot[k] = velo(time[k],zs[:,k])[2]

        sigma_par = zs[4,:]/abs(zs[4,:])
        sigma_theta = thdot/abs(thdot)
        turns_sigma_par = len(np.where(np.diff(np.append(sigma_par,sigma_par[0]))!=0)[0])
        turns_sigma_theta = len(np.where(np.diff(np.append(sigma_theta,sigma_theta[0]))!=0)[0])

        if turns_sigma_par==0 and turns_sigma_theta==0:   orbit_type=0; print('Orbit type: passing')
        elif turns_sigma_par==2 and turns_sigma_theta==2: orbit_type=1; print('Orbit type: banana')
        elif turns_sigma_par==2 and turns_sigma_theta==0: orbit_type=2; print('Orbit type: kidney')
        elif turns_sigma_par==2 and turns_sigma_theta==4: orbit_type=3; print('Orbit type: concave kidney')
        elif turns_sigma_par==0 and turns_sigma_theta==2:
            if np.sign(sigma_par[0]) == 1:                orbit_type=4; print('Orbit type: outer circulating')
            else:                                         orbit_type=5; print('Orbit type: inner circulating')
        else:                                             orbit_type=9; print('No classification possible')
    else:
        orbit_type=9; print('No classification for thin orbits')

    return [omb, Omphi, pphi, H, mu, orbit_type, lambda_0, r_av]


event_bananatip.direction = 1.
event_bananatip.terminal = True
event_circ.terminal = True
nr_bananas = 100


""" Orbit Run """

solution = np.zeros([nr_bananas,8])
init_val = np.zeros([nr_bananas,4])
solution_thin = np.zeros([nr_bananas,8])
plotting = False
th0 = -0.2*np.pi
#r0 = 0.04

for k in range(nr_bananas):
    if k%10 == 0:
        plotting = True
        print(k)

    # Randomize initial conditions
    va = 1.3; vb = 1.5
    v_v0 = (vb - va)*random() + va
    vpa = 0.1; vpb = 0.15
    vpar_v0 = (vpb - vpa)*random() + vpa
    r0a = 0.03; r0b = 0.1
    r0 = (r0b - r0a)*random() + r0a

    # Store [omb, Omphi, pphi, H, mu, orbit_type, lambda_0, r_av] in 'solution'
    solution[k,:] = orbit(r0=r0, th0=th0, v=v_v0, vpar=vpar_v0, plotting=plotting, method='full')
    init_val[k,:] = [r0, th0, v_v0, vpar_v0]
    solution_thin[k,:] = orbit(r0=r0, th0=th0, v=v_v0, vpar=vpar_v0, plotting=False, method='thin')
    plotting = False

plt.show()

#np.savetxt('solution_r_fixed.csv', solution, delimiter=',')
#np.savetxt('solution_thin_r_fixed.csv', solution_thin, delimiter=',')
#np.savetxt('init_val_r_fixed', init_val, delimiter=',')


#""" Magnetic field """
#stest = np.linspace(0,1,20)
#thtest = np.linspace(0,2*np.pi,20)
#ss,tth = np.meshgrid(stest,thtest)
#phi = 0.
#x[2] = 0
#Btest = np.zeros(len(stest))
#
#for k in range(len(stest)):
#    x[0] = stest[k]      # s
#    x[1] = thtest[k]      # theta
#
#    libneo_rt_mc.magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
#
##orbit(plotting=True, vpar=0.9, method='full')
