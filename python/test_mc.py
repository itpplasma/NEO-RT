# %%
import numpy as np
from random import random
import matplotlib.pyplot as plt
from neo_rt_fffi import libneo_rt_mc, parmot_mod
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import Rbf, splrep, splev, sproot, InterpolatedUnivariateSpline
from scipy.integrate import solve_ivp
from IPython import get_ipython


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

print('ro0: {}'.format(parmot_mod.ro0))

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
    """ Passing orbit - full turn at varphi=2pi """
    return z[2] - 2*np.pi

def run_bananas(s):
    # Initial conditions
    z[0] = s
    z[1] = 0.0
    z[2] = -0.7*np.pi
    z[3] = 1.38   # normalized velocity module  v / v_0
    z[4] = 1e-15  # pitch v_\parallel / v, be careful with events !!
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)

    integ = solve_ivp(
        velo, (0, 5e4*dtau), z0, events=(event_bananatip,event_circ), 
        max_step=10*dtau, method='LSODA')

    print(integ.t_events)

    zs = integ.y
    plt.plot(np.sqrt(zs[0,:])*np.cos(zs[2,:]), np.sqrt(zs[0,:])*np.sin(zs[2,:]), ',')
 
    # Om^th = om_b = 2*pi/t_bounce = 2*pi*v0/tau_bounce
    omb = 2*np.pi*v0/integ.t_events[0][0]
    
    # Toroidal precession frequency
    Omphi = omb*integ.y_events[0][0][1]/(2*np.pi)
    
    return [omb, Omphi] 

#parmot_mod.ro0 = rlarm*bmod00  # Actual Larmor radius
#plt.figure(figsize=(4.0, 4.0))
#freqs = []
#srange = np.array([0.005, 0.01, 0.02, 0.03, 0.12, 0.22, 0.32, 0.52, 0.72])
#for sk in srange:
#    freqs.append(run_bananas(sk))
#plt.xlabel(r'$\rho_\mathrm{tor}\,\cos\vartheta_\mathrm{B}$')
#plt.ylabel(r'$\rho_\mathrm{tor}\,\sin\vartheta_\mathrm{B}$')
#plt.title("Drift orbits")
#plt.xlim([-1.0, 1.0])
#plt.ylim([-1.0, 1.0])
#
#parmot_mod.ro0 = 1e-3*rlarm*bmod00  # Smaller Larmor radius
#plt.figure(figsize=(4.0, 4.0))
#freqs_thin = []
#for sk in srange:
#    freqs_thin.append(run_bananas(sk))
#plt.xlabel(r'$\rho_\mathrm{tor}\,\cos\vartheta_\mathrm{B}$')
#plt.ylabel(r'$\rho_\mathrm{tor}\,\sin\vartheta_\mathrm{B}$')
#plt.title("Thin orbit approximation")
#plt.xlim([-1.0, 1.0])
#plt.ylim([-1.0, 1.0])
## %%
#freqs = np.array(freqs)
#freqs_thin = np.array(freqs_thin)
#
#
#plt.figure()
#plt.plot(np.sqrt(srange), freqs_thin[:,0] , 'k--', zorder=10)
#plt.plot(np.sqrt(srange), freqs[:,0], color='tab:red')
#plt.xlabel(r'$\rho_\mathrm{tor}$')
#plt.ylabel(r'$\omega_b / \mathrm{s}^{-1}$')
#plt.title("Bounce frequency")
#plt.xlim([0, 0.9])
#plt.legend(['thin', 'full'])
#
#plt.figure()
## Rescale by factor 1e-3 due to other Larmor radius
#plt.plot(np.sqrt(srange), freqs_thin[:,1]/1e-3, 'k--', zorder=10)
#plt.plot(np.sqrt(srange), freqs[:,1], color='tab:red')
#plt.xlabel(r'$\rho_\mathrm{tor}$')
#plt.ylabel(r'$\Omega_{\mathrm{tB}} / \mathrm{s}^{-1}$')
#plt.title("Toroidal precession frequency (magnetic)")
#plt.xlim([0, 0.9])
#plt.legend(['thin', 'full'])
#


def full_orbit(r=0.7, th=-0.5*np.pi, v=1.38, vpar=0.5, plotting = False):

    # Initial conditions
    z[0] = r        # s
    z[1] = 0.0      # varphi
    z[2] = th0      # vartheta
    z[3] = v        # normalized velocity module  v / v_0
    z[4] = vpar     # pitch v_\parallel / v, be careful with events !!
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64) 
    
    v=v_v0*v0
    vpar = z[4]*v
    vperp = np.sqrt(v**2 - vpar**2)       
        
    # Calculate constants of motion
    libneo_rt_mc.magfie(z[0:3], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    mu = am*p_mass*vperp**2/(2*bmod[0])
    mu = mu/ev
    H = am*p_mass*v**2/2
    H = H/ev
    pphi = 0 
    r0 = z0[0]
        
    # Run FullOrbit:
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64)
    integ = solve_ivp(
        velo, (0, 5e4*dtau), z0, events=(event_bananatip,event_circ), 
        max_step=10*dtau, method='LSODA')
    zs = integ.y
    time_events = integ.t_events
    y_events = integ.y_events
    
    if len(time_events[1])==0:
        # integrate again from tip for a complete orbit
        z02 = zs[:,-1]  
        z02[2] = z02[2]%(2*np.pi)   
        if z02[2] > np.pi: z02[2] = -(2*np.pi - z02[2])     # project angle theta back in (-pi,pi) (as to not wrongly trigger event_circ)
        z02[4] = 1e-15                                      # set vpar to numerical zero at tip
        integ2 = solve_ivp(velo, (0, 5e4*dtau), z02, events=(event_bananatip,event_circ), max_step=1*dtau, method='LSODA')
        time_events = integ2.t_events
        y_events = integ2.y_events
        zs = integ2.y      
    
    orbit_type = 0
    if plotting == True:
        plt.figure()
        plt.plot(np.sqrt(zs[0,:])*np.cos(zs[2,:]), np.sqrt(zs[0,:])*np.sin(zs[2,:]), ',')
        plt.title('Thin Orbits $r_0$ = ' + str(r0) + ', $\\theta_0 = $' + str(z0[2]) + ', v/v0 = ' + str(z0[3]) + ', vpar/v = ' + str(z0[4]))
        plt.show()

    
    if len(time_events[0])>0:
        orbit_type = 1
        omb = 2*np.pi*v0/time_events[0][0]           # Bounce frequency 
        Omphi = omb*y_events[0][0][1]/(2*np.pi)      # Toroidal precession frequency
    
    else:
        omb = 2*np.pi*v0/time_events[1][0]
        Omphi = omb*y_events[1][0][1]/(2*np.pi)
        
    return [omb, Omphi, pphi, H, mu, orbit_type]


   
def zdot_thin(t,y,v,eta):
    ret = np.zeros(5)
    libneo_rt_mc.magfie([y[0],y[1],y[2]], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    ret[0] = 0
    ret[1] = y[4]*v*hctrvr[1]               # varphi
    ret[2] = y[4]*v*hctrvr[2]               # vartheta
    ret[3] = 0  
    ret[4] = -eta*v*hctrvr[2]*bder[2]     # vpar/v
    return ret

def thin_orbit(r=0.7, th=-0.5*np.pi, v=1.38, vpar=0.5, plotting = False):
    
    
    # Initial conditions
    z = libneo_rt_mc._ffi.new('double[5]')
    z[0] = r      # s
    z[1] = 0.0    # varphi
    z[2] = th     # vartheta
    z[3] = v      # normalized velocity module  v / v_0
    z[4] = vpar   # pitch v_\parallel / v, be careful with events !!
    z0 = np.frombuffer(libneo_rt_mc._ffi.buffer(z), dtype=np.float64) 

    v = v_v0*v0
    vpar = z[4]*v
    vperp = np.sqrt(v**2 - vpar**2)       
    eta = vperp**2/(v**2*bmod[0])           # 1/T
    zdot = lambda t, y: zdot_thin(t,y,v,eta)
    
    # Run FullOrbit:
    integ = solve_ivp(zdot, (0, 1e5*dtau/v0), z0, events=(event_bananatip,event_circ), 
        max_step=10*dtau/v0, method='LSODA')
    zs = integ.y
    time_events = integ.t_events
    y_events = integ.y_events
    
    if len(time_events[1])==0:
        # integrate again from tip for a complete orbit
        z02 = zs[:,-1]  
        z02[2] = z02[2]%(2*np.pi)   
        if z02[2] > np.pi: z02[2] = -(2*np.pi - z02[2])     # project angle theta back in (-pi,pi) (as to not wrongly trigger event_circ)
        z02[4] = 1e-15                                      # set vpar to numerical zero at tip
        integ2 = solve_ivp(zdot, (0, 1e5*dtau/v0), z02, events=(event_bananatip,event_circ), max_step=10*dtau/v0, method='LSODA')
        time_events = integ2.t_events
        y_events = integ2.y_events
        zs = integ2.y   

    
    # Calculate constants of motion
    libneo_rt_mc.magfie(z[0:3], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    
    mu = am*p_mass*vperp**2/(2*bmod[0])     # erg/T
    mu = mu/ev                              # eV/T 
    H = am*p_mass*v**2/2                    # erg
    H = H/ev                                # eV  
    pphi = 0
    r0 = z0[0]

    # Plotting 
    if plotting == True:
        plt.figure()
        plt.plot(np.sqrt(zs[0,:])*np.cos(zs[2,:]), np.sqrt(zs[0,:])*np.sin(zs[2,:]), ',')
        plt.title('Thin Orbits $r_0$ = ' + str(r0) + ', $\\theta_0 = $' + str(z0[2]) + ', v/v0 = ' + str(z0[3]) + ', vpar/v = ' + str(z0[4]))
        plt.show()
        
    orbit_type = 0
    if len(time_events[0])>0:
        orbit_type = 1
        omb = 2*np.pi/time_events[0][0]                     # Bounce frequency 
        Omphi = omb*y_events[0][0][1]/(2*np.pi)  # Toroidal precession frequency
    
    else:
        omb = 2*np.pi/time_events[1][0]
        Omphi = omb*y_events[1][0][1]/(2*np.pi)  # Toroidal precession frequency
        

    return [omb, Omphi, pphi, H, mu, orbit_type]



def interpol(solution):
    """
    Interpolate bounce or drift frequency over constants of motion (pph, H)
    for nr random initial values
    """
    x = solution[:,2]       # pphi
    y = solution[:,3]       # H
    z = solution[:,4]       # mu
    w = solution[:,0]       # w_bounce
    
    type_ = solution[:,5]   # orbit type (passing=0/trapped=1)
    idx_pas = np.where(type_ == 0)
    idx_tr = np.where(type_ == 1)
    
    # New values for interpolation
    yi = np.linspace(min(y),max(y),200)
    zi = np.linspace(min(z),max(z),200)
    Yi,Zi = np.meshgrid(yi,zi)
    
    # 2D Rbf interpolation in H and mu
    rbf_lin = Rbf(y,z,w)#/np.sqrt(y))#,function='linear')
    W_rbf_lin = rbf_lin(Yi,Zi)
    
    ### Plotting ###
    
    # 2D plot (H, mu, w_bounce)
    plt.ticklabel_format(style='sci', scilimits=(0,0))
    plt.title('Radial Basis Functions for $\omega_b$ for $r_0$ = ' + str(r0) + ', $th_0$ = ' + str(th0))
    plt.pcolor(Yi,Zi,W_rbf_lin,cmap=cm.jet)
    plt.xlabel('H \ eV')
    plt.ylabel('$\mu$ \ eV/T')
    plt.colorbar()
    plt.show()
    
    # 3D plot (H,mu,w_bounce)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('$H \ eV$')
    ax.set_ylabel('$\mu$ \ eV/T')
    ax.set_zlabel('$\omega_b$ \ Hz')
    ax.ticklabel_format(style='sci', scilimits=(0,0))
    if len(idx_pas[0]) != 0:
        ax.scatter(y[idx_pas], z[idx_pas], w[idx_pas], marker='^', label='Passing Orbits')
    ax.scatter(y[idx_tr], z[idx_tr], w[idx_tr], marker='.', label='Trapped Orbits')
    ax.legend()
    plt.title(plt.title('$r_0$ = ' + str(r0) + ', $th_0$ = ' + str(th0)))
    plt.show()

#    # 3D plot - heatmap (pphi, H, mu, w_bounce)
#    fig = plt.figure()
#    ax = Axes3D(fig)
#    ax.set_xlabel('$p_\phi$')
#    ax.set_ylabel('H \ eV')
#    ax.set_zlabel('$\mu$ \ eV/T')
#    ax.ticklabel_format(style='sci', scilimits=(0,0))
#    if len(idx_pas[0]) != 0:
#        img11 = ax.scatter(x[idx_pas], y[idx_pas], z[idx_pas], c=w[idx_pas], marker='^', cmap='winter', label='Passing Orbits')
#        fig.colorbar(img11)
#    img12 = ax.scatter(x[idx_tr], y[idx_tr], z[idx_tr], c=w[idx_tr],  marker='.',  cmap='autumn', label='Trapped Orbits')
#    fig.colorbar(img12)
#    ax.legend()
#    plt.title(plt.title('$r_0$ = ' + str(r0) + ', $ th_0$ = ' + str(th0)))
#    plt.show()
    
    
    # 2D plot - heatmap (pphi, mu, w_bounce (trapped/passing))
    fig, ax = plt.subplots()
    if len(idx_pas[0]) != 0:
        img = ax.scatter(y[idx_pas], z[idx_pas], c=w[idx_pas], marker='^', cmap='winter', label='Passing Orbits')#/np.sqrt(y[idx_pas]), marker='^', cmap='winter', label='Passing Orbits')
        fig.colorbar(img)
    img2 = ax.scatter(y[idx_tr], z[idx_tr], c=w[idx_tr], marker='.', cmap='autumn', label='Trapped Orbits')#/np.sqrt(y[idx_tr]), marker='.', cmap='autumn', label='Trapped Orbits')
    plt.xlim([min(y) - (max(y)-min(y))*0.1, max(y) + (max(y)-min(y))*0.1])
    plt.ylim([min(z) - (max(z)-min(z))*0.1, max(z) + (max(z)-min(z))*0.1])
    plt.xlabel('H \ eV')
    plt.ylabel('$\mu$ \ eV/T')
    plt.title('$\omega_b$ \ Hz for $r_0$ = ' + str(r0) + ', $\\theta_0$ = ' + str(th0))
    fig.colorbar(img2)
    ax.legend()
    plt.show()
    

event_bananatip.direction = 1.0
event_bananatip.terminal = True
event_circ.terminal = True
nr_bananas = 100
    
# Initial conditions for fixed r0
#    th0a = -np.pi*0.05; th0b = -np.pi*0.5
#    th0 = (th0b - th0a)*random() + th0a
r0 = 0.5
th0 = -0.4*np.pi



""" Orbit Run """
solution = np.zeros([nr_bananas,6])
init_val = np.zeros([nr_bananas,4])
solution_thin = np.zeros([nr_bananas,6])
init_val_thin = np.zeros([nr_bananas,4])
plotting = False
for k in range(nr_bananas):
    if k%20 == 0:
        plotting = True
        print(k)
    va = 1.3; vb = 1.5     
    v_v0 = (vb - va)*random() + va
    vpa = 0.05; vpb = 0.3
    vpar_v0 = (vpb - vpa)*random() + vpa
    
    [omb, Omphi, pphi, H, mu, orbit_type] = full_orbit(r=r0, th=th0, v=v_v0, vpar=vpar_v0, plotting=plotting)
    solution[k,:] = [omb, Omphi, pphi, H, mu, orbit_type]
    init_val[k,:] = [r0, th0, v_v0, vpar_v0]
    
    [omb, Omphi, pphi, H, mu, orbit_type] = thin_orbit(r=r0, th=th0, v=v_v0, vpar=vpar_v0, plotting=plotting)
    solution_thin[k,:] = [omb, Omphi, pphi, H, mu, orbit_type]
    init_val_thin[k,:] = [r0, th0, v_v0, vpar_v0]
    plotting = False
    
    
plt.show()



#get_ipython().run_line_magic('matplotlib', 'qt')
interpol(solution_thin)
interpol(solution)

idx_trapped_full = np.where(solution[:,5]==1)[0]
idx_passing_full = np.where(solution[:,5]==0)[0]
idx_trapped_thin = np.where(solution_thin[:,5]==1)[0]
idx_passing_thin = np.where(solution_thin[:,5]==0)[0]

fig = plt.figure()
plt.plot(solution[idx_trapped_full,4],solution[idx_trapped_full,0],'rx', label='Full Orbits trapped')
plt.plot(solution_thin[idx_trapped_thin,4],solution_thin[idx_trapped_thin,0],'bx', label='Thin Orbits trapped')
plt.plot(solution[idx_passing_full,4],solution[idx_passing_full,0],'r.', label='Full Orbits passing')
plt.plot(solution_thin[idx_passing_thin,4],solution_thin[idx_passing_thin,0],'b.', label='Thin Orbits passing')
plt.title('Bounce frequency $\omega_b$ for full and thin orbits for $r_0$ = ' + str(r0) + ', $\\theta_0$ = ' + str(th0))
plt.xlabel('$\mu$ \ eV/T')
plt.ylabel('$\omega_b$ \ Hz')
plt.legend()
fig.show()

fig2 = plt.figure()
plt.plot(solution[idx_trapped_full,3],solution[idx_trapped_full,0],'rx',label='Full Orbits trapped')
plt.plot(solution_thin[idx_trapped_thin,3],solution_thin[idx_trapped_thin,0],'bx', label='Thin Orbits trapped')
plt.plot(solution[idx_passing_full,3],solution[idx_passing_full,0],'r.',label='Full Orbits passing')
plt.plot(solution_thin[idx_passing_thin,3],solution_thin[idx_passing_thin,0],'b.', label='Thin Orbits passing')
plt.title('Bounce frequency $\omega_b$ for full and thin orbits for $r_0$ = ' + str(r0) + ', $\\theta_0$ = ' + str(th0))
plt.xlabel('H \ eV')
plt.ylabel('$\omega_b$ \ Hz')
plt.legend()
fig2.show()


fig, ax = plt.subplots()
ax.set_xlabel('H \ eV')
ax.set_ylabel('$\mu$ \ eV/T')
ax.ticklabel_format(style='sci', scilimits=(0,0))
ax.scatter(solution[idx_trapped_full,3], solution[idx_trapped_full,4], c=solution[idx_trapped_full,0], marker='x', cmap='Reds', label='trapped')
img=ax.scatter(solution[idx_passing_full,3], solution[idx_passing_full,4], c=solution[idx_passing_full,0], marker='.', cmap='Reds', label='passing')
cbar3 = fig.colorbar(img)
cbar3.set_label('Full Orbits')
ax.scatter(solution_thin[idx_trapped_thin,3],solution_thin[idx_trapped_thin,4],c=solution_thin[idx_trapped_thin,0], marker='x', cmap='Blues')
img2=ax.scatter(solution_thin[idx_passing_thin,3], solution_thin[idx_passing_thin,4], c=solution_thin[idx_passing_thin,0], marker='.', cmap='Blues')
cbar4 = fig.colorbar(img2)
cbar4.set_label('Thin Orbits')
ax.legend()
plt.title('Bounce frequency in Hz for full and thin orbits for $r_0$ = ' + str(r0) + ', $\\theta_0$ = ' + str(th0))
plt.show()

