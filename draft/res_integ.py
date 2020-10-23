# See https://doc.lagout.org/science/0_Computer%20Science/2_Algorithms/Implicit%20Curves%20and%20Surfaces_%20Mathematics%2C%20Data%20Structures%2C%20and%20Algorithms%20%5BGomes%2C%20Voiculescu%2C%20Jorge%2C%20Wyvill%20%26%20Galbraith%202009-05-15%5D.pdf
# https://hplgit.github.io/prog4comp/doc/pub/._p4c-solarized-Matlab016.html
# https://www.sciencedirect.com/science/article/pii/S0021999113002337?via%3Dihub
# https://arxiv.org/pdf/1706.00840.pdf

#%%
from numpy import *
from matplotlib.pyplot import *

def f(x, y): return x**2 + y**2 - 1.0
def gradf(x, y): return array((2.0*x, 2.0*y))

x = linspace(-1.2, 1.2, 30)
y = linspace(-1.2, 1.2, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

figure()
contour(X, Y, Z)
colorbar()

def find_coarse(f, smin, smax, tol, maxit):
    '''
    Coarse root search by bisection
    '''
    z1 = f(smin)
    z2 = f(smax)
    print(z1, z2)
    if z2 < z1: smin, smax = smax, smin  # monotonicity
    for k in range(maxit):
        s = (smin+smax)/2.0
        z = f(s)
        print(f'Iteration {k}: [{smin}, {smax}], s={s}, z={z}')
        if abs(z) < tol: break
        if z > 0.0:
            smax = s
        else:
            smin = s

    return s

find_coarse(lambda s: f(0.2, s), smin=0.5, smax=1.0, tol=1e-3, maxit=100)

#%% Explicit Euler
dt = 0.1
nt = 100

w = empty((nt+1, 2))
w[0, :] = (1.0, 0.0)

for k in range(nt):
    df = gradf(w[k, 0], w[k, 1])
    jac = sqrt(sum(df**2))
    w[k+1, 0] = w[k, 0] - dt*df[1]/jac
    w[k+1, 1] = w[k, 1] + dt*df[0]/jac

plot(w[:, 0], w[:, 1], 'rx')

#%% Implicit midpoint rule
from scipy.optimize import root

def rootfun(w, wprev):
    wmid = (w + wprev)/2.0
    df = gradf(wmid[0], wmid[1])
    jac = sqrt(sum(df**2))
    return np.array((
        w[0] - (wprev[0] - dt*df[1]/jac),
        w[1] - (wprev[1] + dt*df[0]/jac),
        w[2] - (wprev[2] + 1.0*dt/jac)
    ))

w = empty((nt+1, 3))
w[0, :] = (1.0, 0.0, 0.0)
for k in range(nt):
    sol = root(rootfun, x0=w[k, :], args=w[k, :], tol=1e-13)
    w[k+1, :] = sol.x

plot(w[:, 0], w[:, 1], 'ks')
legend(('Expl Euler', 'Midpoint'))

# %%
figure()
t = linspace(0, 100*dt, 101)
plot(w[:, 0]**2 + w[:, 1]**2)
plot((2*pi*w[:, 2]*dt/t)**2)

# %%
jac0 = 2.0
jac0*w[:, 2]/t
#(w[:, 2]/t)/(sqrt(w[:, 0]**2 + (w[:, 1]-0.5)**2))
# %% See https://plotly.com/python/3d-isosurface-plots/
from plotly.graph_objects import *

def f(x, y, z): return x**2 + y**2 + z**2 - 1.0
def gradf(x, y, z): return array((2.0*x, 2.0*y, 2.0*z))

X, Y, Z = mgrid[-1.2:1.2:40j, -1.2:1.2:40j, -1.2:1.2:40j]

F = f(X, Y, Z)

Figure(data = Isosurface(
    x = X.flatten(), y = Y.flatten(), z = Z.flatten(), value=F.flatten(),
    isomin = 0.0, isomax = 0.0, surface_count = 1,
    colorbar_nticks=5, # colorbar ticks correspond to isosurface values
    caps=dict(x_show=False, y_show=False)
))


#%% Implicit midpoint rule

def rootfun(w, wprev):
    wmid = (w[:3] + w[3:])/2.0
    df = gradf(wmid[0], wmid[1], wmid[2])
    dwdu = (w[:3] - wprev)/du
    dwdv = (w[3:] - wprev)/dv
    return np.array((
        dwdu[1]*dwdv[2] - dwdu[2]*dwdv[1] - df[0],
        dwdu[2]*dwdv[0] - dwdu[0]*dwdv[2] - df[1],
        dwdu[0]*dwdv[1] - dwdu[1]*dwdv[0] - df[2],
        #dwdu[0]**2 + dwdu[1]**2 + dwdu[2]**2,
        #dwdv[0]**2 + dwdv[1]**2 + dwdv[2]**2,
        dwdu[0],
        dwdu[1],
        dwdu[0]*dwdv[0] + dwdu[1]*dwdv[1] + dwdu[2]*dwdv[2],
        #dwdu[0]**2 + dwdu[1]**2 + dwdu[2]**2,
        #dwdv[0]**2 + dwdv[1]**2 + dwdv[2]**2,
        #0.0
    ))

du = 0.01
dv = 0.01
nu = 9
nv = 1

w = empty((nu+1, nv+1, 3))
w[0, 0, :] = (1.0, 0.0, 0.0)

for ku in range(nu):
    for kv in range(nv):
        sol = root(
            rootfun,
            x0=np.concatenate((w[ku, kv, :], w[ku, kv, :])),
            args=w[ku, kv, :], tol=1e-13
        )
        print(sol)
        w[ku+1, kv, :] = sol.x[:3]
        w[ku, kv+1, :] = sol.x[3:]

from plotly.express import *

scatter_3d(
    x = w[:,:,0].flatten(),
    y = w[:,:,1].flatten(),
    z = w[:,:,2].flatten()
)

#%% Implicit midpoint rule

def rootfun(w, wprev):
    wmid = (w[:3] + w[3:])/2.0
    df = gradf(wmid[0], wmid[1], wmid[2])
    dwdu = (w[:3] - wprev)
    dwdv = (w[3:] - wprev)
    return np.array((
        dwdu[1]*dwdv[2] - dwdu[2]*dwdv[1] - df[0]*du*dv,
        dwdu[2]*dwdv[0] - dwdu[0]*dwdv[2] - df[1]*du*dv,
        dwdu[0]*dwdv[1] - dwdu[1]*dwdv[0] - df[2]*du*dv,
        dwdu[0]**2 + dwdu[1]**2 + dwdu[2]**2,
        dwdv[0]**2 + dwdv[1]**2 + dwdv[2]**2,
        0.0
        #dwdu[0]
        #dwdu[2] - 1.0,
        #dwdu[0]*dwdv[0] + dwdu[1]*dwdv[1] + dwdu[2]*dwdv[2],
        #dwdu[0]**2 + dwdu[1]**2 + dwdu[2]**2,
        #dwdv[0]**2 + dwdv[1]**2 + dwdv[2]**2,
        #0.0
    ))

du = 0.01
dv = 0.01
nu = 1
nv = 1

w = empty((nu+1, nv+1, 3))
w[0, 0, :] = (1 .0, 0.0, 0.0)

for ku in range(nu):
    for kv in range(nv):
        sol = root(
            rootfun,
            x0=array((1, 0.1, 0, 1, 0, 0.1)),
            #x0=np.concatenate((w[ku, kv, :], w[ku, kv, :])),
            args=w[ku, kv, :], tol=1e-13
        )
        print(sol)
        w[ku+1, kv, :] = sol.x[:3]
        w[ku, kv+1, :] = sol.x[3:]

from plotly.express import *

scatter_3d(
    x = w[:,:,0].flatten(),
    y = w[:,:,1].flatten(),
    z = w[:,:,2].flatten()
)


#%% Step one: find curve on surface by cutting with coordinate plane

x30 = 0.0  # initial level for coordinate plane
du = 2*pi/40
dv = 0.1
nu = 20
nv = 20

def rootfun(w, wprev):
    wmid = (w + wprev)/2.0
    df = gradf(wmid[0], wmid[1], x30)
    jac = sqrt(sum(df**2))
    return np.array((
        w[0] - wprev[0] + du*df[1]/jac,
        w[1] - wprev[1] - du*df[0]/jac
    ))

w = empty((nu+1, 2))
H = empty((nu+1))
H[0] = 0
w[0, :] = (1.0, 0.0)

for ku in range(nu):
    sol = root(
        rootfun,
        x0=w[ku, :],
        args=w[ku, :], tol=1e-13
    )
    w[ku+1, :] = sol.x
    H[ku+1] = H[ku] + du  # To check if we integrate arc length correctly

from plotly.express import *

scatter_3d(
    x = w[:,0].flatten(),
    y = w[:,1].flatten(),
    z = x30*ones(nu+1)
)

# %%
wmid = (w[1:,:] + w[:-1,:])/2.0
eu = (w[1:,:] - w[:-1,:])/du

def rootfun(w, wprev):
    wmid = (w + wprev)/2.0
    df = gradf(wmid[0], wmid[1], wmid[2])
    jac = sqrt(sum(df**2))
    return np.array((
        eu[ku, 1]*(w[2] - wprev[2]) - dv*df[0]/jac,
        -eu[ku, 0]*(w[2] - wprev[2]) - dv*df[1]/jac,
        eu[ku, 0]*(w[1] - wprev[1]) - eu[ku, 1]*(w[0] - wprev[0]) - dv*df[2]/jac
    ))

ku = 1
w2 = empty((nv+1, 3))
w2[0, :] = np.append(wmid[ku], x30)

for kv in range(nv):
    sol = root(
        rootfun,
        x0=w2[kv, :],
        args=w2[kv, :], tol=1e-13
    )
    w2[kv+1, :] = sol.x

Figure(
    (
        Isosurface(
            x = X.flatten(), y = Y.flatten(), z = Z.flatten(), value=F.flatten(),
            isomin = 0.0, isomax = 0.0, surface_count = 1,
            colorbar_nticks=5, # colorbar ticks correspond to isosurface values
            caps=dict(x_show=False, y_show=False)
        ),
        Scatter3d(
            x = w[:,0].flatten(),
            y = w[:,1].flatten(),
            z = x30*ones(nu+1)
        ),
        Scatter3d(
            x = w2[:,0].flatten(),
            y = w2[:,1].flatten(),
            z = w2[:,2].flatten()
        ),
    )
)




# %%
