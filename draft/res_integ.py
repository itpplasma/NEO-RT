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
# %%

# %%
