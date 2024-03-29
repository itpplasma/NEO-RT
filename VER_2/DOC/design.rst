Code design
===========

Bounce averages are required to compute canonical frequencies
and Fourier harmonics. This is implemented in :f:subr:`~orbit/bounce_average`
and relies on a given starting position :code:`zstart`.
A suitable starting position with given invariants :math:`\balpha`
is provided by :code:`TODO`.

1. Read fields and profiles
2. Precompute and interpolate canonical frequencies :math:`\Omega^k(\balpha)`

TODO

Normalization
------------------------
Tracing of guiding-center orbits is implemented in :code:`sub_alpha_lifetime.f90`
in :code:`MC/SRC` as well as :code:`POTATO/VER_X/SRC`

Routines:

* `velo` returns $dz(\bar{t}_s)/d\bar{t}_s$ of phase-variables $z$
* `orbit_timestep` integrates orbits over finite time via RK4/5 with possible random collisions

Normalization is as follows. It differs by the one used internally for symplectic integrators
in SIMPLE by some factors $\sqrt{2}$ (see documentation there).

We introduce thermal velocity

.. math::
  v_{0s}=\sqrt{\frac{2T}{m}},

normalised gyroradius

.. math::
  \rho_{0s}=\frac{mc}{e B_{\mathrm{ref}}}v_{0s},

and cyclotron frequency

.. math::
  \omega_{0s}=\frac{e B_{\mathrm{ref}}}{mc}.

Here the reference temperature $T$ is in units of energy and can also mean the energy of mono-energetic
particles (like fast alphas). The reference $B$ field $B_{\mathrm{ref}}$ is the
physical field in Gauss when field routines give numerical value $B=1$.
So if field routines return Tesla, one should set it to $10^4$, otherwise to $1$.

The actual Larmor radius is given by

.. math::
  \rho &= \frac{v/v_{0s}}{B/B_{\mathrm{ref}}} \rho_{0s} = \frac{\bar v_s}{\bar B_s} \rho_{0s} \\
       &= \frac{\sqrt{H/T}}{B/B_{\mathrm{ref}}} \rho_{0s} = \frac{\sqrt{\bar H_s}}{\bar B_s} \rho_{0s}


Some Python code to set constants properly is as follows:

.. code-block:: py

  import numpy as np

  # Constants
  c = 2.9979e10           # speed of light in cm/s
  qe = 4.8032e-10         # electron charge in franklin ( = 3.336e-10C)
  e_mass = 9.1094e-28     # electron mass in g
  p_mass = 1.6726e-24     # proton mass in g
  ev = 1.6022e-12         # 1 eV in erg ( = 1e-7J)
  am = 2                  # Atomic mass 2 of deuterium ions
  Zb = 1                  # Atomic charge 1 of deuterium ions
  tempi1 = 2e3            # ion temperature in eV
  v0 = np.sqrt(2.0*tempi1*ev/(am*p_mass))  # Reference (thermal) velocity

  # Reference Larmor radius of thermal particles
  m = am*p_mass
  bmod_ref = 1.0                 # What does bmod=1 mean in Gauss?
  ro0 = v0*m*c/(Zb*qe*bmod_ref)  # Larmor radius in bmod_ref

  bmod_test = 1e4                # Field to compute Larmor radius in
  print(f'rlarm = {bmod_ref/bmod_test*ro0:.2f} cm')

Additional variables are

.. math::
  \lambda	=\cos\,\Theta_{p}=\frac{v_{\parallel}}{v}=\frac{v_{\parallel}}{\sqrt{v_{\parallel}^{\,2}+v_{\perp}^{\,2}}}
	  =\sqrt{\frac{m}{2}}\frac{v_{\parallel}}{\sqrt{mv_{\parallel}^{\,2}/2+\mu B}}=\sqrt{\frac{m}{2}}\frac{v_{\parallel}}{\sqrt{H-e\Phi}}.

where the pitch angle $\Theta_{p}$ measures the angle between particle velocity $\boldsymbol{v}$ and magnetic field $\boldsymbol{B}$. Similarly

.. math::
  \lambda^{2}	=\frac{mv_{\parallel}/2}{H-e\Phi}=\frac{H-\mu B-e\Phi}{H-e\Phi}

so

.. math::
  H=\frac{\mu B}{1-\lambda^{2}}+e\Phi.

In the non-relativistic case we have

.. math::
  \bar{v}_{\parallel s}	&=\frac{v_{\parallel}}{v_{0s}}\,,\\
  \bar{H}_{s}	&=\frac{H}{T}\,,\\
  \bar{t}_{s}	&=v_{0s}\,t\,,\\
  \bar{\Phi}_{s}	&=\frac{e}{T}\Phi\,,\\
  \bar{p}_{s}	&=\frac{p}{mv_{0s}}=\frac{v}{v_{0s}}=\bar{v}_{s}\\
  \bar{\mu}_{s}	&=\frac{\bar{p}_{s}^{\,2}(1-\lambda^{2})}{2B}
    =\frac{p^{\,2}(1-\lambda^{2})}{2m^{2}v_{0}^{\,2}B}
    =\frac{p_{\perp}^{\,2}}{4mTB}=\frac{\mu}{2T}\\
  \omega_{c}&=\frac{eB}{mc}=\frac{v_{0s}}{\rho_{0s}}B=\frac{v_{0}}{\rho_{0}}B.

The toroidal canonical momentum $p_{\varphi}$ is normalized to become $\psi^*$,
i.e.  it is multiplied with $c/e$, and has a dimension of poloidal flux in
Gaussian units.

.. math::
  p_{\varphi}
    &=mv_{\parallel}h_{\varphi}+\frac{e}{c}A_{\varphi}\\
    &=mv_{0}\frac{v_{\parallel}}{v}\frac{v}{v_{0}}h_{\varphi}+m\frac{v_{0}}{\rho_{0}}A_{\varphi}\\
    &=\underbrace{\frac{mv_{0}}{\rho_{0}}}_{c/e}
    \underbrace{
      \left(\rho_{0}\frac{v}{v_{0}}\lambda h_{\varphi}+A_{\varphi}\right)
    }_{\bar{p}_\varphi = \psi^*}


Hamiltonian perturbation
------------------------

We compute a toroidal harmonic of the Hamiltonian perturbation according to

.. math::
  H_n = \left(\frac{e_{\alpha}}{m_{\alpha}c}J_{\perp}B_{0}+m_{\alpha}v_{\parallel0}^{2}\right)\frac{B_{n}}{B_{0}}

In :f:subr:`~transport/hpert`, $\bar{H}_n = H_n/T$ is computed with normalization

.. math::
  T	& =\frac{mv_{0}^{\,2}}{2}

  mv_{\parallel}^{\,2}
    &= mv_{0}^{\,2}\frac{v^{2}}{v_{0}^{\,2}}\frac{v_{\parallel}^{\,2}}{v^{\,2}}
     = 2T\frac{v^{2}}{v_{0}^{\,2}}\lambda^{2}

  \mu B	&= \frac{mv_{\perp}^{\,2}}{2}
        = \frac{mv_{0}^{\,2}}{2}\frac{v^{2}}{v_{0}^{\,2}}\frac{v_{\perp}^{\,2}}{v^{2}}
        = T\frac{v^{2}}{v_{0}^{\,2}}(1-\lambda^{2})
