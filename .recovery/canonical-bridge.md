# Raw full-orbit canonical bridge

This packet extends the reviewed finite-drive contract, not its terminal
transport claim. Base: NEO-RT `8c429a5fb2f5182295843fb9e31b6b909ab89e92`.
The analytical S3 work remains a separate inherited candidate.

For a static axisymmetric retained-order guiding-center model, use physical
cylindrical components and the one-form

`Gamma = [(q/c) A0 + m vpar b] . dX + (m c/q) mu d(zeta)`.

SI replaces `q/c` by `q` and `m c/q` by `m/q`. The Hamiltonian is
`E = m vpar^2/2 + mu |B0| + q Phi0`. Background A0 must produce the same B0;
providing an unrelated vector potential invalidates all action claims.

The conserved canonical toroidal momentum is
`Pphi = R [(q/c) A0_phi + m vpar b_phi]`. At fixed Pphi and mu, reducing the
cyclic toroidal coordinate gives the meridional one-form
`[(q/c) A0_R + m vpar b_R] dR + [(q/c) A0_Z + m vpar b_Z] dZ`.
Its primitive-cycle integral divided by 2pi defines Jb. Trapped cycles follow
physical time. Passing cycles use a fixed positive poloidal orientation:
multiply the time integral by the physical winding sigma=+/-1. Accordingly,
`Omega_b = sigma 2pi/T` and `Omega_phi = Delta(phi)/T`; trapped sigma is +1.
Do not mix a signed passing frequency with a time-oriented unsigned action.

For smooth nondegenerate primitive cycles at fixed Pphi and mu,
`dJb/dE = 1/Omega_b` and `dJb/dPphi = -Omega_phi/Omega_b`.
Numerical differentiation must shoot initial conditions at fixed canonical
invariants. Varying radius or pitch independently does not implement these
derivatives. The separate S5 worker prepares this independent oracle.

The six-dimensional symplectic volume is
`m^2 R |b.Bstar| dR dphi dZ dvpar dmu dzeta`, with
`Bstar = B0 + (m c/q) vpar curl(b)` (SI c->1). Equivalently, it is canonical
`dJb dPphi dJg dtheta_b dtheta_phi dzeta`, where
`Jg = (m c/|q|) mu` under the positive gyro-cycle orientation.
Integrating all three angles gives
`(2pi)^2 (m c/|q|) T dE dmu dPphi`. A velocity-space distribution must be
converted to canonical phase density by `f_canonical = f_velocity/m^3`;
otherwise using this measure silently changes the physical particle density.
The independent S5 check compares the full 6x6 one-form determinant to this
density, rather than checking source expressions against themselves.

The first implementation returns complete radial trajectories, the same Xdot
for the drive, energy/Pphi, primitive-cycle action and local symplectic density.
It refuses actions if A0 is absent or the declared endpoint fails closure.
The caller must independently establish that the supplied period is primitive.
The initial oracle uses an independently solved closed trapped orbit through
q=2, and an actual nonzero rational magnetic mode with paired electrostatic
baseline. A0=(0,psi/R,-log R) produces the circular fixture B0 exactly.

Remaining physical transport bridge: find actual FOW resonant roots, use their
canonical action derivatives/Jacobian and physical distribution, and specify a
consistent collisional operator. A fixed Lorentzian applied to off-resonance
bare H is only a fixed-gauge diagnostic: H changes by detuning times chi and
the collision/distribution transformation has not been supplied. On-resonance
delta-function transport avoids that bare-H ambiguity but still requires the
root and population measure. The omitted delta-b symplectic terms remain a
higher-order arbitrary-field limitation of the declared retained-order model.

The thin-island terminal limit additionally needs a true paired raw/ideal
physical projection, the actual island population mask, and quantitative
comparison of w to the radial FOW, resonance and collisional widths. A bounded
finite-orbit H or epsilon~w^2 scan does not prove a uniform physical D/torque
limit; a resonance layer aligned with the island can falsify uniformity.
