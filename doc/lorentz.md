# Bounce-averaged Lorentz response

The `neort_lorentz_native` API solves a complex pitch-diffusion response. It
retains the trapped population and both passing signs with a common physical
phase. It does not replace a Lorentz operator by a scalar resonance width.

Enable this mode in `&params` with `lorentz_pitch_collisions=.true.` and set
`lorentz_pitch_steps` (P1 cells per branch, default 32). `lorentz_rate_source=0`
uses the supplied physical `lorentz_nu` in s^-1. Source 1 requires an explicit
mixture, supplied by `species_file` or the in-memory `neort_prepare_mixture`
API, and sums Maxwellian-background rates at each actual particle speed.
`lorentz_rate_scale` multiplies either source for collisionality scans. These
controls are independent of Krook widths; simultaneous Krook activation is
rejected. `mth_max_abs` applies to the legacy independent-harmonic path; Lorentz
mode always retains the two canonical harmonic records described below.

The supported magnetic well is
`B(theta)=B0*(1-eps*cos(theta)+delta*cos(2*theta))`, with positive `B`,
`eps>4*abs(delta)` and a positive curvature margin resolved above Fourier
roundoff. Its derivative is `B0*sin(theta)*(eps-4*delta*cos(theta))`, so its
only extrema are the nondegenerate minimum at zero and maximum at pi.
`delta=0` preserves the original cosine well. Constant safety factor,
integer `n*q`, thin orbits, linear response and no oscillating drift phase
remain required. The trapped canonical harmonic is zero; both passing
harmonics are `-n*q`. Magnetic shear, multiple wells, finite radial orbit
width, nonlinear attenuation, trapped-only operation and inconsistent
trapped/passing magnetic-drift selections are rejected. These restrictions
certify the common phase and the trapped/passing joining rule.

The standalone format-8 field must contain exactly `(m,n)=(0,0),(1,0)`,
optionally followed by `(2,0)`. This actual-spectrum guard excludes unresolved
harmonics before the analytic single-well certificate is evaluated. Native
orbit integration locates the actual turning points; no cosine turning-angle
formula is used. This finite single-well projection does not admit arbitrary
axisymmetric geometry or noninteger `n*q`.

For physical deflection frequency `nu_D(v)`, the local operator is
`C_L=(nu_D/2)*d_xi[(1-xi^2)*d_xi]`. At fixed speed, the native orbit callback
computes the positive period `tau`, complex Hamiltonian `H`, and
`A=2*nu_D*eta*tau*<1/B-eta>`. The susceptibility solves

```text
-d_eta(A*d_eta z) + tau*i*R*z = tau*H,
R = n*(Omega_E + <Omega_magnetic_drift>).
```

The three pitch intervals share one boundary value at `eta=1/Bmax`. Their weak
form conserves the sum of collisional fluxes at this junction; the physical outer
endpoints have natural zero flux. Interior Gauss quadrature avoids evaluating
the logarithmically singular orbit period at the junction. Native orbit
integration supplies all coefficients, including points outside legacy spline
knots. A pure-diffusion null mode receives no damping floor: compatible sources
use a zero weighted-mean gauge; incompatible sources report no stationary solve.

The nonadiabatic response is
`g=-i*n*c*q/(qi*sign_theta*psi_pr)*(A1+A2*u^2)*z`. This coordinate-consistent
first-order model includes `delta C*F0=-C0[H*d_H F0]`, so a zero physical drive
does not scatter equilibrium. Direct electromagnetic torque uses
`Re(conjg(H)*z)`, with the same causal `1/pi` normalization as the existing
transport kernel. The accompanying D weights are susceptibility-derived legacy
weights; they do not certify a separate particle-flux observable or its relation
to collisional friction.
At zero rotation/drift and zero physical force at every speed, physical torque
is zero while the auxiliary susceptibility may have an incompatible constant
source. Such runs set `transport/susceptibility_valid=false` and explicitly mark
the D columns as undefined placeholders in `<runname>_collision_model.out`.
A nonzero physical source with this incompatibility still fails the solve.
The same marked zero-torque limit applies when the collision rate and detuning
both vanish. Without torque/force initialization, a singular susceptibility is
rejected rather than consulting force values left by an earlier call.

Transport requires a positive physical deflection rate. An enabled zero-rate
calculation with nonzero physical forcing is rejected: a finite pitch graph does
not define the retarded collisionless delta-resonance limit. Use the disabled
Lorentz native path for collisionless resonant calculations. Exact zero forcing
at every speed permits only the marked zero-torque result with undefined D
weights, including when `lorentz_rate_scale=0`.

This is a low-collision orbit-averaged model. It requires collisions slow compared
with bounce/transit motion away from the trapped/passing boundary; it does not
establish an arbitrary-collisionality drift-kinetic approximation. Energy
diffusion and field-particle momentum/energy exchange are absent. Full
test-particle and conserving Fokker–Planck operators require additional velocity
response variables and projections.

The projection follows the transformed conservative operator in
[Brizard et al. (2009), section IV](https://arxiv.org/html/0906.4427) and the
action-space derivation in
[Albert (2020), section 4.9](https://doi.org/10.3217/978-3-85125-746-5).
The response and torque convention follows
[Albert et al. (2016), section III.2](https://arxiv.org/html/1607.04665).
The restricted common-phase graph is the projection implemented here.
