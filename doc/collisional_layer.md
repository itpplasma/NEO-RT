# Collisional boundary-layer closure

`collisional_layer = .true.` enables the local two-boundary closure. Default
is false. Passing/loss coupling is omitted; diffusion and torque use the same
factor. This is a large-aspect-ratio approximation, not a finite-width orbit
solver.

## Reachable distances

The model solves `-G'' + i*z*G = 1`, `G'(-a)=0`, `G(b)=0`, with
`R = Re integral(G)/pi`, `a=(eta_dt-eta)/delta`, `b=(eta-eta_tp)/delta`.
The geometry bounds are `eta_tp=1/Bmax`, `eta_dt=1/Bmin`
([src/magfie.f90:68](../src/magfie.f90#L68)); the trapped search interval is
set by configurable pitch offsets ([src/freq.f90:91](../src/freq.f90#L91)).
The width is `delta=(D_eta/abs(Omega'))**(1/3)`, with
`D_eta=4*nu_d*eps*eta**2*(E(m)/K(m)-(1-m))` and
`m=(eta_dt-eta)/(eta_dt-eta_tp)`
([src/collisional_layer.f90:663](../src/collisional_layer.f90#L663)). Runtime supplies the
root derivative ([src/transport.f90:218](../src/transport.f90#L218)) and
energy-dependent all-species collision rate. Collision densities and
temperatures have positivity/finite checks but no positive lower bound or
finite physical upper bound ([src/collision_rates.f90:47](../src/collision_rates.f90#L47)).
The derivative likewise has no positive lower bound. Thus neither the pitch
search margins nor the velocity cutoff bound the scaled distances uniformly:
physical inputs can drive both toward zero or infinity, and their ratio also
depends on the pitch location. No finite rectangular table covers all inputs.

## Validated evaluation domain

* For `a+b <= 1/5`, use `R=(a+b)**3/(3*pi)`, including either zero-distance
  endpoint. On `x=(z+a)/(a+b)`, the equation is
  `-u''+i*L**3*(x-a/L)*u=1`, `G=L**2*u`, `L=a+b`.
  Dropping the potential gives `u=(1-x**2)/2` and the stated integral.
  The inverse of `-d²/dx²` with these boundary conditions has sup norm `1/2`;
  the real correction starts at second order. Its relative error against the
  leading term is at most `1.5*t**2/(1-t**2)`, `t=L**3/2`, hence below
  `2.401e-5` at `L=1/5`. This bound is uniform in aspect ratio.
* The table includes zero distance on both axes, then `1/16 .. 64` in
  quarter-octaves (42 by 42). Runtime uses tensor cubic Lagrange interpolation
  of `log(R)`. Each coordinate uses `log2(distance)` at or above `1/16`, and
  distance itself between zero and `1/16`. The zero/zero corner is bypassed
  by the narrow-domain formula. Thus asymmetric domains with one arbitrarily
  small distance remain covered without taking `log(0)`.
* A distance beyond 64 uses a matched far-boundary expansion, anchored at
  the table edge. The particular solution is `G=-i/z-2/z**4+...`.
  Integrating its omitted tail gives `2/(3*pi*d**3)`. The reflecting endpoint
  layer contributes `-1/(pi*a**3)`, giving
  `C_a=-1/(3*pi*a**3)-1/(sqrt(2)*pi*a**4.5)+O(a**-6)`;
  the absorbing endpoint layer gives
  `C_b=1/(sqrt(2)*pi*b**1.5)+17/(12*pi*b**3)+O(b**-4.5)`.
  To obtain the next terms, scale the endpoint coordinate by `k=sqrt(+-i*d)`:
  the decaying homogeneous solution is
  `exp(-t)*(1+(t**2+t)/(4*d*k)+...)`. Its integral is
  `(1+3/(4*d*k)+...)/k` and its endpoint derivative is
  `k*(-1+1/(4*d*k)+...)`; combining these with the particular solution
  gives the stated coefficients.
  Runtime adds `C_a(a)-C_a(64)` and/or `C_b(b)-C_b(64)` to the edge value.
  When both distances are at least 64, it uses `1+C_a(a)+C_b(b)` directly.
  This retains the finite opposite boundary or recovers the limit one when
  both distances diverge. It does not hold either distance fixed at 64.
* These branches cover the full nonnegative distance quadrant. Negative and
  nonfinite distances stop explicitly instead of entering interpolation.

## Reproducibility and numerical accuracy

`python3 tools/generate_collisional_layer.py --check` recomputes 1763
nodes by adaptive collocation plus the collapsed zero/zero limit. CTest
`test_collisional_layer_table` independently uses a tridiagonal finite-difference solve with reflected Neumann ghost point,
trapezoidal integration and Richardson extrapolation. Its convergence threshold
is `1e-8` absolute, with a `2e-8` absolute production-node tolerance. All 1681
cell centers and 480 edge/near-axis probes must have less than 1% relative
interpolation error (verified maximum: 0.6601%). These are sampled accuracy
checks, not a proven uniform interpolation bound.

The oracle also checks narrow-domain aspect ratios and endpoints at `1e-4`
relative tolerance, and each far-boundary extension through `(512,256)` at
`2e-8` absolute tolerance (on top of any interpolation error on the finite
axis). It checks thin asymmetric domains and explicit rejection of invalid
inputs in Fast and Debug. At `(0.1,0.1)`, independent collocation and finite
differences give `0.00084882598725`; the narrow-domain formula gives
`0.00084882636316` (relative difference `4.43e-7`).
