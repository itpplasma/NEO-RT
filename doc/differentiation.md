---
title: Analytical and Automatic Differentiation Plan
---

# Analytical and Automatic Differentiation Plan

This page is the planning document for adding analytical differentiation and
automatic differentiation (AD) support to NEO-RT (issue
[#40](https://github.com/itpplasma/NEO-RT/issues/40)). It defines the layered
implementation order, splits the work into atomic issues, records the
differentiability properties of every input class, and fixes the verification
policy for derivative checks.

This is a plan, not a feature description of a single PR. The atomic issues
below are meant to be implemented one at a time, each with its own tests and
pull request, in the order given. The order starts with the lowest-risk work
(profile algebra) and ends with research-level pieces (separatrix handling,
full reverse-mode transport).

## Goal

Add analytical differentiation and automatic differentiation support to NEO-RT
in layers. The first layer should expose cheap profile and spline sensitivities
for uncertainty quantification and optimizer use. Later layers should handle
resonance-root motion, nonlinear attenuation, orbit averages,
magnetic-geometry sensitivities, and the trapped/passing separatrix.

The near-term user-facing result should be:

- linear error propagation for profile uncertainty,
- gradients for profile fitting and optimizer loops,
- a documented boundary between differentiable sensitivities, piecewise
  sensitivities, and non-differentiable class changes,
- a path toward differentiable Fortran kernels through Enzyme or another AD
  tool once the low-level code is ready.

## Why this belongs in NEO-RT

NEO-RT computes NTV torque in resonant transport regimes by a Hamiltonian
action-angle formulation. The torque and transport coefficients are assembled
from:

- local plasma and rotation profiles,
- local thermodynamic forces,
- canonical bounce/transit and toroidal frequencies,
- roots of the resonance condition,
- orbit averages of the perturbation Hamiltonian,
- optional nonlinear attenuation.

Several inputs affect only local amplitudes. Others move the resonance root or
change the orbit class. These cases have different differentiability properties
and should not be treated by one finite-difference wrapper around the
executable.

## Current code map

### Profile and force layer

- `src/profiles.f90`
  - `prepare_plasma_splines` builds spline coefficients for profile columns.
  - `init_plasma_at_s` evaluates `ni1`, `ni2`, `Ti1`, `Ti2`, `Te` and their `s`
    derivatives, then sets `qi`, `mi`, `vth`, `dvthds`, and collision inputs.
  - `prepare_profile_splines` builds the rotation-profile spline.
  - `init_profile_at_s` evaluates `M_t`, `dM_tds`, `Om_tE = vth*M_t/R0`, and
    `dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0`.
  - `init_thermodynamic_forces` sets `A1` and `A2`.

Density mostly scales amplitudes and torque. Temperature affects `vth`,
`dvthds`, collision coefficients, the Maxwellian factor, and the velocity
normalization. Rotation or electric-field input changes `Om_tE`, hence the
resonance condition.

### Frequency and radial-derivative layer

- `src/freq.f90`
  - `init_canon_freq_trapped_spline` and `init_canon_freq_passing_spline`
    sample eta, call `bounce_time` and `bounce_fast`, and build splines for
    `Omth` and optionally `OmtB`.
  - `Om_th` returns `Omega_theta` and derivatives with respect to `v` and `eta`.
  - `Om_tB` returns the bounce-averaged magnetic precession and derivatives with
    respect to `v` and `eta`.
  - `Om_ph` returns the toroidal canonical frequency. Trapped:
    `Om_tE + OmtB`. Passing: `Om_tE + Omth/iota + OmtB`.
  - `d_Om_ds` currently uses a centered finite difference in `s`, temporarily
    mutating the global flux-surface state.

This layer already exposes derivatives with respect to `v` and `eta`. The
missing part is a reliable derivative with respect to profile and geometry
parameters. `d_Om_ds` is also a direct target for replacement by analytical
differentiation of the local geometry and orbit-average formulas, or by AD once
state is explicit.

### Resonance layer

- `src/resonance.f90`
  - `driftorbit_coarse` scans eta on uniform intervals for sign changes of
    `R(v, eta) = mth*Omth + mph*Omph`.
  - `driftorbit_root` refines one root by bisection and returns `(eta_res,
    dR/deta)`.

The hard-root formula is mathematically the Dirac-delta reduction of the eta
integral. It is piecewise differentiable as long as:

- the number of roots does not change,
- the root is simple,
- the orbit class stays fixed,
- the root does not hit the trapped/passing separatrix or an eta-domain
  boundary.

Inside one such cell, implicit differentiation gives

```text
d eta_res / d p = -(d R / d p) / (d R / d eta)
```

and any root-evaluated quantity `F(v, eta_res, p)` has derivative

```text
d F / d p = F_p + F_eta * d eta_res / d p.
```

The current return value `eta_res(2)` is exactly the denominator needed for
this.

### Transport layer

- `src/transport.f90`
  - `compute_transport_integral` loops over normalized speed `ux`.
  - For each `v = ux*vth`, it scans roots, refines roots, calls `bounce_fast`
    at the root, computes `Hmn2`, applies `nonlinear_attenuation`, and
    accumulates `D11`, `D12`, and torque.
  - It divides the contribution by `abs(eta_res(2))`, the hard-root Jacobian
    from the delta-function eta integral.
  - `D11int`, `D12int`, and `Tphi_int` are pure and cheap.
  - `timestep_transport` duplicates orbit-timestep logic and evaluates the
    perturbation Hamiltonian along the orbit.

This is the main place to add linear error propagation after the low-level
derivative pieces exist. It is also where a soft resonance kernel could be added
as an alternative to hard roots.

### Orbit layer

- `src/orbit.f90`
  - `bounce_time` and `bounce_fast` use DVODE.
  - `bounce_integral` detects an orbit turn and treats passing and trapped
    orbits differently.
  - `vpar` and `vperp` contain square roots of `1 - eta*B` and `eta*B`.
  - orbit class is controlled by `eta < etatp`.

Orbit averages are differentiable only inside a fixed orbit class and away from
turning-point and separatrix singularities. Across the trapped/passing boundary
the bounce time has logarithmic behavior and the class label changes. A
derivative of "the trapped bounce frequency" with respect to a perturbation
that makes the orbit passing is not an ordinary Eulerian derivative. This must
be treated as a one-sided or piecewise derivative, or by a smoothed class
mixture in a deliberately regularized model.

### Nonlinear attenuation layer

- `src/nonlin.f90`
  - `nonlinear_attenuation` calls `Om_ph`, `d_Om_ds`, `omega_prime`,
    `coleff`, and `attenuation_factor`.
  - `omega_prime` is pure and algebraic once its inputs are available.
- `src/attenuation_factor.f90`
  - `attenuation_factor` interpolates `log(Theta)` versus `log(D)`.
- `src/polylag_3.f`
  - `plag1d` and `plag3d` already return interpolation values and derivatives.

The interpolation derivative exists, but `attenuation_factor` currently
discards it. Exposing `dTheta/dD` gives a cheap analytical sensitivity of
nonlinear attenuation once `dD/dp` is known.

### Magnetic-field and libneo boundary

- `CMakeLists.txt`
  - `USE_STANDALONE` defaults to `ON`, using the standalone magnetic-field path.
  - `USE_STANDALONE=OFF` links against an external NEO-2 checkout.
- `src/orbit.f90`, `src/transport.f90`
  - orbit RHS calls `do_magfie`.
  - perturbation amplitude calls `do_magfie_pert_amp`.

This boundary determines whether geometry sensitivities belong inside NEO-RT,
libneo, NEO-2, or an adapter. Profile-only sensitivities can be independent of
libneo. Geometry sensitivities should start at the standalone field path, then
be lifted to libneo once the interface is explicit enough.

## Theory map

### Hamiltonian action-angle foundation

Albert (2020), Chapter 3, derives action-angle variables for guiding-center
motion in axisymmetric tokamak geometry. The radial position on a drift surface
is an implicit function of toroidal momentum through the canonical toroidal
momentum relation. The derivatives of the radial position can be expressed as
ratios of local derivatives of the toroidal momentum. This gives a clean
analytical route for differentiating drift-surface quantities, at least in the
thin-orbit ordering and away from singular boundaries.

The canonical poloidal frequency is the bounce or transit frequency. The
toroidal canonical frequency is the bounce-averaged toroidal velocity. In the
thin-orbit approximation, the toroidal frequency splits into the passing-orbit
transit term and a precession term:

```text
Omega_phi = q*omega_b*delta_tp + Omega_t.
```

The `delta_tp` switch is exactly the trapped/passing class discontinuity
reflected in the code through `eta < etatp`.

### Resonant transport and hard roots

Albert (2020), Chapter 4, writes the resonance condition as

```text
m . Omega - omega = 0.
```

For quasistatic magnetic perturbations in NTV, this becomes

```text
(m2 + n*q*delta_tp)*omega_b + n*(Omega_tE + <Omega_tB>) = 0.
```

The quasilinear delta function collapses the eta integral into a sum over
resonance roots. The transport coefficient contains the Jacobian

```text
1 / |m2*d omega_b/d eta + n*d Omega_phi/d eta|.
```

This is the same object returned by `driftorbit_root`. It is also the source of
large sensitivities near weak or merging roots.

### Nonlinear attenuation

Albert (2020), Section 4.9, replaces the quasilinear result by the same
resonant structure multiplied by a nonlinear attenuation factor `Theta(D)`. The
code implements this in `nonlinear_attenuation` by computing an effective
resonant diffusion coefficient and interpolating a precomputed attenuation
table.

The table lookup itself is differentiable inside table cells because `plag1d`
returns a polynomial derivative. The harder part is the derivative of the
resonant diffusion coefficient, which depends on frequency derivatives,
collision coefficients, orbit averages, and `omega_prime`.

### Literature anchors

- Albert et al. 2016, "Evaluation of toroidal torque by non-resonant magnetic
  perturbations in tokamaks for resonant transport regimes using a Hamiltonian
  approach", arXiv:1607.04665. The NEO-RT scope paper; treats
  low-collisional quasilinear resonant NTV regimes including superbanana
  plateau and drift-orbit resonances in a unified way.
- Albert 2020, "Hamiltonian Theory of Resonant Transport Regimes in Tokamaks
  with Perturbed Axisymmetry", DOI `10.3217/978-3-85125-746-5`. Chapter 3 for
  action-angle variables, Section 3.6 for the thin-orbit trapped/passing
  structure, Chapter 4 for quasilinear resonant transport and nonlinear
  attenuation, Appendix B for comparison to Shaing formulas.
- Shaing et al. 2009a/2009b and Shaing 2015. Used in Appendix B of Albert 2020
  for superbanana plateau, bounce-transit, and drift-resonance comparisons.
- Shaing et al. 2010, "Plasma toroidal momentum dissipation by resonant
  magnetic perturbations in tokamaks", DOI `10.1103/PhysRevLett.105.145002`.
  Relevant for radial-electric-field dependence of NTV torque.
- Kasilov et al. 2014, NEO-2 scope paper for quasilinear NTV comparisons and
  Boozer-coordinate perturbation handling.
- Duarte et al. 2019, "Collisional resonance function in discrete-resonance
  quasilinear plasma systems", DOI `10.1063/1.5129260`. Background for soft or
  collisional resonance broadening beyond a pure Dirac-delta root.
- Burby et al. 2022, "Hamiltonian formulations of quasilinear theory for
  magnetized plasmas", DOI `10.3389/fspas.2022.1010133`. Useful for keeping any
  soft-kernel or adjoint extension in canonical variables.
- Enzyme automatic differentiation documentation, https://enzyme.mit.edu/.
  Enzyme differentiates LLVM IR and documents support across languages that
  lower to LLVM IR, including Fortran. Candidate after NEO-RT kernels are made
  explicit enough for AD.

## Differentiability by input class

### Density profile

Density affects:

- `ni1`, `dni1ds`, `ni2`, `dni2ds` from `init_plasma_at_s`,
- `A1 = dni1ds/ni1 + ...`,
- torque amplitude through `Tphi_int`,
- collision coefficients through `loacol_nbi`.

For fixed `Ti`, `M_t`, magnetic geometry, and fixed resonance roots, density
changes should not move orbit classes or resonance roots except through
collision-dependent nonlinear attenuation. Density is the prime target for
first linear error propagation.

### Temperature profile

Temperature affects:

- `Ti1`, `dTi1ds`, `Ti2`, `Te`,
- `vth = sqrt(2*Ti1*ev/mi)`,
- `dvthds`,
- normalized speed `ux`, physical speed `v = ux*vth`,
- Maxwellian factors in `D11int`, `D12int`, and `Tphi_int`,
- `A1` and `A2`,
- collision coefficients,
- nonlinear attenuation.

Temperature can move roots indirectly because the velocity grid maps through
`vth` and because frequency calls use the physical speed. Tractable before
geometry sensitivities, but should come after density and local-force
derivatives.

### Rotation, electric field, and Mach profile

Rotation affects:

- `M_t`, `dM_tds`,
- `Om_tE`,
- `dOm_tEds`,
- `A1` through the electric precession term,
- `Om_ph`,
- resonance roots.

This is the first important case where the root motion is the main sensitivity.
It should use implicit differentiation of `R = 0`. It is also the input class
where a soft resonance kernel could help optimizer behavior near branch
changes.

### Magnetic geometry and perturbation amplitudes

Geometry affects:

- `Bmin`, `Bmax`, `etatp`, `etadt`, `eps`,
- `B0`, `dVds`, `psi_pr`, `q`, `iota`, `dqds`,
- bounce/transit time,
- magnetic precession,
- perturbation Hamiltonian amplitude,
- orbit class.

Perturbation amplitude affects `Hmn2` directly and in the nonlinear attenuation
parameter. For fixed geometry and roots, this is cheap. Full geometry
differentiation is harder because it reaches `do_magfie`,
`do_magfie_pert_amp`, orbit integration, and libneo/NEO-2 interfaces.

## Work Plan

### Level 0: add a derivative safety harness

Files:

- `test/test_profiles.f90` or new profile-sensitivity test.
- `test/test_reslines.f90`.
- `test/test_omega_prime.f90`.
- Possibly a new Python comparison script next to `test/omega_prime.py`.

Tasks:

- Add golden tests for `init_plasma_at_s`, `init_profile_at_s`, and
  `init_thermodynamic_forces`.
- Add fixed-root tests where the root count and root bracket are asserted before
  any derivative is compared.
- Add finite-difference reference checks only for narrow local functions, not
  the full executable.
- Record root-count changes as events, not derivative failures.

Acceptance:

- `make test` passes.
- A test fails if `Om_tE`, `dOm_tEds`, `A1`, `A2`, or `eta_res(2)` changes
  unexpectedly.
- Tests name the fixed cell: orbit class, harmonic, root number, eta interval,
  and input parameter.

### Level 1: expose analytical profile derivatives

Files:

- `src/profiles.f90`
- new test file under `test/`

Tasks:

- Add a small derived type for local profile state and another for first
  derivatives with respect to selected scalar profile inputs.
- Keep existing global/threadprivate variables for now, but compute derivatives
  in pure helper routines.
- Differentiate:
  - spline value with respect to profile ordinate,
  - spline derivative with respect to profile ordinate,
  - `vth`,
  - `dvthds`,
  - `M_t`,
  - `dM_tds`,
  - `Om_tE`,
  - `dOm_tEds`,
  - `A1`,
  - `A2`.
- Include direct derivatives with respect to local values first. Add
  spline-control-point derivatives only after the local algebra is tested.

Notes:

- For `vth = sqrt(2*Ti1*ev/mi)`, `d vth / d Ti1 = vth/(2*Ti1)`.
- For `Om_tE = vth*M_t/R0`,
  `d Om_tE = M_t/R0*d vth + vth/R0*d M_t`.
- For `A2 = dTi1ds/Ti1`,
  `d A2 = d(dTi1ds)/Ti1 - dTi1ds*dTi1/Ti1^2`.
- For `A1`, include the density-gradient quotient, the electric-precession
  term, and the temperature-gradient quotient separately.

Acceptance:

- Analytical profile derivatives agree with centered finite differences for
  local test inputs.
- Density-only derivatives show no root motion when nonlinear attenuation is
  disabled and resonance brackets are held fixed.

### Level 2: expose interpolation derivatives where they already exist

Files:

- `src/attenuation_factor.f90`
- `src/polylag_3.f`
- tests under `test/`

Tasks:

- Add an API returning both `Theta` and `dTheta/dD`.
- Use the existing `plag1d` derivative:

```text
log(Theta) = fun(log(D))
dTheta/dD = Theta * der / D
```

- Add bounds behavior tests for the attenuation table.
- Do not change the default `attenuation_factor(D, Theta)` API until the
  derivative API is tested.

Acceptance:

- The derivative API matches finite differences inside table cells.
- Existing nonlinear attenuation behavior is unchanged.

### Level 3: hard-root implicit differentiation

Files:

- `src/resonance.f90`
- `src/freq.f90`
- `src/transport.f90`
- tests under `test/`

Tasks:

- Add a `resonance_root_t` style result or parallel arrays containing:
  - `eta`,
  - `dR/deta`,
  - `dR/dv`,
  - root bracket,
  - orbit class,
  - root status.
- Add a derivative API for a selected parameter `p`:

```text
deta_dp = -dR_dp / dR_deta
```

- Start with `p = Om_tE` because `dR/dOm_tE = mph` for both trapped and passing
  branches.
- Then add `p = M_t`, `p = Ti1`, and local profile-state parameters through
  Level 1.
- Keep a hard guard for `abs(dR/deta)` below a configured threshold. Return a
  non-differentiable or ill-conditioned status instead of a large silent
  number.

Acceptance:

- For fixed simple roots, implicit `d eta_res / d Om_tE` matches finite
  differences that keep the same root branch.
- The code reports a structured event when the root count changes or a root hits
  a boundary.

### Level 4: differentiate hard-root transport contributions inside one cell

Files:

- `src/transport.f90`
- `src/nonlin.f90`
- tests under `test/`

Tasks:

- Split one root contribution into a pure-ish helper:

```text
root_contribution(ux, eta_res, dR_deta, state) -> D11, D12, T, diagnostics
```

- Add derivative propagation:

```text
d/dp [F(eta_res,p) / |R_eta|]
```

  inside fixed-sign, fixed-root cells.
- Include derivatives of `D11int`, `D12int`, `Tphi_int`, `Hmn2`, and attenuation
  as they become available.
- Keep root-count and class-change events out of the linear derivative. Report
  them separately for UQ.

Acceptance:

- One selected benchmark root has a derivative that agrees with
  branch-preserving finite differences.
- If the finite-difference perturbation changes root count, the test expects an
  event instead of a derivative comparison.

### Level 5: add a soft resonance option

Files:

- `src/resonance.f90`
- `src/transport.f90`
- tests and one diagnostic plot script

Motivation:

The hard-root form is cheap and correct for simple isolated roots. It is brittle
for optimizers because root counts can change and `1/|R_eta|` blows up near weak
roots. A soft kernel keeps the eta integral explicit:

```text
integral d eta F(v, eta, p) K_epsilon(R(v, eta, p))
```

with `K_epsilon` normalized to approximate `delta(R)`.

Candidate kernels:

- Lorentzian:

```text
K_epsilon(R) = (1/pi) * epsilon / (R^2 + epsilon^2)
```

- Gaussian:

```text
K_epsilon(R) = exp(-R^2/(2*sigma^2)) / (sqrt(2*pi)*sigma)
```

- Collisional broadening kernel based on the resonance-function literature, if
  a physical width is available.

Cost:

The concern that the integral along the resonance may remain cheap is valid.
NEO-RT already builds eta splines for frequencies and scans eta intervals. A
soft mode can evaluate a modest eta quadrature per `ux`, harmonic, and orbit
class. It avoids repeated bisection and branch bookkeeping. The expensive part
is not evaluating the scalar kernel; it is obtaining orbit averages. Therefore
the soft mode should reuse the existing frequency and `Hmn` splines or add a
cache for `bounce_fast` results on the eta grid.

Design:

- Add a runtime option:

```text
resonance_mode = "hard_root" | "soft_eta"
resonance_width = ...
```

- In `soft_eta`, integrate over eta using cached/splined `Omth`, `Omph`,
  `Hmn2`, and attenuation inputs.
- Keep `hard_root` as the default.
- Expose derivative formulas:

```text
d/dp K(R) = K'(R) * dR/dp
```

  and include `F_p`.

Acceptance:

- As width tends to zero, `soft_eta` approaches `hard_root` for an isolated
  simple root.
- For a root-count-changing scan, `soft_eta` produces a continuous response
  while `hard_root` reports events.
- Runtime is measured against the hard-root path.

### Level 6: orbit-average sensitivities within a fixed class

Files:

- `src/orbit.f90`
- `src/freq.f90`
- `src/transport.f90`

Tasks:

- Add a fixed-class orbit-average derivative path.
- Start with analytical derivatives of algebraic RHS pieces in
  `poloidal_velocity`.
- Then differentiate the ODE solution with variational equations or forward AD.
- Keep event-time derivatives explicit: the bounce time is an event/root of the
  orbit map.
- Do not cross `eta = etatp` in this level.

Acceptance:

- `dOmth/deta` and `dOmth/dv` from the new path agree with the existing spline
  derivatives away from the separatrix.
- One orbit-average derivative of `Hmn` agrees with a branch-preserving finite
  difference.

### Level 7: replace `d_Om_ds` finite differencing

Files:

- `src/freq.f90`
- `src/orbit.f90`
- magnetic-field interface modules

Tasks:

- Replace the current centered `s` finite difference with:
  - analytical local derivatives where available,
  - variational orbit sensitivities for orbit averages,
  - or AD-generated derivatives of a pure kernel.
- Stop mutating global `s` inside derivative routines. Pass state explicitly,
  at least through the derivative path.
- Add a test that catches state leakage across OpenMP threads.

Acceptance:

- New `d_Om_ds` agrees with the old finite difference on a non-singular point.
- Repeated calls in different order give identical results.
- The derivative path does not change global flux-surface state after return.

### Level 8: AD preparation for Fortran kernels

Files:

- `src/profiles.f90`
- `src/freq.f90`
- `src/orbit.f90`
- `src/transport.f90`
- build files

Tasks:

- Isolate pure kernels with explicit inputs and outputs.
- Remove hidden dependencies on mutable module state inside candidate AD
  kernels.
- Avoid AD through DVODE at first. Differentiate small kernels and
  interpolation routines first.
- Add a build experiment with an LLVM Fortran route and Enzyme only after the
  pure kernels are separated.
- Document unsupported calls: I/O, global state mutation, external solver
  internals, and threadprivate state.

Candidate first AD kernels:

- profile algebra,
- `D11int`, `D12int`, `Tphi_int`,
- `omega_prime`,
- attenuation interpolation,
- local magnetic-field algebra if inputs are already explicit.

Non-candidate first kernels:

- `bounce_fast` through DVODE,
- `bounce_integral` event detection,
- `driftorbit_root` branch logic,
- `do_magfie` until the interface contract is explicit.

Acceptance:

- A small kernel has an AD derivative checked against the analytical derivative
  and finite difference.
- The build path is optional and does not break the default compiler path.

### Level 9: trapped/passing separatrix handling

Files:

- `src/orbit.f90`
- `src/freq.f90`
- `src/resonance.f90`
- `src/transport.f90`

Problem:

The derivative of a trapped-orbit quantity is not defined as an ordinary smooth
derivative when the perturbed input changes the orbit into a passing one. The
code branches on `eta < etatp`, and the theory contains `delta_tp`. Near the
boundary, bounce/transit times and derivatives are singular.

Options:

- Treat sensitivities as piecewise and one-sided. This is physically honest and
  best for UQ reports.
- Add an event measure: probability or interval mass crossing the separatrix
  under input uncertainty.
- Add a smooth class mask only in a deliberately regularized model, not as a
  silent replacement for the hard Hamiltonian model.
- In soft-resonance mode, integrate over eta with a class-aware kernel and
  explicit exclusion or smoothing around `etatp`.

Acceptance:

- The code reports separatrix-crossing events in sensitivity runs.
- Documentation states that class-changing perturbations invalidate linear
  error propagation for class-specific quantities.
- Any smoothed class mode has an option name that marks it as regularized.

### Level 10: full adjoint or reverse-mode transport

Files:

- all major physics kernels, after lower levels are done.

Tasks:

- Define scalar objectives:
  - total torque on a surface,
  - flux-surface scan mismatch,
  - optimizer loss against experiment or another code.
- Build adjoints only after:
  - profile derivatives are tested,
  - root events are explicit,
  - orbit-average derivative policy is fixed,
  - geometry interfaces are explicit.
- Prefer reverse mode only for many input parameters and few objectives. Forward
  mode is simpler for a small number of uncertain profile parameters.

Acceptance:

- One scalar objective has a gradient with respect to many profile control
  points.
- Gradient checks pass on a fixed-root, fixed-class case.
- Non-differentiable events are reported beside the gradient.

## UQ policy

Linear error propagation is valid only inside a fixed differentiability cell:

```text
same orbit class
same root count
same root ordering
simple roots
no eta-boundary crossing
no separatrix crossing
same table cell or smooth interpolation path
```

The sensitivity output should therefore contain both:

- derivative values for the fixed cell,
- event diagnostics for root creation, root annihilation, boundary hit, and
  class change.

For profile uncertainty, start with independent scalar perturbations at one flux
surface. Then move to spline-control-point covariance:

```text
Var[y] = J C J^T
```

where `J` is the sensitivity of the selected output to profile control values
and `C` is the profile covariance.

## Atomic issues

1. `[ad-01] add profile-state derivative tests`
2. `[ad-02] expose analytical derivatives for profile algebra`
3. `[ad-03] expose dTheta/dD from attenuation interpolation`
4. `[ad-04] return structured resonance-root diagnostics`
5. `[ad-05] add implicit derivative of eta_res with respect to Om_tE`
6. `[ad-06] differentiate one hard-root transport contribution`
7. `[ad-07] add soft_eta resonance mode behind an option`
8. `[ad-08] replace d_Om_ds finite difference in a fixed-class case`
9. `[ad-09] prepare first explicit Fortran kernel for Enzyme experiment`
10. `[ad-10] report root and separatrix events in UQ output`

## Non-goals for the first PRs

- Do not differentiate the whole executable by finite differences and call it
  AD.
- Do not run AD through DVODE as the first target.
- Do not hide root-count changes behind a large derivative.
- Do not smooth the trapped/passing boundary without an explicit option and
  documentation.
- Do not make libneo a hard dependency for profile-only sensitivities.

## Verification plan

Every implementation PR should include:

```bash
make test
```

and a focused derivative check. The derivative check should name:

- input parameter,
- orbit class,
- harmonic indices,
- root bracket,
- root number,
- hard or soft resonance mode,
- whether nonlinear attenuation is enabled.

Finite-difference checks must be branch-preserving. If the finite-difference
step changes root count or orbit class, the expected result is a diagnostic
event, not a failed derivative.
