---
title: Linear Krook collisions
---

# Linear Krook collisions

Set `krook_collisions=.true.` and supply a positive `krook_nu` to broaden the
linear resonance. `krook_nu` is a damping rate in s^-1 in the same time units as
the orbit frequencies. Apply no factor of 2π. The default is disabled, with a
zero rate. Disabled collisions or an exactly zero rate use the existing
collisionless transport routine. Positive-width Krook transport requires
`nonlin=.false.`; the nonlinear attenuation model has a separate physical
derivation.

For each normalized speed `u=v/vth` and orbit harmonic, NEO-RT replaces

\[
\sum_{\eta_r:R(\eta_r)=0}\frac{F(u,\eta_r)}{|\partial_\eta R|}
\quad\hbox{by}\quad
\int_{\eta_{\min}}^{\eta_{\max}}F(u,\eta)
\frac{\nu(u)}{\pi[R(u,\eta)^2+\nu(u)^2]}\,d\eta,
\qquad R=m_\theta\Omega_\theta+n\Omega_\phi.
\]

`F` is the existing diffusion or torque integrand, including its bounce period,
orbit Fourier amplitude, Maxwellian and thermodynamic forces. Its prefactor
already includes the delta-resonance normalization, so the Lorentzian contains
exactly one factor `1/pi`. The native pitch coordinate `eta` has units G^-1;
the integration uses `deta` without a normalized-pitch Jacobian. Orbits away
from resonance retain the periodic orbit phase used by the collisionless
Fourier coefficient. Finite positive widths keep nonresonant tails.

## Rate selection and matching

The supplied physical rate is independent of NEO-RT's legacy nonlinear
collision-frequency initialization. For `g=krook_geometry_trapped` or
`krook_geometry_passing`, as appropriate, the models are:

| `krook_model` | Consumed damping rate |
| --- | --- |
| `0` | `krook_nu*g` |
| `1` | `krook_nu*g*u**(-3)*(1+ell_eff**2/4)` |
| `2` | `krook_nu*g*u**(-3)` |

`krook_harmonic=0` uses integer `ell_eff=mth`. Setting it to `1` uses
`ell_eff=mth+n*q` for both passing directions and `ell_eff=mth` for trapped
orbits. The sign of the transit frequency distinguishes the passing
directions. Each geometry factor is positive and defaults to one; no geometry
normalization is inferred from the input rate.

The source contracts used for these choices are:

- MARS-Q commit `e6329433cb3facf343a955368e623c6d0ea4eeda`:
  `MarsQ_2FK/kinetic.f:3749-3767` assigns the direct rate (`NUMODEL=0`) or the
  rate divided by `2*eps` (`NUMODEL=1`). `anisotropic.f:4740-4756` gives the
  constant, energy/harmonic and energy-only denominators. The fixed passing
  paths set `H4=1` at lines 4441-4448 and 4631-4638. Use integer harmonics and
  `g=1` with a supplied rate that already includes any `NUMODEL` scaling.
- PENTRC commit `76723ee33fce4b47a5041d27cba858942360012b`:
  `pentrc/energy.f90:375-399` gives the constant and energy/harmonic models;
  `pentrc/pitch.f90:357-369` uses `ell+n*q` and the input rate for passing,
  and `ell` with input rate divided by `2*eps` for trapped. Matching the raw
  PENTRC input rate therefore uses canonical harmonics,
  `krook_geometry_passing=1` and `krook_geometry_trapped=1/(2*eps)`.

The latest fixed MARS passing detuning omits mean magnetic drift. Its matching
lane uses `magdrift_passing=0` and `drift_phase=.false.`. Matching damping rates
alone does not establish agreement between different orbit or response models.
HDF5 `/config` records the resolved collision controls and species selector;
`qs`, `ms`, `vth` and `Om_tE` retain their existing output locations.

## Integration and verification

`krook_eta_steps` (default 32) sets the coarse derivative probes. Each cell also
probes the interior extremum of its cubic-Hermite derivative proxy, so two
folds can be found inside a cell whose endpoint derivatives have the same
sign. Located stationary points split the residual into monotone branches.
The proxy coefficients and pointwise Lorentzian, rate and mapped-node leaves
are emitted by FortSym from `gen/gen_krook.f90`; their checked-in output is
`src/generated/krook_leaf.inc`.

For narrow widths the integrator uses the residual angle
`theta=atan(R/nu)`. The angle difference is evaluated with scaled `atan2`
operands, retaining tiny same-sign tail spans. Geometric residual bands and a
relative sine/cosine map avoid losing tails or narrow resonances through
absolute-angle rounding. A sine-squared angle map regularizes simple quadratic
fold endpoints. When the residual variation is at most `nu`, a direct pitch
rule integrates the same finite Lorentzian. Both routes share one globally
adaptive Gauss-Legendre quadrature: every piece carries a 16/32-node error
estimate, and the piece with the largest error relative to the requested
tolerance of the whole integral is bisected until the summed estimate meets
`krook_atol + krook_rtol*|integral|` in every component. Pieces whose
contribution is negligible are therefore not refined, and refinement grades
geometrically into endpoint layers. The passing and trapped pitch ranges end
`1e-8*etatp` from the trapped-passing boundary, where the bounce weight varies
on that scale; resolving it needs about 27 bisections, so `krook_max_depth`
(maximum bisections of one piece) defaults to 40. A piece budget of `2**14`
bounds the work for an unresolvable integrand.

The quadrature can only converge if the weight is smooth to its tolerance.
The bounce integral of `H_n` uses an absolute tolerance relative to
`taub*max|H_n|`, so the weight is quadratic in the perturbation amplitude
(`test_transport_weight_amplitude`), and orbit integrations use relative
tolerance `1e-13`: at `1e-9` the weight was noisy by about `1e-3`, and the
frequency splines non-monotone, within `1e-5*etatp` of the trapped-passing
boundary. Integrating parallel velocity independently leaves energy errors
that are amplified by `1/|1-eta*B_max|` near the separatrix, making the weight
noisy even at the tighter ODE tolerance (about `1e-4` at `eta/etatp-1 = 1e-8`
on the trapped side); refining the pitch quadrature cannot remove that noise.
Passing transport therefore evaluates parallel velocity from the energy
invariant `v_parallel**2 = v**2*(1-eta*B)` at each orbit step. Trapped transport
integrates the closed orbit in the angle `chi`, with
`theta = theta_m - Delta*cos(chi)` between the turning points, so that
`dt/dchi` follows from the energy invariant and is regular at the turning
points. Close to them, `1-eta*B` is integrated from `dB/dtheta` instead of
being formed as a difference. A first pass gives the orbit period, which sets
the bounce frequency of the phase `exp(-i*mth*Omth*t)` so that the harmonic
closes on the orbit; the interpolated `Om_th` is off by up to `1e-5` in its
logarithmic extrapolation below `eta/etatp-1 = 1e-6`.
Passing transport likewise integrates exactly one transit in the poloidal
angle travelled, `theta = th0 + direction*phi` with
`dt/dphi = 1/|v_parallel*h^theta|`, and closes the phase
`exp(-i*(mth+q*mph)*Omth*t)` with the period of a first pass. The spline
`Om_th` is off by about `5e-7` in the passing interior for ITER C12 at
`s = 0.09`; a time integration over `2*pi/|Om_th|` with that frequency then
changed the squared weight of the strongly cancelling `mth = -11` harmonic
by a factor of about 200 and left it too noisy for the Krook pitch quadrature.
The two passes take different ODE steps, so their periods differ by about
`1e-11`; closing the second-pass phase with the first-pass period left a
mismatch of `(mth+q*mph)*2*pi` times that. Harmonics that cancel to
`|<H_n>| ~ 1e-7*max|H_n|`, as at the low-rate ITER C12 resonances of slow
passing particles, then had a weight noisy by `1e-4` in eta, and the Krook
quadrature exhausted its piece budget with an error estimate stuck at
`3.5e-6` of the integral. The second pass therefore also integrates
`d/dOmth` of the harmonic and shifts it to the frequency of its own period
to first order (`test_krook_passing_harmonic` checks harmonics cancelling to
`6e-8` against a spectral transit reference to `5e-12*max|H_n|`). The
remaining ODE error, about `1e-12*max|H_n|`, still limits the relative
accuracy of the most strongly cancelling harmonics.

`krook_rtol` defaults to `1e-5`, `krook_atol` to zero and
`krook_max_depth` to 40. The uniform probes cannot see structure of `R` on
the logarithmic scale of `|eta/etatp-1|`: the frequency splines have nodes
graded geometrically towards the trapped-passing boundary and switch to a
logarithmic extrapolation at `|eta/etatp-1| = 1e-6`, and the near cancellation
of `mth*Om_th + mph*Om_ph` there turns their interpolation error into
wiggles with extra extrema (TC24 GPEC: `R` relative amplitude about `2e-7`
at `1-eta/etatp = 1.4e-6`, D at `s = 0.0676`; several in
`eta/etatp-1 = 1e-6..4e-5`, T at `s = 0.79` and `0.86`). When the inverse
residual map meets a derivative of the wrong sign for its branch, the
extrema around that pitch are bracketed by samples graded geometrically
towards it, become knots, and the integration restarts (at most 32 times).
Nonconvergence, a residual left non-monotonic after these repairs, a
non-finite residual or an unresolved stationary weight terminates the
calculation. Finite sampling cannot certify an arbitrary oscillatory callback.
Collisionality scans must check convergence in pitch probes, quadrature
tolerance, velocity steps and harmonic cutoff independently.

`test_krook` compares finite-interval affine detuning with constant/linear
amplitudes against analytic angle and logarithmic moments for both slope
signs, interior/boundary/exterior roots and widths from `1e-12` to `1e3`.
Quadratic and hidden cubic folds use high-precision analytic oracles checked
against independently split quadrature. A `1e-8` endpoint layer is checked
against its closed form, and a jittered tail below `1e-21` of the peak must
converge without refinement to the depth limit. `test_krook_passing_weight`
compares near-separatrix passing Hamiltonian averages in both directions with
an independent angular quadrature, including convergence of that reference.
`test_krook_passing_harmonic` compares cancelling passing harmonics in both
directions with a closed-orbit angular reference for supplied transit times
biased by up to `1e-6`, within `1e-10` of the integrand scale `max|H|`.
`test_krook_trapped_weight` compares trapped bounce harmonics at
`eta/etatp-1 = 1e-8, 1e-7, 1e-6` with a panel quadrature in the turning-point
angle with cumulative bounce time.
`test_krook_transport` checks actual field transport against a direct pitch
quadrature and verifies exact equality of disabled and zero-width dispatches. Configuration and selected-profile
tests reject invalid and nonfinite inputs, including in Fast builds.

Regenerate the leaves with a FortSym installation at the recorded revision:

```sh
cmake -S gen -B build/codegen-cmake -DCMAKE_PREFIX_PATH=/path/to/fortsym
cmake --build build/codegen-cmake --parallel 4
build/codegen-cmake/gen_krook src/generated/krook_leaf.inc FORTSYM_REVISION
```

## Species profiles

`profile_species=1` (default) selects the first ion in `plasma.in`; `2` selects
the second ion and `0` selects electrons. The selected profiles determine
mass, signed charge, thermal speed, thermodynamic forces and torque density.
Electron density follows quasineutrality, `ne=Z1*n1+Z2*n2`; electrons use
`qi=-qe` and `me=9.109382e-28 g`. Linear electrons skip the legacy ion
collision kernel; nonlinear electrons are rejected. Selected density,
temperature, mass and charge must be physically admissible and finite, and
selected profile derivatives must be finite.

Rotation input remains the selected species' Mach number. To hold signed
physical `OmegaE` fixed across species, use `M_t=R0*OmegaE/vth` and undo the
`efac/bfac` scaling when it is present.

For a general mixture, [set `species_file`](species.md) to supply D, T and
electrons together with a common physical rotation and explicit profiles.
NEO-RT then derives each species' Mach number and returns all species results
in one native invocation.


For transport without radial orbit width, a Fourier harmonic can be
smaller than the error of the orbit integral. Its squared weight then cannot
meet a purely relative pitch tolerance. This happens for passing orbits and for
high trapped harmonics, whose `|H|/S` falls to the orbit-integration floor in
the deep well. NEO-RT integrates a nonoscillating
envelope `S = (2-eta Bmin) |pertfile_scale| sum|dB_m|/Bmin` with the same
Lorentzian and transport prefactors. For an analytic perturbation, its amplitude
is `|epsmn|`. The amplitude error estimate is `e S`, where
`e = [1+2 pi |mth+nq|] bounce_rtol` for passing orbits and
`e = [1+2 pi (|mth|+|nq|)] bounce_rtol` for trapped orbits, whose phase holds the
secular `mth Omega_theta t` and the bounded `nq theta`. The squared weight then
errs by at most `2 e S |H| + e^2 S^2`. Both terms are integrated with the same
Lorentzian and prefactors: `S^2` to the relative tolerance, the cross term
`S |H|` to 1 % with its own floor `e S^2` as absolute tolerance. Their sum is the
component-wise absolute budget of the pitch quadrature, which therefore scales
quadratically with the perturbation. The earlier Young bound
`(1+1/rtol) e^2 S^2` exceeded the value itself where the harmonic is resolved
near the trapped-passing boundary. These are numerical error estimates, not
rigorous ODE enclosures. The requested relative tolerance and quadrature work
limits are unchanged; unresolved quadrature errors beyond this budget still abort.
`test_krook_c12_cancellation` retains the exact opposite-rotation, constant-Krook,
multiplier-0.01 deck and its failing T counter-passing mth=+6 velocity node.
It checks termination, the bounded tiny-channel reference and amplitude scaling.
`test_krook_trapped_harmonic` retains the circular large-aspect-ratio field
and its failing trapped mth=-40 node, and compares it with fixed panels in
`log(eta-etatp)` within the 2e-3 orbit-error budget of that node;
the separate MARS reference retains its strict 1e-6 bound.
