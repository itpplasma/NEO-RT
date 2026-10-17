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
rule integrates the same finite Lorentzian. Both routes use adaptive
Gauss-Legendre quadrature and sum their error estimates before checking the
requested total tolerance.

`krook_rtol` defaults to `1e-5`, `krook_atol` to zero and
`krook_max_depth` to 20. Nonconvergence, a detected wrong-sign inverse
derivative or an unresolved stationary weight terminates the calculation.
The monotonicity probes cover the native frequency-spline use case and simple
folds; finite sampling cannot certify an arbitrary oscillatory callback.
Collisionality scans must check convergence in pitch probes, quadrature
tolerance, velocity steps and harmonic cutoff independently.

`test_krook` compares finite-interval affine detuning with constant/linear
amplitudes against analytic angle and logarithmic moments for both slope
signs, interior/boundary/exterior roots and widths from `1e-12` to `1e3`.
Quadratic and hidden cubic folds use high-precision analytic oracles checked
against independently split quadrature. `test_krook_transport` checks actual
field transport against a direct pitch quadrature and verifies exact equality
of disabled and zero-width dispatches. Configuration and selected-profile
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
