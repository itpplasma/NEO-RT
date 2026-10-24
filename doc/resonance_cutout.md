---
title: Resonance cutout
---

# Resonance cutout

`resonance_cutout` is an optional regularization of the collisionless
resonance weight. The default `0` leaves the collisionless quadrature
unchanged, bit for bit. The option acts only on the collisionless path
([src/transport.f90](../src/transport.f90), `compute_collisionless_sample`) and
is rejected together with finite Krook broadening, which already replaces the
delta function by a Lorentzian ([Krook collisions](krook.md)), and together
with the [collisional boundary-layer closure](collisional_layer.md), whose
trapped-boundary layer overlaps the cutout layer next to `eta_tp`.

## Where the collisionless weight diverges

For every speed `v` the collisionless quadrature replaces
`delta(R(eta))`, with `R = mth*Omega_theta + mph*Omega_phi`, by `1/|R_eta|` at
each root of `R` in the pitch interval of one orbit class. At a fold of `R`,
where two roots merge and `R_eta = 0`, this weight diverges. Near the fold the
roots exist on one side of a critical speed `v_f`, their distance from the
fold scales as `sqrt(|v - v_f|)` and the weight as `1/sqrt(|v - v_f|)`. The
speed integral therefore remains finite, but a velocity node that lands close
to `v_f` gives a spike whose size depends on the node position. The same
breakdown occurs next to the trapped-passing boundary, where `Omega_theta`
vanishes logarithmically: there `|R_eta/R_eta_eta|` is of the order of the
distance `|eta - eta_tp|`, and resonances with `mth != 0` at slow rotation sit
exactly there because they need a small `Omega_theta`. A flat
resonance function (`R_eta` and `R_eta_eta` both small, as for transit
resonances close to a rational surface at small `Omega_tE`) is not covered by
this criterion and is left unchanged.

## Layer width from pitch-angle scattering

With pitch-angle scattering, the resonant response solves
`-D_eta g'' + i R(eta) g = S` locally. Balancing the two terms gives the
characteristic width scalings

* for a simple root (`R` linear), the classical collisional boundary-layer
  width `delta_lin = (D_eta/|R_eta|)**(1/3)`, which is the width used by the
  [collisional boundary-layer closure](collisional_layer.md);
* at a fold (`R = R_f + R_eta_eta*(eta - eta_f)**2/2`), the width
  `delta_fold = (D_eta/|R_eta_eta|)**(1/4)`.

The delta-function limit requires that `R` stays linear across `delta_lin`,
that is `delta_lin < d` with `d = |R_eta/R_eta_eta|` the distance of the root
from the fold. Since `R_eta = R_eta_eta*d` there, `delta_lin < d` is the same
inequality as `d > delta_fold`, so a root lies inside the resonance layer
when `d < delta_fold`. This linear-layer-fit criterion fixes the order-unity
normalization used here. The width is neither fitted nor obtained by solving
the collisional boundary-layer equation, whose quadratic Taylor coefficient
contains an additional factor of one half. Next to `eta_tp`, `d` is the local linearity
scale of the logarithmic frequency, and the criterion marks roots inside the
collisional layer of the boundary. `R_eta_eta` is obtained by central
differences of the analytic `R_eta` with a step of `1e-4` of the orbit-class
pitch interval, but at most a quarter of the distance to the nearer interval
end, so that the stencil resolves the logarithmic curvature and never crosses
`eta_tp`.

`D_eta = 2*nu_d*eta*<(1 - eta*B)/B>` is the orbit-averaged pitch diffusion
coefficient in the large-aspect-ratio field `B = B0*(1 - eps*cos(theta))`,
evaluated with complete elliptic integrals for trapped and passing orbits
(`bounce_pitch_diffusion` in [src/resonance_cutout.f90](../src/resonance_cutout.f90)).
Its trapped branch coincides with `D_eta_ba` of the collisional boundary-layer
closure. `nu_d` is the energy-dependent deflection frequency of the selected
species against the full background (`deflection_frequency` in
[src/collisional_layer.f90](../src/collisional_layer.f90)). `test_resonance_cutout`
checks both branches against direct theta quadrature to within `2*eps`, which
is the size of the shared large-aspect-ratio modulus approximation.

## Variants

| `resonance_cutout` | Root inside the layer (`d < w*delta_fold`) |
| --- | --- |
| `0` | Off, default: weight `1/|R_eta|` everywhere. |
| `1` | Zero: the root contributes nothing. |
| `2` | Boundary-value fill: the weight is frozen at its value on the layer boundary, `1/(|R_eta_eta|*w*delta_fold)`. It is continuous there and never exceeds `1/|R_eta|`. |

`resonance_cutout_width = w >= 1` (default `1`) widens the layer. It applies to
both variants, and roots outside the widened layer are unchanged.

## Choosing a wider cutoff

The layer above contains only collisional broadening. Finite orbit width adds
a broadening of the resonance frequency, `delta_R_fow ~ |dR/dr|*Delta_r`,
with the radial orbit width `Delta_r ~ q*rho/sqrt(eps)` for trapped and
`~ q*rho` for passing orbits. On a fold this frequency spread corresponds to
the pitch width `delta_fow = sqrt(2*delta_R_fow/|R_eta_eta|)`, and a cutoff
covering both effects is `w = max(1, delta_fow/delta_fold)`. NEO-RT does not
evaluate `dR/dr`, so `w` is supplied by the user as a bound. The ratio of the
collisional fold width to the linear width at the layer edge is one, so
`w = 2` and `w = 4` correspond to frequency widths larger by `4` and `16`.

## Related regularizations elsewhere

* NEO-RT Krook (`krook_collisions`): Lorentzian of width `nu` in `R` in place
  of the delta function ([doc/krook.md](krook.md)).
* PENTRC (GPEC): resonance denominator `i*(l*omega_b*sqrt(x) + n*(omega_E +
  omega_D*x)) - nu` with `nutype` `zero`, `small` (`1e-5*omega_E`), `krook` or
  `harmonic`, and an optional complex energy contour `ximag`
  (`pentrc/energy.f90`, `xintgrnd`; Logan et al., Phys. Plasmas 20, 122507, 2013).
* GPEC: the ideal displacement is multiplied by
  `singfac**2/(singfac**2 + reg_spot**2)`, `singfac = m - n*q`, default
  `reg_spot = 0.05` (`gpec/gpeq.f`, `gpec/gpec.f`). This is a smooth cutout of
  the field input at rational surfaces rather than of the kinetic resonance.
* MARS-K: effective collision frequency in the resonance denominator and
  analytic subtraction of the pitch-angle singularities (`KLAM0`, `KSF0`,
  `KSGH0` in `kinetic.f`).

## Generated kernels

The pointwise expressions are derived in `gen/gen_resonance_cutout.f90` and
committed as `src/generated/resonance_cutout_leaf.inc`. Ordinary builds need no
symbolic algebra package. With FortSym available, build the generator through
`gen/CMakeLists.txt` and run `gen_resonance_cutout
src/generated/resonance_cutout_leaf.inc FORTSYM_REVISION`; regenerating with
the recorded FortSym revision must reproduce the include file.
