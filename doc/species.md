# Shared species mixtures

Set `species_file='species.in'` in `&params` to solve every member of a neutral
mixture in one `neo_rt.x <runname>` invocation. The magnetic field, perturbation,
surface and physical electric rotation are shared. Each species has its own
signed charge, mass, profiles, thermal speed, Mach number and transport result.
An empty `species_file` retains the existing isolated-species input paths.

The file contains two namelists. All values and derivatives below are explicit;
electron density and temperature are not inferred from ion inputs.

```fortran
&species_count
 n_species=3
/
&species
 reference_s=0.09
 omega_e=-700.0
 domega_e_ds=0.0
 coulomb_log=20.0
 species_name='D','T','e'
 mass_kg=3.3435837724e-27,5.0073567446e-27,9.1093837139e-31
 charge_e=1.0,1.0,-1.0
 density_m3=5.0e15,5.0e15,1.0e16
 temperature_ev=2000.0,2000.0,2000.0
 density_gradient_m3=0.0,0.0,0.0
 temperature_gradient_ev=0.0,0.0,0.0
/
```

Mass is in kg, density in m^-3 and temperature in eV. `charge_e` is signed charge
in elementary-charge units. `omega_e` is the signed angular electric rotation
frequency in s^-1, with no factor of 2π. Every derivative is with respect to
dimensionless normalized toroidal flux `s`. `coulomb_log` is an explicit positive
common Coulomb logarithm for collision models that consume pairwise rates.

The file defines an **affine local expansion** about `reference_s`:
`n(s)=n(reference_s)+density_gradient_m3*(s-reference_s)`, with analogous
temperature and rotation expansions. For a surface scan, evaluate the desired
profiles and derivatives at each surface and set `reference_s` to that surface.
This input does not fit or extrapolate a global profile from a single point.
Both density and its gradient must obey charge neutrality. Names must be unique;
all physical inputs must be finite; masses, densities and temperatures must be
positive, and charges nonzero. An affine expansion that produces invalid values
at the requested surface is rejected.

NEO-RT converts SI data to its native cgs units once at species selection. The
common physical rotation and derivative are multiplied by the existing
`efac/bfac` scaling. Species Mach numbers then follow `M_t=R0*OmegaE/vth`,
including the thermal-speed derivative in `dM_t/ds`. The general mixture path
requires `nonlin=.false.`; the older nonlinear attenuation kernel accepts only
its historical ion backgrounds.

The library exposes `species_state_t`, `mixture_t`, `neort_prepare_mixture` and
`neort_compute_mixture_at_s`. Initialize the field with `neort_init`, prepare an
explicit in-memory mixture, then compute at the selected surface. The evaluated
background list is available to collision models. Independent species solves
use the existing surface-parallel state architecture and OpenMP; output is
written serially after all species finish. These data and run interfaces do not
by themselves change a diagonal Krook model into a conserving, coupled
Fokker–Planck operator.

Set `krook_collisions=.true.` and select the physical width source:

| `krook_rate_source` | Width before the population geometry factor |
|---|---|
| 0 (default) | Supplied `krook_nu` in s^-1, unchanged from PR201 |
| 1 | `krook_rate_scale*sum_b nu_D_ab(vth_a)`, then `krook_model` 0/1/2 |
| 2 | `krook_rate_scale*sum_b nu_D_ab(v)` at the actual speed; requires model 0 |

Sources 1/2 require an explicit mixture and use every background, including
electron–electron, ion–ion and electron–ion pairs, with the supplied common
Coulomb logarithm. `krook_rate_scale` is finite and nonnegative; zero restores
the exact collisionless dispatch. Source 1 can use constant, inverse-cube or
harmonic Krook shapes. Source 2 uses the Maxwellian deflection-rate shape
directly, with no additional inverse-cube or harmonic factor. Population
geometry and passing harmonic conventions retain their explicit PR201 controls.
The thermal deflection rate and distinct self-reference amplitude
`Gamma_aa/vth_a^3` are exported per species before scaling and geometry factors.

These widths relax the nonadiabatic response
`g=delta_f-delta_H*dF0/dH`, as in the compared MARS/PENTRC Krook numerators.
They do not damp plain `delta_f`. This choice matters for finite-width
electromagnetic torque away from resonance. Deflection-shaped scalar damping
still has no pitch diffusion, energy diffusion or field-particle coupling.
Lorentz, test-particle FP and conserving multispecies FP require a new coupled
complex response solver; they are not implemented by these selectors.

With `output_format='hdf5'` or `'both'`, one `<runname>.h5` uses output layout
version 2. `/species/<name>` contains explicit SI metadata and the native
`config`, `magfie`, `transport`, `torque` and `harmonics` groups. `/mixture`
records the common reference state, profile model and summed physical torque.
`torque_by_population` contains co-passing, counter-passing and trapped species
sums; `torque_total` sums all populations. Torque is electromagnetic perturbation work for the stated nonadiabatic-g
Krook model, evaluated with the legacy flux-force-equivalent formula
`dTphi_int/ds` in erg. Broadened diffusion coefficients are not independently
qualified physical particle fluxes outside the collisionless limit. Diffusion coefficients have
species-dependent normalizations and are therefore reported per species.

Text output writes the existing files with `<runname>_<species>` prefixes,
plus `<runname>_species.out` with one metadata/torque row per species and
`<runname>_mixture_torque.out` with the species-summed torque. Isolated-species
HDF5 output retains layout version 1.
