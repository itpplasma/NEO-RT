# POTATO

POTATO computes guiding-center orbit resonances and the resonant NTV torque for
an axisymmetric tokamak equilibrium (EFIT g-eqdsk). The resonance run
(`itest_type=3`) builds a canonical-frequency grid per orbit class and solves
`m*Omega_b + n*Omega_phi = 0` for every poloidal harmonic.

## Build

POTATO uses CMake with the Ninja generator (the module dependency order across
the fetched libneo and fortnum needs Ninja's dyndep):

```bash
cmake -S . -B build -G Ninja
cmake --build build -j
```

The default architecture flags target the build host. To build a portable
binary, override the CMake list `POTATO_ARCH_FLAGS`, for example
`'-DPOTATO_ARCH_FLAGS=-march=x86-64-v3;-mtune=generic'` on x86-64, or
`-DPOTATO_ARCH_FLAGS=` to omit POTATO's architecture flags. Dependencies retain
their own build settings and must also support the destination CPU.

`build/potato.x` reads its inputs from the working directory: `potato.in` (the
`&potato_nml` namelist), the g-eqdsk and convex-wall file named in
`field_divB0.inp`, and the profile files.

Resonant-torque runs are Maxwellian-integrated by default. Set
`monoenergetic_x=E/T_ref` to evaluate one fixed total-energy slice instead,
where `T_ref` and the reference potential are taken on the magnetic axis. The
energy measure is normalized as `delta(x-x0)`, so no finite energy-bin width is
applied. For finite-width orbits through varying profiles this remains an
energy-resolved nonlocal torque, not a local `D11`; at zero potential and in the
local-profile limit it supplies the corresponding fixed-energy reference. See
the finite-orbit formulation and AUG application in
[Buchholz et al. (2022)](https://doi.org/10.1088/1742-6596/2397/1/012012).

Thermal runs integrate all `nenerg` cells with the midpoint rule:
`H_i = Phi_min + (i - 1/2)*DeltaH`, `DeltaH = (thermen_max*T_ref + Phi_max - Phi_min)/nenerg`,
for `i = 1..nenerg`. Cell 1 is an interior node, not a zero-energy endpoint;
its Maxwellian and thermodynamic-force weight can contribute with either sign.
Dropping it changes the physical integral. At zero potential, a one-cell run
with `thermen_max=3` samples `H/T_ref=1.5` and must equal three times the
delta-normalized fixed-energy torque at `monoenergetic_x=1.5`.

The g-eqdsk and convex-wall file are a matched equilibrium pair. The convex
wall is the outer boundary of the complete computational domain, not the
separatrix, and must enclose the LCFS with enough margin for every orbit being
traced. POTATO uses this one configured wall for both the field-coordinate
mapping and orbit-loss detection; no second hard-coded `convexwall.dat` is
needed. By default, `edge_extension = .true.` permits orbit integration across
the LCFS into the scrape-off layer. Set it explicitly to `.false.` for a
closed-flux-only run. Edge extension does not turn the convex wall into the
LCFS or remove the outer-domain requirement.

For a MARS field already converted to Boozer harmonics, generate POTATO's
single-`n` cylindrical perturbation without the legacy converter's radial
Gaussian filter:

```bash
python tools/boozer_npz_to_bmod_n.py chartmap.nc components.npz bmod_n.dat \
    --component total --n-tor=+3 --s-max=0.95
```

Here `chartmap.nc` supplies the accepted Boozer-surface geometry and
`components.npz` is the provenance product from `rmp_torque mars_to_boozer`.
The signed `n` must also be used as `n_tor` in `potato.in`. The converter writes
a JSON sidecar, performs no smoothing or fit, uses `s_tor` explicitly, writes
zero outside the outer mapped surface, and adds a zero-valued rectangular
margin so the POTATO spline is not normally evaluated at a clamped nonzero
boundary. Production target orbits must remain inside that mapped surface.
The displayed `n=+3` is the current TC24 MARS-to-CCW result:
`phi_CCW=-phi_MARS` and `n_CCW=-RNTOR=+3`. It is not a universal MARS
default; another field must supply its own signed laboratory-coordinate map.

The perturbation converter accepts both NetCDF3 and HDF5 chartmaps. Geometry
and magnetic harmonics use the `rho` nodes, with `s_tor=rho**2`; the separate
`s` array belongs to the covariant vector potential. Legacy HDF5 maps without
`rho` retain their explicitly supplied `s` geometry axis.

The matching profile converter also keeps the coordinate and electric-field
conventions explicit:

```bash
python tools/neo_rt_profiles_to_potato.py profile.in plasma.in components.npz \
    profile_poly.in --r0-cm=641.0903942 --psi-span-tm2=11.88279543 \
    --relation-sign=1
```

For profiles whose edge cannot be represented by a global polynomial, add
`--format=table` and point `profile_file` at that output. The native table starts
with `POTATO_PROFILE_TABLE_V1`, then a row count and four columns:
`s_pol, n[cm^-3], Ti[eV], Phi[statvolt]`. POTATO interpolates these physical
samples with the same natural cubic spline library used by NEO-RT and computes
the derivatives from that interpolation. Table coordinates must increase and
cover `s_pol=0..1`; density and temperature must be positive. A query beyond
table support stops with an explicit error by default. The existing polynomial format
remains the default. Validate values and gradients over the orbit domain;
a small global polynomial residual can hide a large relative edge error.

`potato_profile_readback.x FILE Ealpha Z escale fscale psi_axis psi_sep QUERIES`
reads either profile format and writes `s_pol,n,T_hat,Phi_hat,dn/dpsi,
dT_hat/dpsi,dPhi_hat/dpsi`. The query file contains a count followed by one
`s_pol` per row. This provides a native value/derivative check with the actual
energy, charge, scale factors and signed poloidal-flux span.

It selects the physical ion using charge and nonzero density, maps
`s_tor -> s_pol=rho_pol^2`, and records the polynomial residuals. The JSON
sidecar also records the signed relation
`dPhi/dpsi_pol = relation_sign*Omega_E/c`; opposite signs are separate
convention diagnostics, never an unrecorded curve flip. For the displayed
right-handed direct-EQDSK chart, the native field equations give
`Omega_phi,E=c*dPhi/dpsi_pol`, so `relation_sign=+1` maps an already-CCW
`Omega_E` profile. The current TC24 NEO profile has first been serialized as
`Omega_E,CCW=-Omega_E,MARS`; the converter does not apply that MARS-to-CCW
change a second time. The first value in the NEO `plasma.in` header is a
validated radial-row count, not a required value of 50. The compiled NEO-RT
profile schema needs only `s_tor,M_t`; the converter derives
`v_th=sqrt(2*T_i/m_i)` from the charge-selected `plasma.in` species. If a
legacy third `v_th` column is present, it is checked against that source rather
than used as an independent temperature profile.

## Optional absorbing profile boundary

`profile_exit_absorb = .true.` in `potato_nml` assigns zero resonant weight
when either physical bounce trace leaves the tabulated profile support. The
profile evaluator returns constant finite edge values while recording the
excursion; returning inside support does not restore the orbit. Each orbit
resets this thread-private state. Box counting independently rejects a trace
that leaves support and deposits zero torque for it. Polynomial profiles are
unaffected because they have no tabulated support boundary.

The option defaults to false, preserving the existing outside-support error.
The canonical-axis closure is checked only after physical orbit admission;
the pre-existing canonical SOL exclusion still skips the orbit with zero
weight. Undefined canonical queries also receive zero weight in absorbing mode.
For each energy the log reports `absorbed`, `total` and `box_trace_absorbed`.
The first two count roots on the retained J_perp sampling grid (including
zero-weight SOL roots in `total`); they do not measure excluded absolute weight.
`box_trace_absorbed` counts additional exits during the independent box trace.

## Optional canonical-axis distribution closure

`canonical_axis_constant_extension = .true.` selects a model in which the
canonical profiles n(psi*), T(psi*) and Phi(psi*) are continuous constants
below the magnetic-axis label. At fixed total energy H this gives
`d ln F0(H,psi*)/d psi* = 0`, so those canonical resonances have zero drive.
The signed normalized label `(psi*-psi_axis)/(psi_edge-psi_axis)` determines
this region for either flux orientation. Both physical bounce traces must
pass before the closure is applied. Physical-orbit profile support checks,
in-domain weights and the separate SOL exclusion remain unchanged.

The option defaults to false: undefined canonical profile queries still fail.
This is an explicit distribution-model boundary condition, not a collision
operator, finite-orbit convergence result or normalization adjustment.

## Running with OpenMP

The grid build and the per-mode root search run in parallel with OpenMP. Use one
thread per physical core:

```bash
./run_potato.sh            # sets the thread count and binding, then runs build/potato.x
```

`run_potato.sh` sets `OMP_NUM_THREADS` to the physical-core count and pins one
thread per core. The manual equivalent:

```bash
OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close ./build/potato.x
```

Do not add the SMT (hyperthread) siblings. The solver is floating-point bound, so
a second thread per core fights for the FPU instead of adding throughput. On a
16-core/32-thread host the rho_pol=0.99 #30835 sweep scales as 620 s (1 thread),
290 s (4), 191 s (8), 148 s (16): 4.2x at 16. All 32 logical threads take ~187 s
on an idle machine and 400-600 s under any concurrent load, because they
oversubscribe the 16 cores.

The reslines output is identical for every thread count.

## Perturbation in Boozer coordinates

The default perturbation reader uses `bmod_n.dat`. To evaluate the original BC
spectrum with [NEO-RT](https://github.com/itpplasma/NEO-RT)'s reader instead, add
these entries to `potato_nml`:

```fortran
use_boozer_perturbation = .true.
boozer_geometry_file = 'boozer_geometry.nc'
boozer_perturbation_file = 'in_file_pert'
```

Create the geometry chart once from the corresponding axisymmetric equilibrium
BC file (requires NumPy and SciPy):

```bash
python3 tools/boozer_geometry.py equilibrium.bc boozer_geometry.nc
```

The chart samples only smooth equilibrium geometry. libneo's inverse returns
normalized toroidal flux `s`, `theta_B`, and `phi_B` at geometric phi=0; its radial
metadata must be `rho_convention=psi_tor_norm`. POTATO evaluates the unchanged
NEO-RT BC spline/Fourier sum at `(s,theta_B)` and multiplies by
`exp(i*n_tor*(phi_B-phi_geom))`. The result is in Gauss, with no conjugation or
factor two. The signed BC toroidal mode must equal `n_tor`. NEO-RT's radial
clamping at the perturbation file endpoints is retained. Coordinates outside the
geometry chart fail explicitly.

The inverse uses the equilibrium BC flux label rather than deriving it from
POTATO's EQDSK psi. A single psi-to-s table requires the two equilibria to have
identical surfaces; a mismatch can displace narrow perturbation layers. Geometry
interpolation accuracy must also be checked against the original BC at those
layers. The chart grid controls geometry accuracy, without resampling the
perturbation. Shared splines load once; each OpenMP worker owns its reader cache
and previous inverse coordinates.

`potato_perturbation_readback.x query.bin result.bin` compares the selected POTATO
reader with NEO-RT at supplied flux coordinates. Input is a native-endian int32
point count followed by six float64 columns (Fortran order): R, Z [cm], s,
theta_B, real and imaginary parts of the toroidal phase. Output contains two
complex128 arrays: POTATO amplitudes and NEO-RT reference amplitudes, in Gauss.

## SIMPLE-compatible invariant handoff

Single-orbit (`itest_type=4`) and radial-frequency-scan (`itest_type=5`) runs
write `potato_invariants.dat` without changing the existing `fort.100` or
`freq_scan.dat` schemas.  The versioned file contains the local state and the
three axisymmetric invariants needed by SIMPLE:

```text
H0       = p^2 + phi_elec
J_perp   = p^2 (1 - xi^2)/B
psi_star = psi + rho0 p xi h_phi = (c/q) P_phi
```

Here `p=v/v0`, `xi=v_parallel/v`, `v0=sqrt(2 E_ref/m)` in cm/s, `B` is in
POTATO's native gauss convention, `phi_elec=q Phi/E_ref`, and `psi_star`,
`psi_axis`, and `psi_edge` use the same native g-eqdsk flux gauge and units
(G cm^2).  Time remains POTATO's `tau=v0*t`; SIMPLE's private symplectic
sqrt(2) scaling is not part of this file.

SIMPLE currently supports only the magnetic Hamiltonian.  Therefore status is
zero only when `phi_elec` is zero to numerical tolerance; a nonzero-potential
row still records POTATO's correct total `H0`, but has status 2 and must not be
converted.  Status 1 means invalid normalization input and status 3 means the
field evaluation failed.

For a scan, rows retain the requested increasing-`R` order (HFS to LFS when
the input bounds are ordered that way), matching POTATO's `R_c` convention.
No `rho_pol`/`rho_tor` identification or outer/inner orbit selection is
implied by the handoff.

## Fixed-energy invariant contour scan

`potato_resonance_contour.x` scans a rectangular `(rho_pol, J_perp)` grid at
fixed kinetic energy. For a 5 keV deuteron scan, set `E_alpha = 5d3` and
`contour_enkin = 1d0`; `contour_enkin` is normalized to `E_alpha`.

```text
contour_rho_min = 0.05d0
contour_rho_max = 0.95d0
contour_nrho = 21
contour_jperp_min = 1d-6
contour_jperp_max = 8d-5
contour_njperp = 21
contour_enkin = 1d0
contour_sigma = 1
```

Run it from a directory containing the normal POTATO equilibrium, wall, profile,
and `potato.in` inputs:

```bash
fo exec --cwd /path/to/run potato_resonance_contour.x
```

The resulting `potato_resonance_contour.dat` retains every requested grid point
and contains `H0`, `J_perp`, `psi_star`, the guiding-centre start, `taub`,
`delphi`, physical angular frequencies, and `ierr`. Keep failed rows for the
rectangular grid and mask them using `ierr`; do not interpolate across them.
