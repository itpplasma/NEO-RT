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

`build/potato.x` reads its inputs from the working directory: `potato.in` (the
`&potato_nml` namelist), the g-eqdsk and convex-wall file named in
`field_divB0.inp`, and the profile files.

The g-eqdsk and convex-wall file are a matched equilibrium pair. The convex
wall is the outer boundary of the complete computational domain, not the
separatrix, and must enclose the LCFS with enough margin for every orbit being
traced. POTATO uses this one configured wall for both the field-coordinate
mapping and orbit-loss detection; no second hard-coded `convexwall.dat` is
needed. Setting `edge_extension = .true.` permits orbit integration across the
LCFS into the scrape-off layer; it does not turn the convex wall into the LCFS
or remove the outer-domain requirement.

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

The matching profile converter also keeps the coordinate and electric-field
conventions explicit:

```bash
python tools/neo_rt_profiles_to_potato.py profile.in plasma.in components.npz \
    profile_poly.in --r0-cm=641.0903942 --psi-span-tm2=11.88279543 \
    --relation-sign=1
```

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
