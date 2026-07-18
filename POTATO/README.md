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
`&potato_nml` namelist), the g-eqdsk named in `field_divB0.inp`, `convexwall.dat`,
and the profile files.

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
