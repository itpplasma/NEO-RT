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
