# Fortran baseline (measured)

Host: linux x86-64, gfortran (f95), -O3 march=x86-64-v2, Release.

| metric | value | method |
|---|---|---|
| clean rebuild (project + spline + vode + fortplot, 287 files) | 16.1 s wall (75.4 s user, -j) | cmake --build build --clean-first |
| ctest (3 Fortran unit/integration tests) | 0.88 s, 3/3 pass | ctest |
| run neo_rt.x driftorbit (s=0.5 base case) | 0.16 s wall | /usr/bin/time -p |

Golden gate: test/golden_record compares 5 output files at rtol=1e-8, atol=1e-15.
ODE: DVODE Adams (method_flag=10), rtol=1e-9, atol=1e-10, root-finding events (nevents=2).
Golden-path external deps actually used: spline (111 lines) + DVODE. No BLAS/LAPACK/NetCDF
beyond LAPACK dptsv called inside spline_coeff.

## Determinism (critical for the gate)

The golden gate is only meaningful single-threaded. Measured:
- OMP_NUM_THREADS=1: Fortran reproduces its own golden.h5 bit-exact
  (max abs diff 0.0 across all 5 arrays, all 9 cases s=0.1..0.9). 9/9 pass.
- Multithreaded (default OpenMP): fresh-run vs golden differs at ~1e-4 relative
  (output/torque/integral), far above rtol=1e-8 -- run-configuration dependent.
  This is why golden_record is NOT part of ctest; only test_timestep,
  test_parallel, test_neort_lib are.

Therefore every cross-language port is gated with OMP_NUM_THREADS=1 against the
single-threaded golden.h5. Gate runtime: 58.4 s wall for 9 cases single-threaded
(golden cases use inp_swi=9, pertfile=T, vsteps=512 -- heavier than the
examples/base inp_swi=8 case which runs in 0.16 s).

golden.h5 is gitignored (Makefile .gitignore) and persists in the working tree
across branch checkouts, so the same reference gates port/c, port/rust,
port/julia. Run: port_c/bench/run_gate.sh /path/to/neo_rt_executable
