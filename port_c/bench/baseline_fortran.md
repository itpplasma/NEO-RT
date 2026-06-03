# Fortran baseline (measured)

Host: linux x86-64, gfortran (f95), -O3 march=x86-64-v2, Release.

| metric | value | method |
|---|---|---|
| clean rebuild (project + spline + vode + fortplot, 287 files) | 16.1 s wall (75.4 s user, -j) | cmake --build build --clean-first |
| ctest (3 Fortran unit/integration tests) | 0.88 s, 3/3 pass | ctest |
| run neo_rt.x driftorbit (s=0.5 base case) | 0.16 s wall | /usr/bin/time -p |

Golden gate: test/golden_record compares 5 output files at rtol=1e-8, atol=1e-15.
ODE: DVODE Adams (method_flag=10), rtol=1e-9, atol=1e-10, root-finding events (nevents=2).
Golden-path external deps actually used: spline (111 lines) + DVODE. No BLAS/LAPACK/NetCDF.
