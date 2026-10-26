#!/usr/bin/env bash
set -euo pipefail
cd /home/ert/runs/claude-recovery-20261002/rawdrive/tiny-fullgc
export OMPI_MCA_plm=isolated OMPI_MCA_hwloc_base_binding_policy=none
export FC=gfortran CC=gcc FO_JOBS=1 CMAKE_BUILD_PARALLEL_LEVEL=1
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
export CTEST_PARALLEL_LEVEL=1 FO_BUILD_TIMEOUT=100 FO_TEST_TIMEOUT=120
export FO_CMAKE_ARGS='-DCMAKE_BUILD_TYPE=Debug'
test -s inputs.sha256
sha256sum --check inputs.sha256 > preflight.log
/home/ert/.local/bin/fo info > fo-info.log 2>&1
/home/ert/.local/bin/fo build > build.log 2>&1
/home/ert/.local/bin/fo test test_raw_oneform > focused-oneform.log 2>&1
cp build/Testing/Temporary/LastTest.log oracle-oneform.log
/home/ert/.local/bin/fo test test_raw_gc > focused-gc.log 2>&1
cp build/Testing/Temporary/LastTest.log oracle-gc.log
printf 'RAW_FULL_GC_AND_GAUGE_ORACLE_PASS\n' > result.txt
