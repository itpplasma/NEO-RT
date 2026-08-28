#!/bin/bash
# Deterministic golden-record gate for any NEO-RT build (Fortran or a port).
# Usage: run_gate.sh /path/to/neo_rt_executable
#
# The gate compares 5 output files per case against test/golden_record/golden.h5
# at rtol=1e-8. It MUST run single-threaded: the Fortran code is bit-exact
# run-to-run only under OMP_NUM_THREADS=1; under OpenMP, results vary at the
# ~1e-4 level (run-configuration dependent), so the multithreaded gate is not
# meaningful for cross-implementation comparison.
#
# golden.h5 is gitignored and persists in the working tree across branch
# checkouts, so the identical reference gates every port branch. Regenerate it
# from the Fortran build with:
#   OMP_NUM_THREADS=1 python regenerate_golden.py /path/to/fortran/neo_rt.x
set -e
EXE="${1:?usage: run_gate.sh /path/to/neo_rt_executable}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO/test"
OMP_NUM_THREADS=1 NEO_RT_EXE="$EXE" \
  uv run --quiet --with f90nml --with h5py --with numpy --with pytest \
  python3 -m pytest golden_record/test_golden_record.py -q
