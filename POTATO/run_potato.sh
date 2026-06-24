#!/bin/bash
# Run potato.x with one OpenMP thread per physical core.
#
# The resonance solver (grid build + per-mode root search) is floating-point
# bound, so the SMT sibling threads contend for the FPU instead of adding
# throughput. On a 16-core/32-thread host the rho_pol=0.99 #30835 sweep takes
# ~148 s on 16 threads bound to cores, but ~187 s on all 32 logical threads when
# idle and 400-600 s under any concurrent load, because 32 threads oversubscribe
# the 16 cores. One thread per physical core is both fastest and load-robust.
#
# Export OMP_NUM_THREADS / OMP_PLACES / OMP_PROC_BIND before calling to override.
set -eu

if [ -z "${OMP_NUM_THREADS:-}" ]; then
  cores=$(lscpu -p=Core 2>/dev/null | grep -v '^#' | sort -u | wc -l)
  case "$cores" in ''|*[!0-9]*) cores=$(nproc) ;; esac
  [ "$cores" -ge 1 ] || cores=$(nproc)
  export OMP_NUM_THREADS="$cores"
fi
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

exe="${POTATO_EXE:-$(dirname "$0")/build/potato.x}"
exec "$exe" "$@"
