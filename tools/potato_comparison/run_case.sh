#!/bin/bash
# Array tasks 0..9: POTATO; 10/11: NEO-RT profiles at s=.09/.36.
set -euo pipefail
TASK=$1
ROOT=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-/home/ert/iter_tc24_runs/_venv_ntv_py312/bin/python}
NCPU=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$ROOT/out"
if [ "$TASK" -lt 10 ]; then
    ID=$("$PY" -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["cases"][int(sys.argv[2])]["id"])' \
        "$ROOT/inputs/scan.json" "$TASK")
    W=$ROOT/out/p$ID
    TARGET=potato.x
    SUBDIR=POTATO
else
    SEL=$([ "$TASK" = 10 ] && echo 0.09 || echo 0.36)
    W=$ROOT/out/neort_s$SEL
    TARGET=neo_rt.x
    SUBDIR=.
fi
if [ -e "$W" ]; then echo "Output already exists: $W" >&2; exit 2; fi
mkdir -p "$W/src"
trap 'rc=$?; if [ "$rc" != 0 ]; then echo "rc $rc" > "$W/FAILED"; fi' EXIT
tar -xzf "$ROOT/neort_source.tgz" -C "$W/src" --strip-components=1
cmake -S "$W/src/$SUBDIR" -B "$W/build" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    > "$W/configure.log" 2>&1
cmake --build "$W/build" -j "$NCPU" --target "$TARGET" > "$W/build.log" 2>&1
EXE=$W/build/$TARGET
if [ "$TASK" -lt 10 ]; then
    "$PY" -I "$ROOT/make_cases.py" "$ROOT" potato "$ID" "$W/run"
    cd "$W/run"
    export OMP_NUM_THREADS=$NCPU OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
    export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    # Retain failures through the LCFS-noise filter, and keep the dropped count.
    /usr/bin/time -v "$EXE" 2> stderr.log | \
        awk '/magfie: point is outside the LCFS/ {n++; next} {print} END {print n+0 > "lcfs_lines_dropped.txt"}' \
        > stdout.log
    test -s integral_torque.dat
    test -s boxcounted_torque.dat
else
    "$PY" -I -c '
import json, sys
for c in json.load(open(sys.argv[1]))["cases"]:
    if c["s"] == float(sys.argv[2]) and c["kind"] != "energy_resolution":
        for r in c["rho_grid"]:
            for m in ("thin", "dph", "radial", "both"):
                print(c["id"], f"{r:.2f}", m)
' "$ROOT/inputs/scan.json" "$SEL" > "$W/list.txt"
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    run1() {
        local d=$ROOT/out/neort/$1_$2_$3
        "$PY" -I "$ROOT/make_cases.py" "$ROOT" neort "$1" "$d" "$2" "$3" && \
            (cd "$d" && "$EXE" driftorbit > stdout.log 2>&1) && \
            test -s "$d/driftorbit_torque.out"
    }
    export -f run1
    export ROOT PY EXE
    xargs -P "$NCPU" -L 1 bash -c 'run1 "$0" "$1" "$2"' < "$W/list.txt"
fi
echo done > "$W/DONE"
