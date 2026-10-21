#!/bin/bash
# The same native test selection runs in a build container or a heavy Slurm tier.
set -euo pipefail
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$src"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    root=$(dirname "$src")
    tier=$(cat "$root/tier.txt")
    . "$root/task.env"
    . /home/ert/opt/kin6d-env.sh
    export PATH=/home/ert/.local/bin:$PATH
else
    tier=${1:?usage: ci/run.sh TIER}
    root=$(mktemp -d "${RUNNER_TEMP:-/tmp}/native-ci.XXXXXX")
    export FC=/usr/bin/gfortran CC=/usr/bin/gcc CXX=/usr/bin/g++
    export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
fi
unset DISPLAY
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
if [ -d "$root/mirrors" ]; then
    export GIT_CONFIG_GLOBAL=$root/gitconfig
    : > "$GIT_CONFIG_GLOBAL"
    git config --global protocol.file.allow always
    shopt -s nullglob
    for mirror in "$root"/mirrors/*/*.git; do
        slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
        git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
    done
fi
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export SCI_CPUS=${SCI_CPUS:-8}
export CMAKE_BUILD_PARALLEL_LEVEL=$SCI_CPUS CTEST_PARALLEL_LEVEL=$SCI_CPUS
# The native parallel oracle requires at least two OpenMP workers.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP" ci-out
{ hostname; /usr/bin/gfortran --version | head -1; git rev-parse HEAD; } | tee ci-out/toolchain.txt
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
case "$tier" in
  test-standalone-off)
    export SCI_RUN_GOLDEN=${SCI_RUN_GOLDEN:-false}
    echo 'Check OFF configure diagnostic without fetching dependencies'
    (
if cmake -S . -B build-off -DUSE_STANDALONE=OFF > configure-off.log 2>&1; then
  cat configure-off.log
  echo "Expected USE_STANDALONE=OFF to fail"
  exit 1
fi
cat configure-off.log
grep -F "USE_STANDALONE=OFF is currently unsupported" configure-off.log
test ! -d build-off/OBJS/_deps

    )
    load_job_environment
    ;;
  test-build)
    python3 -m venv "$root/venv"
    . "$root/venv/bin/activate"
    export SCI_RUN_GOLDEN=${SCI_RUN_GOLDEN:-false}
    echo 'Install Python test dependencies'
    (
pip install pytest numpy h5py f90nml scipy

    )
    load_job_environment
    echo 'Build code'
    (
make CONFIG=${SCI_MATRIX_CONFIG}${SCI_LIBNEO_REF:+ LIBNEO_REF=$SCI_LIBNEO_REF}

    )
    load_job_environment
    echo 'Run ctest'
    (
make ctest

    )
    load_job_environment
    echo 'Run pytest (all tests for Fast)'
    if [ "$SCI_RUN_GOLDEN" = true ] && [ "${SCI_MATRIX_CONFIG}" = Fast ]; then
    (
make pytest

    )
    fi
    load_job_environment
    echo 'Run pytest (skip golden record for Debug)'
    if [ "$SCI_RUN_GOLDEN" = true ] && [ "${SCI_MATRIX_CONFIG}" != Fast ]; then
    (
pytest -m "not fast_only"

    )
    fi
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
