#!/usr/bin/env bash
# Execute inside a scheduler allocation; fo owns target building and execution.
set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" && -z "${PBS_JOBID:-}" ]]; then
    echo 'Run line-drive diagnostics inside a Slurm or PBS allocation.' >&2
    exit 2
fi
if [[ $# -ne 1 ]]; then
    echo 'Usage: run_diagnostics.sh OUTPUT_DIRECTORY' >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/../.." && pwd)
mkdir -p -- "$1"
output_dir=$(cd -- "$1" && pwd)
export LINE_DRIVE_DIAGNOSTICS_DIR="$output_dir"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
cd -- "$repo_dir"
fo exec --cwd "$repo_dir/examples/base" diagnose_line_drive.x
"${PYTHON:-python3}" "$script_dir/plot_diagnostics.py" "$output_dir"
"${PYTHON:-python3}" "$script_dir/plot_diagnostics.py" "$output_dir" --grayscale
