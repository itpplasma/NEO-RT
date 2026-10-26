#!/usr/bin/env bash
#SBATCH --job-name=raw-oneform
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/home/ert/runs/claude-recovery-20261002/rawdrive/tiny-oneform-v2/slurm-%j.log
set -euo pipefail
cd /home/ert/runs/claude-recovery-20261002/rawdrive/tiny-oneform-v2
test "$(sha256sum /home/ert/runs/claude-recovery-20261002-resources/bounded-run-v2.py | cut -d' ' -f1)" = \
    49d12f9c6c954d5e48f6f1b088f3f97e5d71a41f59db3d64bbd5fe6185c3c262
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
    python3 /home/ert/runs/claude-recovery-20261002-resources/bounded-run-v2.py \
    --cpu-use 1 --memory-gib 1 --seconds 180 -- bash ./payload.sh
