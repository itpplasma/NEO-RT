#!/bin/bash
#SBATCH --job-name=tc24-gpec-kinetic
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --time=08:00:00

set -euo pipefail

campaign_root="${CAMPAIGN_ROOT:?CAMPAIGN_ROOT is required}"
source_root="${campaign_root}/source"
dependency_root="${campaign_root}/deps/usr"
netcdf_fortran_root="${dependency_root}/lib/x86_64-linux-gnu"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export LD_LIBRARY_PATH="${netcdf_fortran_root}:${LD_LIBRARY_PATH:-}"

make -C "${source_root}/install" -j"${OMP_NUM_THREADS}" \
    pentrc dcon coil gpec \
    FC=gfortran \
    CC=gcc \
    LAPACKHOME=/usr/lib/x86_64-linux-gnu \
    NETCDFHOME="${netcdf_fortran_root}" \
    NETCDFINC="${dependency_root}/include"

run_case() {
    case_name="$1"
    input_root="${campaign_root}/input-${case_name}"
    run_root="${campaign_root}/run-${case_name}"

    mkdir -p "${run_root}"
    cp \
        "${input_root}/coil.in" \
        "${input_root}/dcon.in" \
        "${input_root}/equil.in" \
        "${input_root}/gpec.in" \
        "${input_root}/pentrc.in" \
        "${input_root}/profiles.kin" \
        "${input_root}/vac.in" \
        "${run_root}/"

    env --chdir="${run_root}" \
        /usr/bin/time -v "${source_root}/dcon/dcon" \
        >"${run_root}/dcon.log" 2>"${run_root}/dcon.resource.log"

    env --chdir="${run_root}" \
        /usr/bin/time -v "${source_root}/gpec/gpec" \
        >"${run_root}/gpec.log" 2>"${run_root}/gpec.resource.log"

    sha256sum \
        "${source_root}/dcon/dcon" \
        "${source_root}/gpec/gpec" \
        "${run_root}/coil.in" \
        "${run_root}/dcon.in" \
        "${run_root}/equil.in" \
        "${run_root}/gpec.in" \
        "${run_root}/pentrc.in" \
        "${run_root}/profiles.kin" \
        "${run_root}/vac.in" \
        "${run_root}/dcon.log" \
        "${run_root}/gpec.log" \
        >"${run_root}/SHA256SUMS"
}

run_case ideal
run_case kinetic
