#!/usr/bin/env python3
"""
Regenerate golden record data for NEO-RT regression tests.

This script runs neo_rt.x for 9 evenly spaced s values (0.1 to 0.9),
collects the output files, and creates a new golden.h5.

Usage:
    python regenerate_golden.py [path_to_neo_rt.x]

If no executable path is provided, defaults to ../../build/neo_rt.x
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import f90nml
import h5py
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR = SCRIPT_DIR / "input"
GOLDEN_H5 = SCRIPT_DIR / "golden.h5"


def get_runs(nruns: int = 10) -> list[str]:
    """Generate case names for evenly spaced s values from 0.1 to 0.9."""
    runs = []
    for i in range(1, nruns):
        s = i / nruns
        run = f"{s:.3f}".replace(".", "p")
        runs.append(run)
    return runs


def s_from_case_name(case_name: str) -> float:
    """Convert case name to s value: '0p100' -> 0.1"""
    return float(case_name.replace("p", "."))


def load_output(input_dir: Path, case: str) -> np.ndarray:
    """Load main output file (diffusion coefficients)."""
    return np.loadtxt(input_dir / f"{case}.out", skiprows=1)


def load_torque(input_dir: Path, case: str) -> np.ndarray:
    """Load torque output file."""
    return np.loadtxt(input_dir / f"{case}_torque.out", skiprows=1)


def load_integral(input_dir: Path, case: str) -> np.ndarray:
    """Load integral output file (transport by harmonic)."""
    return np.loadtxt(input_dir / f"{case}_integral.out")


def load_torque_integral(input_dir: Path, case: str) -> np.ndarray:
    """Load torque integral output file."""
    return np.loadtxt(input_dir / f"{case}_torque_integral.out")


def load_magfie(input_dir: Path, case: str) -> np.ndarray:
    """Load magnetic field profile data."""
    return np.loadtxt(input_dir / f"{case}_magfie.out")


def load_magfie_param(input_dir: Path, case: str) -> dict[str, float]:
    """Load magnetic field parameters as key-value pairs."""
    path = input_dir / f"{case}_magfie_param.out"
    params = {}
    with open(path) as f:
        for line in f:
            if ":" in line and "=" in line:
                parts = line.split(":")
                if len(parts) >= 2:
                    key_val = parts[1].split("=")
                    if len(key_val) == 2:
                        key = key_val[0].strip()
                        val_str = key_val[1].strip()
                        try:
                            if val_str in ("T", ".true.", ".TRUE."):
                                params[key] = 1.0
                            elif val_str in ("F", ".false.", ".FALSE."):
                                params[key] = 0.0
                            else:
                                params[key] = float(
                                    val_str.replace("D", "E").replace("d", "e")
                                )
                        except ValueError:
                            pass
    return params


def create_golden_h5(input_dir: Path, output_h5: Path, cases: list[str]) -> None:
    """Create HDF5 file with all golden record data."""
    print(f"Creating {output_h5} with {len(cases)} test cases...")

    with h5py.File(output_h5, "w") as f:
        f.attrs["description"] = "Golden record data for NEO-RT regression tests"
        f.attrs["n_cases"] = len(cases)
        f.attrs["output_columns"] = (
            "M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12"
        )
        f.attrs["torque_columns"] = "s dVds M_t Tco Tctr Tt"
        f.attrs["integral_columns"] = (
            "M_t mth Dco_trap Dco_pass Dco Dctr_trap Dctr_pass "
            "Dctr Dt_trap Dt_pass Dt etaco etactr etatp etadt"
        )
        f.attrs["torque_integral_columns"] = "mth Tco Tctr Tt"

        for case in cases:
            print(f"  Processing {case}...")
            grp = f.create_group(case)
            grp.attrs["s"] = s_from_case_name(case)
            grp.create_dataset("output", data=load_output(input_dir, case))
            grp.create_dataset("torque", data=load_torque(input_dir, case))
            grp.create_dataset("integral", data=load_integral(input_dir, case))
            grp.create_dataset(
                "torque_integral", data=load_torque_integral(input_dir, case)
            )
            grp.create_dataset("magfie", data=load_magfie(input_dir, case))

            for key, val in load_magfie_param(input_dir, case).items():
                grp.attrs[f"magfie_{key}"] = val


def main() -> None:
    if len(sys.argv) > 1:
        neort_exe = Path(sys.argv[1]).resolve()
    else:
        neort_exe = SCRIPT_DIR.parent.parent / "build" / "neo_rt.x"

    if not neort_exe.exists():
        print(f"Error: Executable not found: {neort_exe}")
        sys.exit(1)

    print(f"Using executable: {neort_exe}")

    runs = get_runs(10)
    print(f"Will generate golden records for: {runs}")

    template = f90nml.read(INPUT_DIR / "template.in")

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        print(f"Working directory: {out_dir}")

        for filename in ["in_file", "in_file_pert", "plasma.in", "profile.in"]:
            src = INPUT_DIR / filename
            if not src.exists():
                print(f"Error: Required input file missing: {src}")
                sys.exit(1)
            (out_dir / filename).symlink_to(src)

        for run in runs:
            nml = template.copy()
            s = float(run.replace("p", "."))
            nml["params"]["s"] = s
            nml.write(out_dir / f"{run}.in", force=True)

        n_parallel: int = os.cpu_count() or 4
        processes: list[subprocess.Popen] = []

        for run in runs:
            if len(processes) >= n_parallel:
                processes[0].wait()
                processes.pop(0)
            print(f"Starting NEO-RT for run {run}...")
            p = subprocess.Popen(
                [str(neort_exe), run],
                cwd=out_dir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            processes.append(p)

        for p in processes:
            p.wait()

        print("All runs completed.")
        create_golden_h5(out_dir, GOLDEN_H5, runs)

    print("Done!")


if __name__ == "__main__":
    main()
