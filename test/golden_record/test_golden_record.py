#!/usr/bin/env python3
"""
Golden record regression tests for NEO-RT.

Runs neo_rt.x for each test case and compares outputs against golden record data
stored in golden.h5.

Test cases are discovered from the HDF5 file groups. Case names follow the pattern
'0pXXX' where XXX represents the s-value (e.g., '0p100' -> s=0.1).

The golden.h5 file is automatically generated from the main branch if it does not
exist. This ensures tests compare against a known-good reference.
"""

import os
import shutil
import subprocess
from pathlib import Path

import f90nml
import h5py
import numpy as np
import pytest

from ensure_golden import ensure_golden

RTOL = 1e-8
ATOL = 1e-15

GOLDEN_DIR = Path(__file__).resolve().parent
INPUT_DIR = GOLDEN_DIR / "input"

# Ensure golden.h5 exists before test discovery
GOLDEN_H5 = ensure_golden()


def s_from_case_name(case_name: str) -> float:
    """Convert case name to s value: '0p100' -> 0.1"""
    return float(case_name.replace("p", "."))


def case_name_from_s(s: float) -> str:
    """Convert s value to case name: 0.1 -> '0p100'"""
    return f"{s:.3f}".replace(".", "p")


def get_test_cases() -> list[str]:
    """Discover test cases from HDF5 groups."""
    with h5py.File(GOLDEN_H5, "r") as f:
        return sorted([key for key in f.keys() if key.startswith("0p")])


def find_executable() -> Path:
    """Find neo_rt.x executable."""
    # Try common locations relative to the test directory
    candidates = [
        GOLDEN_DIR.parent.parent / "build" / "neo_rt.x",
        Path(os.environ.get("NEO_RT_EXE", "")),
    ]

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    # Try finding in PATH
    result = shutil.which("neo_rt.x")
    if result:
        return Path(result)

    pytest.skip("neo_rt.x executable not found")


@pytest.fixture(scope="module")
def executable() -> Path:
    """Fixture providing path to neo_rt.x."""
    return find_executable()


def generate_input_file(work_dir: Path, case: str) -> Path:
    """Generate input file for a test case from template."""

    template_path = INPUT_DIR / "template.in"
    template = f90nml.read(template_path)

    s = s_from_case_name(case)
    template["params"]["s"] = s

    input_path = work_dir / f"{case}.in"
    template.write(input_path, force=True)

    return input_path


def setup_work_dir(work_dir: Path) -> None:
    """Set up symlinks to required input files."""
    for filename in ["in_file", "in_file_pert", "plasma.in", "profile.in"]:
        src = INPUT_DIR / filename
        if not src.exists():
            pytest.fail(
                f"Required input file missing: {src}\n"
                f"Context:\n\t{GOLDEN_DIR=}\n\t{GOLDEN_H5=}\n\t{INPUT_DIR=}"
            )
        dst = work_dir / filename
        if not dst.exists():
            dst.symlink_to(src)


def run_neort(executable: Path, work_dir: Path, case: str) -> None:
    """Run neo_rt.x for a test case."""
    result = subprocess.run(
        [str(executable), case],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"neo_rt.x failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )


def load_output(work_dir: Path, case: str) -> np.ndarray:
    """Load main output file."""
    return np.loadtxt(work_dir / f"{case}.out", skiprows=1)


def load_torque(work_dir: Path, case: str) -> np.ndarray:
    """Load torque output file."""
    return np.loadtxt(work_dir / f"{case}_torque.out", skiprows=1)


def load_integral(work_dir: Path, case: str) -> np.ndarray:
    """Load integral output file."""
    return np.loadtxt(work_dir / f"{case}_integral.out")


def load_torque_integral(work_dir: Path, case: str) -> np.ndarray:
    """Load torque integral output file."""
    return np.loadtxt(work_dir / f"{case}_torque_integral.out")


def load_magfie(work_dir: Path, case: str) -> np.ndarray:
    """Load magnetic field data."""
    return np.loadtxt(work_dir / f"{case}_magfie.out")


@pytest.mark.parametrize("case", get_test_cases())
def test_golden(case: str, executable: Path, tmp_path: Path) -> None:
    """Test that neo_rt.x output matches golden record."""
    # Set up work directory
    work_dir = tmp_path / case
    work_dir.mkdir()
    setup_work_dir(work_dir)

    # Generate input and run
    generate_input_file(work_dir, case)
    run_neort(executable, work_dir, case)

    # Load golden data
    with h5py.File(GOLDEN_H5, "r") as f:
        grp = f[case]
        golden_output = grp["output"][:]
        golden_torque = grp["torque"][:]
        golden_integral = grp["integral"][:]
        golden_torque_integral = grp["torque_integral"][:]
        golden_magfie = grp["magfie"][:]

    # Load actual outputs
    actual_output = load_output(work_dir, case)
    actual_torque = load_torque(work_dir, case)
    actual_integral = load_integral(work_dir, case)
    actual_torque_integral = load_torque_integral(work_dir, case)
    actual_magfie = load_magfie(work_dir, case)

    # Compare all outputs
    assert np.allclose(
        actual_output, golden_output, rtol=RTOL, atol=ATOL
    ), f"Output mismatch for {case}"
    assert np.allclose(
        actual_torque, golden_torque, rtol=RTOL, atol=ATOL
    ), f"Torque mismatch for {case}"
    assert np.allclose(
        actual_integral, golden_integral, rtol=RTOL, atol=ATOL
    ), f"Integral mismatch for {case}"
    assert np.allclose(
        actual_torque_integral, golden_torque_integral, rtol=RTOL, atol=ATOL
    ), f"Torque integral mismatch for {case}"
    assert np.allclose(
        actual_magfie, golden_magfie, rtol=RTOL, atol=ATOL
    ), f"Magfie mismatch for {case}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
