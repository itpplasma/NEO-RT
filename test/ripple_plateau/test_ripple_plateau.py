#!/usr/bin/env python3
"""
Ripple plateau test for NEO-RT.

Tests that the D12/D11 ratio matches the theoretical value of 3.0
for ripple plateau transport.
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).resolve().parent
CASE_NAME = "driftorbit_test"
EPSMN = 1e-3
EXPECTED_D12_D11_RATIO = 3.0
RATIO_TOLERANCE = 0.08
D11_TOLERANCE = 0.05


def find_executable() -> Path:
    """Find neo_rt.x executable."""
    candidates = [
        TEST_DIR.parent.parent / "build" / "neo_rt.x",
        Path(os.environ.get("NEO_RT_EXE", "")),
    ]

    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate

    result = shutil.which("neo_rt.x")
    if result:
        return Path(result)

    pytest.skip("neo_rt.x executable not found")


@pytest.fixture(scope="module")
def executable() -> Path:
    """Fixture providing path to neo_rt.x."""
    return find_executable()


def setup_work_dir(work_dir: Path) -> None:
    """Set up symlinks to required input files."""
    input_file = TEST_DIR / f"{CASE_NAME}.in"
    if not input_file.exists():
        pytest.fail(f"Required input file missing: {input_file}")
    (work_dir / f"{CASE_NAME}.in").symlink_to(input_file)

    in_file_src = TEST_DIR / "in_file"
    if not in_file_src.exists():
        pytest.fail(f"Required input file missing: {in_file_src}")
    (work_dir / "in_file").symlink_to(in_file_src.resolve())


def run_neort(executable: Path, work_dir: Path) -> None:
    """Run neo_rt.x for the test case."""
    result = subprocess.run(
        [str(executable), CASE_NAME],
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"neo_rt.x failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )


def load_output(work_dir: Path) -> np.ndarray:
    """Load main output file."""
    return np.loadtxt(work_dir / f"{CASE_NAME}.out")


def parse_drp_from_magfie_param(work_dir: Path) -> float:
    """Extract Drp value from magfie_param output file."""
    param_file = work_dir / f"{CASE_NAME}_magfie_param.out"
    with open(param_file, "r") as f:
        for line in f:
            if "Drp" in line:
                return float(line.split()[-1])
    pytest.fail(f"Drp not found in {param_file}")


def test_ripple_plateau(executable: Path, tmp_path: Path) -> None:
    """Test that D12/D11 ratio matches theoretical value of 3.0."""
    work_dir = tmp_path / "ripple_plateau"
    work_dir.mkdir()
    setup_work_dir(work_dir)

    run_neort(executable, work_dir)

    data = load_output(work_dir)
    D11 = data[4]
    D12 = data[8]
    Drp = parse_drp_from_magfie_param(work_dir)

    # Test D11/eps^2 against reference Drp
    D11_normalized = D11 / EPSMN**2
    d11_relative_error = abs(D11_normalized - Drp) / Drp
    assert d11_relative_error < D11_TOLERANCE, (
        f"D11/eps^2 = {D11_normalized:.3e} deviates from reference "
        f"Drp = {Drp:.3e} by {d11_relative_error*100:.1f}% "
        f"(tolerance: {D11_TOLERANCE*100:.0f}%)"
    )

    # Test D12/D11 ratio against theoretical value
    actual_ratio = D12 / D11
    ratio_relative_error = abs(actual_ratio - EXPECTED_D12_D11_RATIO) / EXPECTED_D12_D11_RATIO
    assert ratio_relative_error < RATIO_TOLERANCE, (
        f"D12/D11 ratio {actual_ratio:.3e} deviates from expected "
        f"{EXPECTED_D12_D11_RATIO:.3e} by {ratio_relative_error*100:.1f}% "
        f"(tolerance: {RATIO_TOLERANCE*100:.0f}%)"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
