#!/usr/bin/env python3
"""Adaptive velocity integral of the collisionless resonance components.

Circular tokamak of examples/circ.bc, s = 0.5, n = 3, delta B/B = 1e-3,
M_t = 0.01, full resonant transport. Folds of the resonance lines make the
velocity integrand singular like 1/sqrt(u - u_fold); there the vsteps-node
midpoint rule converges like sqrt(du) and erratically: D11 is 6.36e-4,
7.01e-4, 8.69e-4 and 9.87e-4 at 512, 2048, 8192 and 32768 nodes.

The adaptive Gauss-Kronrod integral (vtol > 0) must agree with itself over
two decades of tolerance, and with the 32768-node midpoint value within its
remaining erratic error.
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).resolve().parent
EXAMPLES = TEST_DIR.parent.parent / "examples"
CASE_NAME = "adaptive"

# Midpoint rule with 32768 velocity nodes on this case (154 s): its error is
# erratic in sign, about 3 % at this resolution.
D11_MIDPOINT_32768 = 9.8734994302942061e-04

INPUT = """&params
    s = 0.5
    M_t = 0.01
    qs = 1.0
    ms = 2.014
    vth = 40000000.0
    epsmn = 0.001
    m0 = 0
    mph = 3
    supban = .false.
    magdrift = .true.
    nopassing = .false.
    noshear = .false.
    pertfile = .false.
    nonlin = .false.
    comptorque = .false.
    inp_swi = 8
    vsteps = 16
    vtol = {vtol}
    log_level = -1
/
"""


def find_executable() -> Path:
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


def d11(executable: Path, work_dir: Path, vtol: float) -> float:
    work_dir.mkdir()
    (work_dir / "in_file").symlink_to((EXAMPLES / "circ.bc").resolve())
    (work_dir / f"{CASE_NAME}.in").write_text(INPUT.format(vtol=vtol))
    result = subprocess.run([str(executable), CASE_NAME], cwd=work_dir,
                            capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"neo_rt.x failed:\n{result.stdout}\n{result.stderr}")
    return float(np.loadtxt(work_dir / f"{CASE_NAME}.out")[4])


def test_adaptive_velocity_integral_converges(tmp_path):
    executable = find_executable()
    values = {vtol: d11(executable, tmp_path / f"vtol{vtol:g}", vtol)
              for vtol in (1e-3, 1e-4, 1e-5)}
    reference = values[1e-5]
    for vtol, value in values.items():
        assert abs(value / reference - 1) < 2e-3, (vtol, value, reference)
    assert abs(reference / D11_MIDPOINT_32768 - 1) < 0.04, reference


# Ripple-plateau case of test/ripple_plateau with the plasma of examples/base:
# no magnetic drift, harmonics up to |m| = 41. Pitch and speed adaptive
# integrals at vtol = 1e-5 agree on D11 to 2e-6; the 256-node midpoint rule
# is off by 3e-4.
RIPPLE_DIR = TEST_DIR.parent / "ripple_plateau"
D11_RIPPLE = 2.1583872e-02


def test_adaptive_velocity_shares_error_budget_over_harmonics(tmp_path):
    """Strongly cancelling high harmonics must not be refined to their own
    relative tolerance: this case then took minutes instead of seconds."""
    executable = find_executable()
    work_dir = tmp_path / "ripple"
    work_dir.mkdir()
    text = (RIPPLE_DIR / "driftorbit_test.in").read_text()
    text = text.replace("log_level = -1", "vtol = 1e-3\n    log_level = -1")
    (work_dir / "driftorbit_test.in").write_text(text)
    (work_dir / "in_file").symlink_to((RIPPLE_DIR / "in_file").resolve())
    for name in ("plasma.in", "profile.in"):
        (work_dir / name).symlink_to((EXAMPLES / "base" / name).resolve())
    result = subprocess.run([str(executable), "driftorbit_test"], cwd=work_dir,
                            capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        pytest.fail(f"neo_rt.x failed:\n{result.stdout}\n{result.stderr}")
    value = float(np.loadtxt(work_dir / "driftorbit_test.out")[4])
    assert abs(value / D11_RIPPLE - 1) < 1e-5, value
