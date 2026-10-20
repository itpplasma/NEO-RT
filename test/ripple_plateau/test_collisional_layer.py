#!/usr/bin/env python3
"""Integration tests for the optional collisional boundary-layer factor.

Covers:
  (1) collisional_layer=.false. is bit-identical to the flag being absent
      (default-off back-compatibility on the ripple-plateau case);
  (2) with real (energy-dependent) collisionality and an interior resonance,
      the transport integral approaches the collisionless delta result as
      collisionality -> 0; large density exercises the independently validated
      narrow-domain cubic suppression limit.

The physics tests reuse the ripple-plateau magnetic field but load the real
plasma/rotation profiles from examples/base so that the deflection frequency is
nonzero.  Collisionality is scanned by scaling the plasma density (nu_d ~ n,
layer width delta ~ n^(1/3)); comparing on/off at the same density isolates the
layer multiplier, which is independent of density in the D coefficients.
"""

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent.parent
BASE_DIR = REPO_ROOT / "examples" / "base"
CASE_NAME = "driftorbit_test"


def find_executable() -> Path:
    candidates = [
        Path(os.environ.get("NEO_RT_EXE", "")),
        REPO_ROOT / "build" / "neo_rt.x",
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
    return find_executable()


def _write_input(work_dir: Path, extra_lines: str, nopassing: bool) -> None:
    text = (TEST_DIR / f"{CASE_NAME}.in").read_text()
    if nopassing:
        text = text.replace("nopassing = .false.", "nopassing = .true.")
    text = text.replace("log_level = -1", f"{extra_lines}log_level = -1")
    (work_dir / f"{CASE_NAME}.in").write_text(text)
    (work_dir / "in_file").symlink_to((TEST_DIR / "in_file").resolve())


def _scale_density(dst: Path, factor: float) -> None:
    """Copy examples/base plasma.in, scaling the two ion density columns."""
    lines = (BASE_DIR / "plasma.in").read_text().splitlines()
    out = []
    for i, line in enumerate(lines):
        if i >= 3 and line.strip():
            cols = line.split()
            cols[1] = repr(float(cols[1]) * factor)
            cols[2] = repr(float(cols[2]) * factor)
            out.append("   " + "   ".join(cols))
        else:
            out.append(line)
    dst.write_text("\n".join(out) + "\n")


def _run(executable: Path, work_dir: Path) -> None:
    result = subprocess.run(
        [str(executable), CASE_NAME], cwd=work_dir, capture_output=True, text=True,
        timeout=120
    )
    if result.returncode != 0:
        pytest.fail(f"neo_rt.x failed:\nstdout: {result.stdout}\nstderr: {result.stderr}")


def _d11(work_dir: Path) -> float:
    return float(np.loadtxt(work_dir / f"{CASE_NAME}.out")[4])


def _run_case(
    executable: Path, work_dir: Path, flag: str, factor: float, nopassing: bool
) -> float:
    work_dir.mkdir()
    extra = f"collisional_layer = {flag}\n    " if flag else ""
    _write_input(work_dir, extra, nopassing)
    _scale_density(work_dir / "plasma.in", factor)
    (work_dir / "profile.in").symlink_to((BASE_DIR / "profile.in").resolve())
    _run(executable, work_dir)
    return _d11(work_dir)


def test_off_is_bit_identical(executable: Path, tmp_path: Path) -> None:
    """collisional_layer=.false. reproduces the flag-absent output byte-for-byte."""
    baseline = tmp_path / "baseline"
    explicit_off = tmp_path / "explicit_off"
    for work_dir, flag in ((baseline, ""), (explicit_off, ".false.")):
        work_dir.mkdir()
        _write_input(work_dir, f"collisional_layer = {flag}\n    " if flag else "", False)
        _scale_density(work_dir / "plasma.in", 1.0)
        (work_dir / "profile.in").symlink_to((BASE_DIR / "profile.in").resolve())
        _run(executable, work_dir)

    a = (baseline / f"{CASE_NAME}.out").read_bytes()
    b = (explicit_off / f"{CASE_NAME}.out").read_bytes()
    assert a == b, "collisional_layer=.false. is not bit-identical to the default"


@pytest.mark.fast_only
def test_collisionless_limit(executable: Path, tmp_path: Path) -> None:
    """Vanishing collisionality recovers the delta result within 2%."""
    d_off = _run_case(executable, tmp_path / "off", ".false.", 1e-9, nopassing=True)
    d_on = _run_case(executable, tmp_path / "on", ".true.", 1e-9, nopassing=True)
    ratio = d_on / d_off
    assert abs(ratio - 1.0) < 0.02, f"collisionless limit ratio {ratio:.4f} not within 2%"


@pytest.mark.fast_only
def test_high_collisionality_suppression(executable: Path, tmp_path: Path) -> None:
    """The validated narrow-domain limit suppresses large-density transport."""
    d_off = _run_case(executable, tmp_path / "off", ".false.", 1e9, nopassing=True)
    d_on = _run_case(executable, tmp_path / "on", ".true.", 1e9, nopassing=True)
    ratio = d_on / d_off
    assert ratio < 0.6, f"high-collisionality ratio {ratio:.4f} not suppressed below 0.6"


@pytest.mark.parametrize("extra, message", [
    ("krook_collisions=.true., krook_nu=1,", "cannot be combined"),
    ("", "requires positive plasma profiles"),
])
def test_rejects_undefined_layer(executable, tmp_path, extra, message):
    """The layer must not silently vanish under Krook or missing profiles."""
    _write_input(tmp_path, "collisional_layer=.true., " + extra, True)
    result = subprocess.run([str(executable), CASE_NAME], cwd=tmp_path,
                            capture_output=True, text=True, timeout=120)
    assert result.returncode != 0
    assert message in result.stderr


@pytest.mark.fast_only
def test_monoenergetic_layer(executable, tmp_path):
    """The layer multiplies monoenergetic D11 and its exact x*D11 D12."""
    data = []
    for flag in (".false.", ".true."):
        work = tmp_path / flag
        work.mkdir()
        _write_input(work, f"collisional_layer={flag}, monoenergetic_x=9, comptorque=.true., ", True)
        _scale_density(work / "plasma.in", 1e6)
        (work / "profile.in").symlink_to((BASE_DIR / "profile.in").resolve())
        ledger = work / "roots.txt"
        result = subprocess.run([str(executable), CASE_NAME], cwd=work,
            env=os.environ | {"NEORT_RESONANCE_LEDGER": str(ledger)},
            capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stderr
        contributions = {}
        for line in ledger.read_text().splitlines()[1:]:
            fields = line.split()
            assert fields[4] == 'trap'
            harmonic = int(fields[2])
            contributions[harmonic] = contributions.get(harmonic, 0.0) + float(fields[15])
        native = np.atleast_2d(np.loadtxt(work / f"{CASE_NAME}_torque_integral.out"))
        assert any(value != 0.0 for value in contributions.values())
        for harmonic in native:
            np.testing.assert_allclose(contributions.get(int(harmonic[0]), 0.0),
                                       harmonic[3], rtol=2e-13, atol=1e-280)
        row = np.loadtxt(work / f"{CASE_NAME}.out")
        assert row[4] > 0
        np.testing.assert_allclose(row[8], 9.0*row[4], rtol=2e-14)
        data.append(row[4])
    assert data[1] < 0.6*data[0]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
