#!/usr/bin/env python3
"""Public-safe golden-record gate for the POTATO OpenMP resonance solver.

Why this exists
---------------
The resonance root search was parallelized with OpenMP (per-mode loop in
integrate_class_resonances). One regression on that path -- privatizing the
get_rescond by-products psiast_res/taub_res/delphi_res as host-local privates,
which gfortran's dummy-procedure trampoline does not honor -- silently corrupted
the resonance weights and roots. That corruption is equilibrium-INDEPENDENT, so
a synthetic concentric-circle tokamak reproduces it without shipping any
experimental AUG #30835 g-file or kinetic profiles (which must not enter public
NEO-RT). The circular case is fully reproducible from gen_circular_eqdsk.py.

What it checks
--------------
1. potato.x runs the circular case to completion (exit 0).
2. The per-energy-slice resonance-line counts (lines of fort.31415 grouped by
   total energy, the first column) match golden_reslines.json within +-2%.
   The buggy commit 363fe03 yields per-slice counts roughly 0.6x the golden
   values, far outside the tolerance, so the gate fires on the regression.
3. Determinism: the sorted fort.31415 from OMP_NUM_THREADS=4 is byte-identical
   to OMP_NUM_THREADS=1 (the reslines are written from a critical section in
   thread-completion order, so only the sorted set is order-invariant).

Invoke directly with pytest, or via ctest (POTATO/CMakeLists.txt registers it as
the test `resonance_golden` when Python + numpy are available).
"""
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
GOLDEN = json.loads((HERE / "golden_reslines.json").read_text())

# Counts are small integers; allow cross-platform root jitter but stay far below
# the correct-vs-buggy gap (golden ~133/200 vs buggy ~82/120).
REL_TOL = 0.02
ABS_FLOOR = 3  # tolerate at least +-3 lines on the smallest slices

INPUT_FILES = (
    "potato.in",
    "field_divB0.inp",
    "convexwall.dat",
    "profile_poly.in",
    "circ.eqdsk",
    "bmod_n.dat",
)


def find_executable() -> Path:
    # POTATO_EXE wins so a caller can point the gate at a specific build
    # (e.g. comparing a candidate against a reference) without a rebuild.
    candidates = [
        Path(os.environ.get("POTATO_EXE", "")),
        HERE.parent.parent / "build" / "potato.x",
    ]
    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    found = shutil.which("potato.x")
    if found:
        return Path(found)
    pytest.skip("potato.x executable not found (build POTATO first)")


def ensure_inputs() -> None:
    """Regenerate the synthetic inputs if any are missing."""
    if all((HERE / f).exists() for f in INPUT_FILES):
        return
    subprocess.run(
        ["python3", str(HERE / "gen_circular_eqdsk.py")],
        cwd=HERE,
        check=True,
    )


def run_case(exe: Path, work_dir: Path, threads: int) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    for name in INPUT_FILES:
        shutil.copy(HERE / name, work_dir / name)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    result = subprocess.run(
        [str(exe)],
        cwd=work_dir,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-25:])
        pytest.fail(
            f"potato.x failed (OMP_NUM_THREADS={threads}, "
            f"rc={result.returncode}):\n{tail}\nstderr:\n{result.stderr}"
        )
    res = work_dir / "fort.31415"
    if not res.is_file():
        pytest.fail(f"fort.31415 not produced (OMP_NUM_THREADS={threads})")
    return res


def per_slice_counts(lines: list[str]) -> list[tuple[float, int]]:
    """Group fort.31415 lines by total energy (first column), sorted by energy.

    Returns [(toten, count), ...] so it can be compared against the golden
    slices without depending on the (thread-dependent) write order.
    """
    counts: dict[str, int] = {}
    keys: dict[str, float] = {}
    for line in lines:
        parts = line.split()
        if not parts:
            continue
        key = parts[0]
        counts[key] = counts.get(key, 0) + 1
        keys[key] = float(parts[0])
    return [(keys[k], counts[k]) for k in sorted(counts, key=lambda k: keys[k])]


@pytest.fixture(scope="module")
def reslines_runs(tmp_path_factory) -> dict[int, list[str]]:
    """Run the circular case once at OMP=1 and once at OMP=4; share the output.

    Module-scoped so potato.x is invoked exactly twice for the whole gate
    (~40 s total), not once per test.
    """
    ensure_inputs()
    exe = find_executable()
    base = tmp_path_factory.mktemp("resonance_golden")
    runs = {}
    for threads in (1, 4):
        res = run_case(exe, base / f"omp{threads}", threads=threads)
        runs[threads] = res.read_text().splitlines()
    return runs


def test_resonance_counts_match_golden(reslines_runs) -> None:
    """Per-slice resonance-line counts match the golden record within tolerance.

    Fires on commit 363fe03, whose corrupted weights drop the counts well below
    the +-2% band.
    """
    actual = per_slice_counts(reslines_runs[1])
    golden = GOLDEN["slices"]

    assert len(actual) == len(golden), (
        f"slice count mismatch: got {len(actual)} energy slices "
        f"{[c for _, c in actual]}, expected {len(golden)} "
        f"{[s['count'] for s in golden]}"
    )
    for (toten, count), ref in zip(actual, golden):
        tol = max(ABS_FLOOR, REL_TOL * ref["count"])
        assert abs(count - ref["count"]) <= tol, (
            f"resonance-line count for slice toten~{toten:.6g}: "
            f"got {count}, golden {ref['count']} (tol +-{tol:.1f}). "
            "A large drop indicates the get_rescond weight-corruption "
            "regression."
        )


def test_omp_determinism(reslines_runs) -> None:
    """Sorted fort.31415 is identical for OMP_NUM_THREADS=1 and =4."""
    sorted1 = sorted(reslines_runs[1])
    sorted4 = sorted(reslines_runs[4])
    assert sorted1 == sorted4, (
        "OMP_NUM_THREADS=4 produced a different sorted resonance-line set "
        f"than OMP_NUM_THREADS=1 ({len(sorted1)} vs {len(sorted4)} lines); "
        "the parallel solver is not deterministic."
    )


def test_torque_weight_and_residence_ledgers_close(tmp_path) -> None:
    """Independent sums recover the accepted scalar, modes, and box measure."""
    ensure_inputs()
    run_dir = tmp_path / "ledger"
    run_case(find_executable(), run_dir, threads=1)

    torque = float(np.loadtxt(run_dir / "integral_torque.dat"))
    weights = np.loadtxt(run_dir / "potato_resonance_weights.dat", ndmin=2)
    boxes = np.loadtxt(run_dir / "boxcounted_torque.dat", ndmin=2)
    outside = np.loadtxt(
        run_dir / "boxcounted_torque_outside.dat", ndmin=2
    )
    completeness = np.loadtxt(
        run_dir / "potato_completeness.dat", dtype=int, ndmin=2
    )
    mode_integral = np.loadtxt(
        run_dir / "subint_ofH0int_104_vsJperp_adapt.dat", ndmin=2
    )

    # Ledger columns are independently parsed from the executable schema.
    assert weights.shape[1] == 14
    assert np.isclose(weights[:, 12].sum(), torque, rtol=2e-13, atol=1e-10)
    assert np.isclose(
        boxes[:, 0].sum() + outside[:, 0].sum(),
        torque,
        rtol=2e-13,
        atol=1e-10,
    )

    root_disposition = weights[:, 6].astype(int)
    box_disposition = weights[:, 7].astype(int)
    assert set(root_disposition) <= {0, 1}
    assert set(box_disposition) <= {-1, 0}
    assert np.all(weights[root_disposition == 1, 12] == 0.0)
    assert np.all(weights[box_disposition == -1, 12] == 0.0)
    traced_fraction = weights[box_disposition == 0, 13]
    assert np.all(traced_fraction >= -1e-12)
    assert np.all(traced_fraction <= 1.0 + 1e-12)

    # One row per accepted root and no failed box disposition in any energy.
    for row in completeness:
        energy_index = row[0]
        status = row[1]
        root_succeeded = row[13]
        box_attempted, box_succeeded, box_failed, box_skipped = row[16:20]
        energy_rows = weights[weights[:, 0] == energy_index]
        assert status == 0
        assert energy_rows.shape[0] == root_succeeded
        assert box_failed == 0
        assert box_attempted == box_succeeded
        assert box_attempted + box_skipped == root_succeeded

    # Final adaptive-integral columns are ordered by canonical m2.
    final_modes = mode_integral[-1, 2:]
    mode_labels = np.arange(-3, 4)
    weight_modes = np.array(
        [weights[weights[:, 4] == mode, 12].sum() for mode in mode_labels]
    )
    assert np.allclose(weight_modes, final_modes, rtol=2e-13, atol=1e-10)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
