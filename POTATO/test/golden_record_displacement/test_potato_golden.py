#!/usr/bin/env python3
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import numpy as np
except ImportError:
    sys.exit(77)


SKIP = 77
DEFAULT_CASE = Path("potato_benchmarks/rung4_torque_30835/run/potato_zero_mid")


# The historical cross-platform integral-torque spread was 6e-8. Keep margin
# for rounding while checking every output value and its row association.
RTOL = 1.0e-6


def compare_rows(data, reference, key_columns, scale_groups=()):
    """Compare canonical rows, preserving coordinates and mode associations."""
    if data.shape != reference.shape:
        return f"shape {data.shape} != {reference.shape}"
    if not np.all(np.isfinite(data)):
        return "nonfinite output"
    order = lambda rows: np.lexsort(tuple(rows[:, i] for i in reversed(key_columns)))
    data = data[order(data)]
    reference = reference[order(reference)]
    scale = np.abs(reference)
    for group in scale_groups:
        # High torque modes approach roundoff. Use the largest torque in the
        # same row, rather than a sum over unrelated rows or physical units.
        scale[:, group] = np.max(scale[:, group], axis=1, keepdims=True)
    error = np.abs(data - reference)
    mismatch = error > RTOL * scale
    if np.any(mismatch):
        row, col = np.argwhere(mismatch)[0]
        return (
            f"row {row} column {col}: {data[row, col]!r} vs "
            f"{reference[row, col]!r} (error {error[row, col]:.2e} > "
            f"tolerance {RTOL * scale[row, col]:.2e})"
        )
    return None


def stage_inputs(golden_path, inputs, dst):
    """Stage a committed public case named by the golden record."""
    src = golden_path.parent / inputs["dir"]
    dst.mkdir()
    for name in inputs["files"]:
        shutil.copy(src / name, dst / name)
    shutil.copy(golden_path.parent / inputs["potato_in"], dst / "potato.in")


def copy_case(src, dst):
    shutil.copytree(src, dst)
    for pattern in (
        "fort.*",
        "*torque*.dat",
        "subint_ofH0int_104_vsJperp_*.dat",
        "potato.log",
    ):
        for path in dst.glob(pattern):
            path.unlink()


def main():
    if len(sys.argv) != 3:
        print("usage: test_potato_golden.py <potato.x> <golden.json>", file=sys.stderr)
        return 2

    exe = Path(sys.argv[1])
    golden_path = Path(sys.argv[2])
    golden = json.loads(golden_path.read_text())
    # A golden record that names committed public inputs always runs. The
    # private AUG record keeps its skip when its case directory is absent.
    inputs = golden.get("inputs")
    case_dir = Path(os.environ.get("POTATO_GOLDEN_CASE", DEFAULT_CASE))
    threads = os.environ.get("POTATO_GOLDEN_THREADS", "16")

    if not exe.is_file():
        print(f"potato.x not found: {exe}", file=sys.stderr)
        return 1 if inputs else SKIP
    if inputs is None and not case_dir.is_dir():
        return SKIP

    with tempfile.TemporaryDirectory(prefix="potato-golden-") as tmp:
        work = Path(tmp) / "case"
        if inputs:
            stage_inputs(golden_path, inputs, work)
        else:
            copy_case(case_dir, work)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = threads
        result = subprocess.run(
            [str(exe)],
            cwd=work,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=int(os.environ.get("POTATO_GOLDEN_TIMEOUT", "1800")),
        )
        if result.returncode != 0:
            print(result.stdout, file=sys.stdout)
            print(result.stderr, file=sys.stderr)
            return result.returncode

        torque = float(np.loadtxt(work / "integral_torque.dat"))
        ref_torque = golden["integral_torque"]
        if not abs(torque - ref_torque) <= RTOL * abs(ref_torque):
            print(
                f"integral_torque mismatch: {torque!r} vs {ref_torque!r} "
                f"(relative {abs(torque / ref_torque - 1.0):.2e} > {RTOL:.0e})"
            )
            return 1

        with np.load(golden_path.parent / golden["rows_file"]) as references:
            for name, expected in golden["files"].items():
                error = compare_rows(
                    np.loadtxt(work / name, ndmin=2), references[name],
                    expected["key_columns"], expected.get("scale_groups", []),
                )
                if error:
                    print(f"{name} mismatch: {error}")
                    return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
