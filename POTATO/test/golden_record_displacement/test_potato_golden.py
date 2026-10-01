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


# The outputs agree across platforms only to rounding (macOS arm64 gives a
# relative integral-torque difference of 6e-8 against the Linux baseline), so
# compare with a relative tolerance rather than bitwise. 1e-6 leaves a factor
# of about 15 above that spread and still catches any physics change.
RTOL = 1.0e-6


def column_stats(path):
    """Shape and per-column sum and sum of |x|; independent of row order."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return list(data.shape), data.sum(axis=0), np.abs(data).sum(axis=0)


def column_scale(abs_sum, groups):
    """Tolerance scale per column: its own sum of |x|, or the summed sum of |x|
    of its group for columns that share a physical unit (torque per mode)."""
    scale = np.array(abs_sum, dtype=float)
    for group in groups:
        scale[group] = scale[group].sum()
    return scale


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
    case_dir = Path(os.environ.get("POTATO_GOLDEN_CASE", DEFAULT_CASE))
    threads = os.environ.get("POTATO_GOLDEN_THREADS", "16")

    if not exe.is_file() or not case_dir.is_dir():
        return SKIP

    golden = json.loads(golden_path.read_text())
    with tempfile.TemporaryDirectory(prefix="potato-golden-") as tmp:
        work = Path(tmp) / "case"
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
            timeout=int(os.environ.get("POTATO_GOLDEN_TIMEOUT", "180")),
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

        for name, expected in golden["files"].items():
            shape, col_sum, abs_sum = column_stats(work / name)
            if shape != expected["shape"]:
                print(f"{name} shape {shape} != {expected['shape']}")
                return 1
            scale = column_scale(expected["abs_sum"], expected.get("scale_groups", []))
            dev = np.abs(col_sum - np.array(expected["sum"])) / np.where(
                scale > 0.0, scale, 1.0
            )
            if not np.all(dev <= RTOL):
                col = int(np.argmax(dev))
                print(
                    f"{name} column {col} sum {col_sum[col]!r} vs "
                    f"{expected['sum'][col]!r} (relative {dev[col]:.2e} > {RTOL:.0e})"
                )
                return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
