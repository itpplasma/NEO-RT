#!/usr/bin/env python3
import hashlib
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
DEFAULT_CASE = Path(
    "/home/ert/proj/potato_benchmarks/rung4_torque_30835/run/potato_zero_mid"
)


def numeric_hash(path, sort_rows):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if sort_rows:
        order = np.lexsort(tuple(data[:, i] for i in range(data.shape[1] - 1, -1, -1)))
        data = data[order]
    payload = "\n".join(
        " ".join(f"{value:.17e}" for value in row) for row in data
    ) + "\n"
    return list(data.shape), hashlib.sha256(payload.encode("ascii")).hexdigest()


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
        if torque != golden["integral_torque"]:
            print(f"integral_torque mismatch: {torque} != {golden['integral_torque']}")
            return 1

        for name, expected in golden["files"].items():
            shape, digest = numeric_hash(work / name, expected.get("sort_rows", False))
            if shape != expected["shape"] or digest != expected["sha256"]:
                print(f"{name} mismatch: shape={shape} sha256={digest}")
                return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
