#!/usr/bin/env python3
"""Compare a NEO-RT port executable's output to the Fortran golden record.

Usage: compare_golden.py <executable> [rtol]

Runs the executable for every golden case single-threaded and compares the 5
output arrays. Default rtol 1e-3 is the accepted cross-language tolerance:
different correct Adams integrators (DVODE vs CVODE) diverge most on the
near-separatrix trapped resonances, which drives totals to ~1e-3. The dominant
transport and the field/magfie output match far tighter (7 digits / bit-exact).
Run with NEO_RT cases under OMP_NUM_THREADS=1.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import f90nml
import h5py
import numpy as np

GD = Path(__file__).resolve().parents[2] / "test" / "golden_record"
INP = GD / "input"
ARRAYS = [("output", ".out", 1), ("torque", "_torque.out", 1),
          ("integral", "_integral.out", 0),
          ("torque_integral", "_torque_integral.out", 0),
          ("magfie", "_magfie.out", 0)]


def main():
    exe = Path(sys.argv[1]).resolve()
    rtol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
    with h5py.File(GD / "golden.h5", "r") as f:
        cases = sorted(k for k in f.keys() if k.startswith("0p"))
    env = dict(os.environ, OMP_NUM_THREADS="1")
    worst = 0.0
    failed = []
    for case in cases:
        s = float(case.replace("p", "."))
        with tempfile.TemporaryDirectory() as td:
            w = Path(td)
            for fn in ["in_file", "in_file_pert", "plasma.in", "profile.in"]:
                (w / fn).symlink_to(INP / fn)
            nml = f90nml.read(INP / "template.in")
            nml["params"]["s"] = s
            nml.write(w / f"{case}.in", force=True)
            r = subprocess.run([str(exe), case], cwd=w, env=env,
                               capture_output=True, text=True)
            if r.returncode != 0:
                failed.append((case, "run failed: " + r.stderr[-200:]))
                continue
            with h5py.File(GD / "golden.h5", "r") as f:
                g = f[case]
                for name, suf, skip in ARRAYS:
                    a = np.loadtxt(w / f"{case}{suf}", skiprows=skip)
                    b = g[name][:]
                    if a.shape != b.shape:
                        failed.append((case, f"{name} shape {a.shape} vs {b.shape}"))
                        continue
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
                    mr = float(np.nanmax(np.where(np.abs(b) > 1e-290, rel, 0.0)))
                    worst = max(worst, mr)
                    if not np.allclose(a, b, rtol=rtol, atol=1e-12):
                        failed.append((case, f"{name} max rel {mr:.2e} > {rtol:.0e}"))
    print(f"worst rel diff across all cases/arrays: {worst:.2e} (rtol {rtol:.0e})")
    for c, m in failed:
        print(f"  FAIL {c}: {m}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
