#!/usr/bin/env python3
"""Run the line-drive physics mutation checks in an isolated checkout.

Use inside an allocated scheduler job with bounded FO_JOBS and thread counts.
The source is restored after each mutation, including on failure. Build errors
do not count as a detected physics mutation.
"""
import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fo", default="fo")
    parser.add_argument("--log-dir", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    source = root / "src" / "line_drive.f90"
    original = source.read_text()
    args.log_dir.mkdir(parents=True, exist_ok=True)
    line = "loc(1) = hmag - (qi/c)*sum((lf%A + gauge_fix_eval(gf, th))*xdot)"
    mutations = [
        ("omit_symplectic", line, "loc(1) = hmag", "test_line_drive",
         "collapse |H_line-H_boozer|/|H_boozer|"),
        ("double_count", line, line + " + wboozer*dBL/lf%bmod",
         "test_line_drive", "collapse |H_line-H_boozer|/|H_boozer|"),
        ("omit_radial_covariant_field", "lf%beta*lf%sgB(1) + ", "",
         "test_line_prestudy", "Cartesian curl dBE"),
    ]
    try:
        for name, old, new, test, marker in mutations:
            if original.count(old) != 1:
                raise RuntimeError(f"{name}: source anchor is not unique")
            source.write_text(original.replace(old, new, 1))
            last_test = root / "build" / "Testing" / "Temporary" / "LastTest.log"
            last_test.unlink(missing_ok=True)
            result = subprocess.run([args.fo, "test", test], cwd=root,
                                    text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
            (args.log_dir / f"{name}.log").write_text(result.stdout)
            test_output = last_test.read_text() if last_test.exists() else ""
            (args.log_dir / f"{name}-LastTest.log").write_text(test_output)
            if result.returncode == 0:
                raise RuntimeError(f"{name}: the broken physics passed")
            if test_output.count(" Testing: ") != 1 or f" Testing: {test}\n" not in test_output:
                raise RuntimeError(f"{name}: the expected named oracle did not execute")
            failures = [row for row in test_output.splitlines()
                        if "FAIL" in row and marker in row]
            if not failures:
                raise RuntimeError(f"{name}: failure did not reach its physics oracle")
            print(f"{name}: detected by {failures[0].strip()}", flush=True)
            source.write_text(original)
    finally:
        source.write_text(original)
    for test in ("test_line_drive", "test_line_prestudy"):
        subprocess.run([args.fo, "test", test], cwd=root, check=True)
        output = (root / "build" / "Testing" / "Temporary" / "LastTest.log").read_text()
        (args.log_dir / f"restored-{test}-LastTest.log").write_text(output)
        if output.count(" Testing: ") != 1 or f" Testing: {test}\n" not in output:
            raise RuntimeError(f"restored {test}: expected named test did not execute")
        if "Test Passed." not in output or f"{test}: all checks passed" not in output:
            raise RuntimeError(f"restored {test}: scientific oracle did not pass")


if __name__ == "__main__":
    main()
