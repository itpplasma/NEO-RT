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
            result = subprocess.run([args.fo, "test", test], cwd=root,
                                    text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
            (args.log_dir / f"{name}.log").write_text(result.stdout)
            if result.returncode == 0:
                raise RuntimeError(f"{name}: the broken physics passed")
            failures = [row for row in result.stdout.splitlines()
                        if "FAIL" in row and marker in row]
            if not failures:
                raise RuntimeError(f"{name}: failure did not reach its physics oracle")
            print(f"{name}: detected by {failures[0].strip()}", flush=True)
            source.write_text(original)
    finally:
        source.write_text(original)
    subprocess.run([args.fo, "test", "test_line_(drive|prestudy)"], cwd=root,
                   check=True)


if __name__ == "__main__":
    main()
