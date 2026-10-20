#!/usr/bin/env python3
"""Reproduce R_TABLE: -G'' + i*z*G = 1, G'(-a)=0, G(b)=0.

R = Re integral[-a,b] G(z) dz / pi.  SciPy's adaptive fourth-order
collocation solves the real/imaginary system; its piecewise polynomial is
integrated exactly.  No old table entries enter the solve.  Run from any cwd:
  python3 tools/generate_collisional_layer.py          # update production table
  python3 tools/generate_collisional_layer.py --check  # verify reproducibility
The independent CTest oracle uses a finite-difference matrix instead.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

SOURCE = Path(__file__).resolve().parents[1] / "src/collisional_layer.f90"
BEGIN = "    ! BEGIN GENERATED R_TABLE\n"
END = "    ! END GENERATED R_TABLE\n"


def layer_ratio(a, b):
    def rhs(z, y):
        return np.array([y[1], -z*y[2] - 1, y[3], z*y[0]])

    def boundary(left, right):
        return np.array([left[1], left[3], right[0], right[2]])

    z = np.linspace(-a, b, max(100, int((a + b)*16)))
    solution = solve_bvp(rhs, boundary, z, np.zeros((4, len(z))),
                         tol=2e-10, max_nodes=60000)
    if not solution.success:
        raise RuntimeError(f"R({a},{b}): {solution.message}")
    return solution.sol.integrate(-a, b)[0]/np.pi


def generated_table():
    nodes = 2.0**np.arange(-1, 7)
    lines = [BEGIN, "    real(dp), parameter :: R_TABLE(NGRID, NGRID) = reshape([ &\n"]
    for a in nodes:
        row = [f"{layer_ratio(a, b):.12e}_dp" for b in nodes]
        for i in range(0, 6, 3):
            lines.append("        " + ", ".join(row[i:i+3]) + ", &\n")
        tail = ", &" if a != nodes[-1] else "], &"
        lines.append("        " + ", ".join(row[6:]) + tail + f"  ! d_dt/delta = {a:g}\n")
    lines.extend(["        [NGRID, NGRID], order=[2, 1])\n", END])
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = SOURCE.read_text()
    start = source.index(BEGIN)
    end = source.index(END, start) + len(END)
    replacement = generated_table()
    if args.check:
        if source[start:end] != replacement:
            raise SystemExit("R_TABLE differs; regenerate with tools/generate_collisional_layer.py")
        print("R_TABLE reproduces all 64 collocation solves")
    else:
        SOURCE.write_text(source[:start] + replacement + source[end:])


if __name__ == "__main__":
    main()
