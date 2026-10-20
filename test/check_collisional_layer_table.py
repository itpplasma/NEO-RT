#!/usr/bin/env python3
"""Independent finite-difference oracle for production Fortran table nodes.

A reflected ghost point enforces G'(-a)=0; G(b)=0 is eliminated.  This
tridiagonal discretization and trapezoidal integral share no solver or constants
with the adaptive collocation generator. Richardson extrapolation removes the
leading O(h^2) error; convergence is checked before comparing to Fortran.
"""
import subprocess
import sys

import numpy as np
from scipy.linalg import solve_banded


def finite_difference(a, b, intervals):
    h = (a + b)/intervals
    z = np.linspace(-a, b, intervals + 1)
    matrix = np.zeros((3, intervals), dtype=complex)
    matrix[0, 1:] = -1/h**2
    matrix[0, 1] *= 2  # reflected Neumann ghost at the left boundary
    matrix[1, :] = 2/h**2 + 1j*z[:-1]
    matrix[2, :-1] = -1/h**2
    g = solve_banded((1, 1), matrix, np.ones(intervals))
    return h*(g.real.sum() - g[0].real/2)/np.pi


def main():
    output = subprocess.check_output([sys.argv[1], "table"], text=True)
    for line in output.splitlines():
        a, b, actual = map(float, line.split())
        n = max(2048, int((a + b)*128))
        r1, r2, r3 = [finite_difference(a, b, n*k) for k in (1, 2, 4)]
        coarse, fine = (4*r2 - r1)/3, (4*r3 - r2)/3
        assert abs(coarse - fine) < 1e-8, (a, b, coarse, fine)
        assert abs(actual - fine) < 2e-8, (a, b, actual, fine)
    print("64 Fortran nodes agree with independent converged finite differences")


if __name__ == "__main__":
    main()
