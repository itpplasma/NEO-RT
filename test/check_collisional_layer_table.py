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


def converged_ratio(a, b):
    if a + b == 0:
        return 0.0
    n = max(2048, int((a + b)*128))
    r1, r2, r3 = [finite_difference(a, b, n*k) for k in (1, 2, 4)]
    coarse, fine = (4*r2 - r1)/3, (4*r3 - r2)/3
    assert abs(coarse - fine) < 1e-8, (a, b, coarse, fine)
    return fine


def main():
    executable = sys.argv[1]
    for mode, relative, absolute in [('table', 0, 2e-8), ('interior', 0.01, 0),
                                     ('edges', 0.01, 0)]:
        output = subprocess.check_output([executable, mode], text=True)
        worst = 0
        for line in output.splitlines():
            a, b, actual = map(float, line.split())
            expected = converged_ratio(a, b)
            error = abs(actual - expected)
            assert error <= absolute + relative*abs(expected), (a, b, actual, expected)
            if expected != 0:
                worst = max(worst, error/abs(expected))
        print(f"{mode}: {len(output.splitlines())} independent checks, max relative error {worst:.6g}")

    def actual(a, b):
        return float(subprocess.check_output([executable, 'ratio', str(a), str(b)], text=True))

    # Narrow-domain aspect ratios, matching boundary and original blocker.
    for a, b in [(0.1, 0.1), (0, .25), (.25, 0), (.001, .249),
                 (.249, .001), (1e-5, 2e-5), (.125, .125)]:
        expected = converged_ratio(a, b)
        assert abs(actual(a, b)/expected - 1) < 1e-4, (a, b)

    # Thin asymmetric domains that the old minimum-distance clamp mishandled.
    for a, b in [(.001, 1), (1, .001), (1e-9, 2), (2, 1e-9)]:
        expected = converged_ratio(a, b)
        assert abs(actual(a, b)/expected - 1) < 0.01, (a, b)

    # Both asymptotic axes, including one finite boundary and the matching edge.
    for a, b in [(0, 128), (128, 0), (.0625, 64), (64, .0625), (.0625, 128), (128, .0625),
                 (1, 128), (128, 1), (64, 64), (128, 128), (512, 256)]:
        expected = converged_ratio(a, b)
        assert abs(actual(a, b) - expected) < 2e-8, (a, b, actual(a, b), expected)

    # Invalid distances must fail closed, including Fast nonfinite input.
    for a, b, message in [(-1, 1, 'nonnegative'), (1, -1, 'nonnegative'), ('NaN', 1, 'finite'),
                          (1, 'Infinity', 'finite')]:
        result = subprocess.run([executable, 'ratio', str(a), str(b)],
                                capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr, (a, b, result)
    print("Narrow/far limits and explicit invalid-input rejection verified")


if __name__ == "__main__":
    main()
