"""Native profile values/gradients against analytic and independent SciPy oracles."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline


def main(binary):
    with tempfile.TemporaryDirectory(prefix="potato-profile-test-") as temporary:
        root = Path(temporary)
        x = np.linspace(0, 1, 201)**1.2
        table = np.column_stack((x, 2e13*(1-.3*x),
            50+4000*(1-x)**3+300*np.exp(-((x-.8)/.1)**2), 1+2*x))
        profile = root / "profile.table"

        def write_table(rows):
            with profile.open("w") as stream:
                stream.write(f"POTATO_PROFILE_TABLE_V1\n{len(rows)}\n")
                np.savetxt(stream, rows, fmt="%.16e")

        write_table(table)
        q = np.concatenate((x, .5*(x[1:]+x[:-1])))
        # Reverse queries to exercise the mutable interval cache on every thread.
        q = q[::-1]
        query = root / "query.in"

        def run(threads, coordinates=q):
            with query.open("w") as stream:
                stream.write(f"{len(coordinates)}\n")
                np.savetxt(stream, coordinates, fmt="%.16e")
            env = dict(os.environ, OMP_NUM_THREADS=str(threads))
            return subprocess.run([binary, str(profile), "5000", "-1", "2", "3", "-.2", "-2.2", str(query)],
                env=env, capture_output=True, text=True)

        expected = np.zeros((len(q), 7))
        expected[:, 0] = q
        for column, scale in ((1, 1), (2, 1/(2*5000)), (3, -4.8032e-10/(3*5000*1.6022e-12))):
            oracle = CubicSpline(x, table[:, column]*scale, bc_type="natural")
            expected[:, column] = oracle(q)
            expected[:, column+3] = oracle(q, 1)/-2
        for threads in (1, 4):
            result = run(threads)
            assert result.returncode == 0, result.stderr
            actual = np.loadtxt(result.stdout.splitlines())
            np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-10)
        # Affine density and Phi have exact values and derivatives independently
        # of a spline implementation and include the negative psi span/charge.
        np.testing.assert_allclose(actual[:, 1], 2e13*(1-.3*q), rtol=2e-15)
        np.testing.assert_allclose(actual[:, 4], 3e12, rtol=2e-12)
        # Preserve the old polynomial branch and its unit/scaling factors.
        with profile.open("w") as stream:
            stream.write("legacy polynomial\ncolumns\n")
            for coef in ([0]*8+[-6e12, 2e13], [0]*10,
                         [0]*8+[-500, 1000], [0]*8+[2, 1]):
                stream.write(" ".join(map(str, coef))+"\n")
        result = run(4)
        assert result.returncode == 0, result.stderr
        legacy = np.loadtxt(result.stdout.splitlines())
        np.testing.assert_allclose(legacy[:, 1], 2e13*(1-.3*q), rtol=2e-15)
        np.testing.assert_allclose(legacy[:, 2], (1000-500*q)/10000, rtol=2e-15)
        phi_scale = -4.8032e-10/(3*5000*1.6022e-12)
        np.testing.assert_allclose(legacy[:, 3], (1+2*q)*phi_scale, rtol=2e-15)
        np.testing.assert_allclose(legacy[:, 4], 3e12, rtol=2e-15)
        np.testing.assert_allclose(legacy[:, 5], .025, rtol=2e-15)
        np.testing.assert_allclose(legacy[:, 6], -phi_scale, rtol=2e-15)
        write_table(table)
        assert run(1, np.array([1.01])).returncode != 0, "outside-support query accepted"
        broken = table.copy(); broken[1, 0] = broken[0, 0]
        write_table(broken)
        assert run(1).returncode != 0, "duplicate coordinate accepted"
        broken = table.copy(); broken[-1, 2] = -1
        write_table(broken)
        assert run(1).returncode != 0, "negative temperature accepted"
    print("native tabulated/legacy profile values, gradients, threads and rejected invalid inputs PASS")


if __name__ == "__main__":
    main(str(Path(sys.argv[1]).resolve()))
