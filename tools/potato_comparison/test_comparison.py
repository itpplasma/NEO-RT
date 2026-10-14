import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analyze import MODES, deposited, main, potato_window, volume, window
from make_cases import write_profiles


class ComparisonPhysics(unittest.TestCase):
    def test_absolute_normalization_is_not_fitted_away(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "inputs").mkdir()
            cases = []
            for energy in (20, 200, 2000):
                cases.append(dict(id=f"e{energy}", T0_eV=energy, s=.09,
                                  kind="ladder", rho_grid=[.12, .2, .3, .4, .48]))
            cases += [dict(cases[0], id="half", kind="half_amplitude", base="e20"),
                      dict(cases[0], id="resolution", kind="energy_resolution", base="e20")]
            (root / "inputs/scan.json").write_text(json.dumps({"cases": cases}))
            edges = np.linspace(0, 1, 101)
            for case in cases:
                torque = -case["T0_eV"] * 1e-6 * (.25 if case["id"] == "half" else 1)
                d = root / "out" / ("p" + case["id"]) / "run"
                d.mkdir(parents=True)
                tb = torque * np.diff(volume(edges)) / 1e-7
                np.savetxt(d / "boxcounted_torque.dat",
                           np.column_stack((tb, edges[:-1]**2, edges[1:]**2)))
                np.savetxt(d / "integral_torque.dat", [sum(tb)])
                if case["kind"] != "energy_resolution":
                    for rho in case["rho_grid"]:
                        for mode in MODES:
                            n = root / "out/neort" / f"{case['id']}_{rho:.2f}_{mode}"
                            n.mkdir(parents=True)
                            np.savetxt(n / "driftorbit_torque.out",
                                       [[0, 2, 0, 20 * torque, 0, 0]])
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(root), 0)
            # A relative missing 1/2 gives a quarter of the torque. It must fail
            # the absolute cold comparison, rather than changing normalization.
            d = root / "out/pe20/run"
            boxes = np.loadtxt(d / "boxcounted_torque.dat")
            boxes[:, 0] /= 4
            np.savetxt(d / "boxcounted_torque.dat", boxes)
            np.savetxt(d / "integral_torque.dat", [sum(boxes[:, 0])])
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(root), 1)
            result = json.loads((root / "out/analysis.json").read_text())
            self.assertAlmostEqual(result["gates"]["e20_cold"]["ratio"], .25)
            self.assertFalse(result["gates"]["e20_cold"]["pass"])

    def test_shared_profiles_and_thermal_energy(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            for energy in (20, 200, 2000):
                write_profiles({"T0_eV": energy, "s": .09}, d)
                p = np.loadtxt(d / "plasma.in", skiprows=3)
                poly = np.loadtxt(d / "profile_poly.in", skiprows=2)
                v = np.loadtxt(d / "profile.in")
                # Independent physical oracle, including the absolute gradient.
                expected = energy * (1.045 - .5 * p[:, 0])
                np.testing.assert_allclose(p[:, 3], expected, rtol=1e-14)
                np.testing.assert_allclose(np.polyval(poly[2], p[:, 0]), expected, rtol=1e-14)
                mass = 2.014 * 1.6726e-24
                np.testing.assert_allclose(.5 * mass * v[:, 2]**2 / 1.602176e-12,
                                           expected, rtol=1e-14)
                self.assertAlmostEqual(np.polyder(poly[2])[-1], -energy / 2)

    def test_volume_against_pappus(self):
        rho = np.array([.12, .3, .6, .8])
        c = 6.2 - np.sqrt(6.2**2 - .62**2)
        r = np.sqrt(rho**2 * c * (12.4 - rho**2 * c))
        np.testing.assert_allclose(volume(rho), (np.pi * r**2) * (2 * np.pi * 6.2), rtol=1e-12)

    def test_constant_torque_window(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            edges = np.linspace(0, 1, 101)
            torque = -3.7e-5
            boxes = np.column_stack((torque * np.diff(volume(edges)) / 1e-7,
                                     edges[:-1]**2, edges[1:]**2))
            np.savetxt(d / "boxcounted_torque.dat", boxes)
            self.assertAlmostEqual(potato_window(d, .3, .06), torque, places=14)
        grid = np.linspace(.1, .8, 71)
        self.assertAlmostEqual(window(grid, np.full_like(grid, torque), .3, .06), torque, places=14)

    def test_deposition_conserves_volume_weighted_torque(self):
        # For N(rho) V'(rho) = constant, normalized symmetric convolution leaves
        # the source exactly unchanged inside a complete window.
        grid = np.linspace(.1, .8, 7001)
        c = 6.2 - np.sqrt(6.2**2 - .62**2)
        dv = 8 * np.pi**2 * 6.2 * c * grid * (6.2 - grid**2 * c)
        torque = -1 / dv
        expected = -.12 / (volume(.36) - volume(.24))
        for kernel in ("tophat", "arcsine"):
            actual = deposited(grid, torque, .3, .06, .04, kernel)
            self.assertAlmostEqual(actual / expected, 1, delta=2e-6)


if __name__ == "__main__":
    unittest.main()
