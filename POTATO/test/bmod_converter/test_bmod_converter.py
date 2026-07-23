#!/usr/bin/env python3
"""Self-contained regression for the no-filter Boozer-to-POTATO converter."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator


TOOL = Path(__file__).parents[2] / "tools" / "boozer_npz_to_bmod_n.py"
SPEC = importlib.util.spec_from_file_location("boozer_npz_to_bmod_n", TOOL)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        chartmap = root / "chartmap.nc"
        components = root / "components.npz"
        output = root / "bmod_n.dat"

        s = np.linspace(0.01, 0.81, 9)
        theta = np.linspace(0.0, 2.0*np.pi, 32, endpoint=False)
        zeta = np.array([0.0, np.pi])
        radius = 20.0*np.sqrt(s)
        r = 160.0 + radius[:, None]*np.cos(theta)
        z_plane = radius[:, None]*np.sin(theta)
        toroidal_shift = 0.08*np.sin(theta)
        x = np.empty((2, theta.size, s.size))
        y = np.empty_like(x)
        z = np.empty_like(x)
        for index, angle in enumerate(zeta):
            # The chart is at constant Boozer zeta, but geometric phi differs
            # by a manufactured, poloidally varying toroidal shift.
            phi_geom = angle-toroidal_shift
            x[index] = (r*np.cos(phi_geom)).T
            y[index] = (r*np.sin(phi_geom)).T
            z[index] = z_plane.T
        with h5py.File(chartmap, "w") as handle:
            handle.attrs["zeta_convention"] = "boozer"
            handle["s"] = s
            handle["rho"] = np.sqrt(s)
            handle["theta"] = theta
            handle["zeta"] = zeta
            handle["x"] = x
            handle["y"] = y
            handle["z"] = z

        modes = np.array([-2, 2])
        coefficients = np.column_stack((0.002*s, (0.01+0.003j)*s))
        np.savez(
            components, boozer_s=s, boozer_m=modes,
            boozer_total=coefficients,
        )
        metadata = MODULE.convert(
            chartmap, components, output, component="total", n_tor=-3,
            s_max=0.81, nrad=257, nzet=257, ntheta=256,
            margin_fraction=0.1,
        )
        rad, zet, values = MODULE.read_bmod_n(output)
        if metadata["signed_toroidal_mode"] != -3:
            raise AssertionError("signed toroidal harmonic was not preserved")
        if json.loads(output.with_suffix(".dat.json").read_text()) != metadata:
            raise AssertionError("metadata sidecar differs from returned provenance")
        if values.shape != (257, 257) or not np.all(np.isfinite(values)):
            raise AssertionError("invalid POTATO output grid")
        surface_errors = metadata["gridding_relative_l2_by_surface"]
        if not np.array_equal(surface_errors["s_tor"], s):
            raise AssertionError("surface-error provenance lost the mapped radial grid")
        if (
            len(surface_errors["relative_l2"]) != s.size
            or not np.all(np.isfinite(surface_errors["relative_l2"]))
            or len(surface_errors["absolute_l2_tesla"]) != s.size
            or len(surface_errors["reference_l2_tesla"]) != s.size
        ):
            raise AssertionError("invalid per-surface gridding-error provenance")

        # Check the resolved interior surfaces after the R-Z gridding step.
        sample = RegularGridInterpolator((rad, zet), values)
        theta_test = np.linspace(0.0, 2.0*np.pi, 128, endpoint=False)
        errors = []
        references = []
        for surface in (2, 4, 6):
            points = np.column_stack((
                160.0 + radius[surface]*np.cos(theta_test),
                radius[surface]*np.sin(theta_test),
            ))
            expected = np.sum(
                coefficients[surface, :, None]
                * np.exp(1j*np.outer(modes, theta_test)), axis=0,
            )
            expected *= np.exp(1j*(-3)*0.08*np.sin(theta_test))
            actual = sample(points)
            errors.append(actual-expected)
            references.append(expected)
        relative_l2 = np.linalg.norm(np.concatenate(errors))/np.linalg.norm(
            np.concatenate(references)
        )
        if relative_l2 > 2.0e-2:
            raise AssertionError(f"R-Z reconstruction error is {relative_l2:.3e}")
        if metadata["toroidal_angle_transform"] != (
                "A_RZ=A_B*exp(i*n*(phi_B-phi_geom))"):
            raise AssertionError("missing Boozer-to-geometric toroidal-angle provenance")
        if not np.isclose(metadata["toroidal_shift_radians_max_abs"], 0.08):
            raise AssertionError("manufactured toroidal shift was not recorded")
        if values[0, 0] != 0.0j or values[-1, -1] != 0.0j:
            raise AssertionError("outside-map grid margin must be zero")
    return 0


if __name__ == "__main__":
    sys.exit(main())
