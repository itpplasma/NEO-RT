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


def manufactured_amplitude_identity() -> None:
    """Check Eq. 10's signed pair against POTATO's one-real-field mode.

    For a physical field ``Re[A exp(i*n*phi)]``, the mathematical Fourier pair
    is ``H_n=A/2`` and ``H_-n=conj(A)/2``.  Thus the direct signed-pair power
    is ``|A|^2/2`` and the Eq. 10 prefactor ``pi**(3/2)/2`` becomes the
    executable single-mode prefactor ``pi**(3/2)/4``.
    """
    amplitude = 0.7 - 0.4j
    n_tor = -3
    phi = np.linspace(-np.pi, np.pi, 257, endpoint=False)
    physical = np.real(amplitude*np.exp(1j*n_tor*phi))
    positive = amplitude/2.0
    negative = np.conj(amplitude)/2.0
    reconstructed = positive*np.exp(1j*n_tor*phi) + \
        negative*np.exp(-1j*n_tor*phi)
    if not np.allclose(physical, reconstructed.real):
        raise AssertionError("signed Fourier pair does not reconstruct real field")

    pair_power = MODULE.pair_power_from_two_sided_coefficients(positive, negative)
    direct_eq10 = 0.5*pair_power
    executable_single_n = 0.5*MODULE.pair_power_from_real_field_amplitude(amplitude)
    if not np.isclose(direct_eq10, executable_single_n):
        raise AssertionError("Eq. 10 and the one-real-field /4 normalization differ")
    if not np.allclose(
            MODULE._to_real_field_single_n(positive, MODULE.TWO_SIDED_COMPLEX),
            amplitude):
        raise AssertionError("two-sided coefficient was not converted to A")
    if not np.allclose(
            MODULE._to_real_field_single_n(amplitude, MODULE.REAL_FIELD_SINGLE_N),
            amplitude):
        raise AssertionError("real-field amplitude was unexpectedly rescaled")
    for ambiguous in (None, "", "unspecified"):
        try:
            MODULE._to_real_field_single_n(amplitude, ambiguous)
        except ValueError:
            continue
        raise AssertionError("ambiguous amplitude convention was accepted")


def main() -> int:
    manufactured_amplitude_identity()
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
            amplitude_convention=MODULE.REAL_FIELD_SINGLE_N,
        )
        rad, zet, values = MODULE.read_bmod_n(output)
        if metadata["signed_toroidal_mode"] != -3:
            raise AssertionError("signed toroidal harmonic was not preserved")
        if metadata["input_amplitude_convention"] != MODULE.REAL_FIELD_SINGLE_N:
            raise AssertionError("input amplitude convention was not recorded")
        if metadata["output_amplitude_convention"] != MODULE.REAL_FIELD_SINGLE_N:
            raise AssertionError("output amplitude convention was not recorded")
        if metadata["conjugate_harmonic"] != "implicit; no explicit -n is stored":
            raise AssertionError("conjugate-harmonic policy was not recorded")
        if "full real-field amplitude" not in metadata["source_fourier_convention"]:
            raise AssertionError("real-field source convention is not descriptive")
        if metadata["pair_accounting"]["implicit_single_n"] != "0.5*abs(A)**2":
            raise AssertionError("implicit-pair accounting was not recorded")
        if json.loads(output.with_suffix(".dat.json").read_text()) != metadata:
            raise AssertionError("metadata sidecar differs from returned provenance")
        if values.shape != (257, 257) or not np.all(np.isfinite(values)):
            raise AssertionError("invalid POTATO output grid")

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
        two_sided_components = root / "components_two_sided.npz"
        two_sided_output = root / "two_sided.dat"
        np.savez(
            two_sided_components, boozer_s=s, boozer_m=modes,
            boozer_total=coefficients/2.0,
        )
        two_sided_metadata = MODULE.convert(
            chartmap, two_sided_components, two_sided_output,
            component="total", n_tor=-3, s_max=0.81, nrad=257, nzet=257,
            ntheta=256, margin_fraction=0.1,
            amplitude_convention=MODULE.TWO_SIDED_COMPLEX,
        )
        _, _, two_sided_values = MODULE.read_bmod_n(two_sided_output)
        if not np.allclose(values, two_sided_values):
            raise AssertionError("equivalent amplitude representations disagree")
        if two_sided_metadata["input_amplitude_convention"] != (
                MODULE.TWO_SIDED_COMPLEX):
            raise AssertionError("two-sided source convention was not recorded")
        if "two-sided" not in two_sided_metadata["source_fourier_convention"]:
            raise AssertionError("two-sided source convention is not descriptive")

        try:
            MODULE.convert(
                chartmap, components, root / "ambiguous.dat", component="total",
                n_tor=-3, s_max=0.81, nrad=9, nzet=9, ntheta=16,
                margin_fraction=0.1,
            )
        except ValueError as error:
            if "ambiguous" not in str(error):
                raise AssertionError("ambiguous-input rejection is not explicit") from error
        else:
            raise AssertionError("missing amplitude provenance was accepted")
        if not np.isclose(metadata["toroidal_shift_radians_max_abs"], 0.08):
            raise AssertionError("manufactured toroidal shift was not recorded")
        if values[0, 0] != 0.0j or values[-1, -1] != 0.0j:
            raise AssertionError("outside-map grid margin must be zero")
    return 0


if __name__ == "__main__":
    sys.exit(main())
