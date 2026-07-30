from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np


PATH = (
    Path(__file__).parents[1]
    / "python"
    / "mars_vector_to_boozer_netcdf.py"
)
SPEC = spec_from_file_location("mars_vector_to_boozer_netcdf", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def circular_equilibrium() -> SimpleNamespace:
    radial = np.array([0.2, 0.8])
    major_radius = 6.0
    minor_radius = 1.5
    # R=R0+a*s*cos(theta), Z=a*s*sin(theta), zero toroidal shift.
    rmnc = np.column_stack(
        (np.full_like(radial, major_radius), minor_radius * radial)
    )
    zeros = np.zeros_like(rmnc)
    zmns = np.column_stack((np.zeros_like(radial), minor_radius * radial))
    return SimpleNamespace(
        s=radial,
        iota=np.full_like(radial, 0.4),
        m=np.array([0, 1]),
        n=np.array([0, 0]),
        nper=1,
        rmnc=rmnc,
        rmns=zeros,
        zmnc=zeros,
        zmns=zmns,
        vmnc=zeros,
        vmns=zeros,
    )


def test_axis_basis_matches_analytic_circular_chart() -> None:
    equilibrium = circular_equilibrium()
    radial = equilibrium.s
    theta = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    geometry = MODULE.evaluate_axis_geometry(equilibrium, radial, theta)
    s_grid = radial[:, None]
    a = 1.5

    assert np.allclose(
        geometry["radius"], 6.0 + a * s_grid * np.cos(theta), atol=2.0e-15
    )
    assert np.allclose(geometry["radius_s"], a * np.cos(theta), atol=2.0e-15)
    assert np.allclose(
        geometry["radius_theta"], -a * s_grid * np.sin(theta), atol=2.0e-15
    )
    assert np.allclose(
        geometry["height"], a * s_grid * np.sin(theta), atol=2.0e-15
    )
    assert np.allclose(geometry["height_s"], a * np.sin(theta), atol=2.0e-15)
    assert np.allclose(
        geometry["height_theta"], a * s_grid * np.cos(theta), atol=2.0e-15
    )
    assert np.max(np.abs(geometry["shift"])) == 0.0


def test_covariant_projection_matches_physical_tangent_oracle() -> None:
    equilibrium = circular_equilibrium()
    radial = equilibrium.s
    modes = np.arange(-4, 5)
    coefficients = np.zeros((3, len(radial), len(modes)), dtype=complex)
    amplitude = 0.75 - 0.2j
    coefficients[1, :, np.flatnonzero(modes == 0)[0]] = amplitude
    ntheta = 64

    projected = MODULE.projection_coefficients(
        equilibrium, radial, modes, coefficients, ntheta
    )
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    radius = 6.0 + 1.5 * radial[:, None] * np.cos(theta)
    tangent_norm = np.sqrt(radius**2 + (0.4 * 1.5 * radial[:, None]) ** 2)
    physical_oracle = amplitude / tangent_norm
    full_modes = np.fft.fftshift(
        np.fft.fftfreq(ntheta, d=1.0 / ntheta)
    ).astype(int)
    oracle = np.fft.fftshift(
        np.fft.fft(physical_oracle, axis=1) / ntheta, axes=1
    )[:, np.isin(full_modes, modes)]

    assert np.allclose(projected, oracle, rtol=0.0, atol=2.0e-15)


def test_periodic_interpolation_preserves_fourier_mode() -> None:
    x = np.linspace(-np.pi, np.pi, 32, endpoint=False)
    target = np.linspace(-2.9, 3.0, 41)
    values = np.exp(3j * x)[None, :]
    interpolated = MODULE.periodic_interp(x, values, target)[0]
    # The independent analytic oracle includes the cubic interpolation error
    # for 32 samples of the m=3 mode (max 3.42e-4 on this target grid).
    assert np.allclose(interpolated, np.exp(3j * target), rtol=0.0, atol=3.5e-4)
