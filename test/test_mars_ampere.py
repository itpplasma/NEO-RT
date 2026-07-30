from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np


PATH = Path(__file__).parents[1] / "python" / "mars_ampere.py"
SPEC = spec_from_file_location("mars_ampere", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_lower_mars_density_matches_physical_dot_products() -> None:
    radial = np.array([0.25, 0.7])[:, None]
    chi = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)[None, :]
    major_radius = 5.8
    minor_radius = 1.3
    radius = major_radius + minor_radius * radial * np.cos(chi)
    radius_s = np.broadcast_to(minor_radius * np.cos(chi), radius.shape)
    height_s = np.broadcast_to(minor_radius * np.sin(chi), radius.shape)
    radius_chi = -minor_radius * radial * np.sin(chi)
    height_chi = minor_radius * radial * np.cos(chi)
    jacobian = radius * (
        radius_s * height_chi - radius_chi * height_s
    )
    c1 = np.broadcast_to((0.4 - 0.1j) * (1.0 + 0.2 * np.cos(chi)), radius.shape)
    c2 = np.broadcast_to((-0.3 + 0.5j) * (1.0 + 0.1 * np.sin(chi)), radius.shape)
    c3 = np.broadcast_to((0.2 + 0.6j) * np.ones_like(chi), radius.shape)

    b1, b2, b3 = MODULE.lower_mars_density(
        radius,
        radius_s,
        height_s,
        radius_chi,
        height_chi,
        jacobian,
        c1,
        c2,
        c3,
    )

    # Independent oracle: construct B=(C^i/J)e_i in cylindrical components,
    # then take physical dot products with each coordinate basis vector.
    b_radius = (c1 * radius_s + c2 * radius_chi) / jacobian
    b_height = (c1 * height_s + c2 * height_chi) / jacobian
    b_toroidal = c3 * radius / jacobian
    oracle1 = b_radius * radius_s + b_height * height_s
    oracle2 = b_radius * radius_chi + b_height * height_chi
    oracle3 = b_toroidal * radius
    assert np.allclose(b1, oracle1, rtol=2.0e-15, atol=2.0e-15)
    assert np.allclose(b2, oracle2, rtol=2.0e-15, atol=2.0e-15)
    assert np.allclose(b3, oracle3, rtol=2.0e-15, atol=2.0e-15)


def test_staggered_curl_matches_analytic_covariant_field() -> None:
    radial = np.linspace(0.15, 0.95, 17) ** 1.3
    half = 0.5 * (radial[:-1] + radial[1:])
    chi = np.linspace(0.0, 2.0 * np.pi, 64, endpoint=False)
    toroidal_mode = -3
    phase2 = np.exp(2j * chi)
    phase4 = np.exp(-4j * chi)
    b1 = radial[:, None] ** 2 * phase2 + (0.2j - 0.1) * phase4
    b2 = half[:, None] ** 3 * phase4
    b3 = (0.5 + half[:, None] ** 2) * phase2

    j1, j2, j3 = MODULE.curl_covariant_staggered(
        radial, chi, toroidal_mode, b1, b2, b3
    )
    oracle1 = 2j * b3 - 1j * toroidal_mode * b2
    oracle2 = (
        1j * toroidal_mode * b1[1:-1]
        - 2.0 * radial[1:-1, None] * phase2
    )
    oracle3 = (
        3.0 * radial[1:-1, None] ** 2 * phase4
        - 2j * radial[1:-1, None] ** 2 * phase2
        + 4j * (0.2j - 0.1) * phase4
    )
    assert np.allclose(j1, oracle1, rtol=2.0e-13, atol=2.0e-13)
    assert np.allclose(j2[1:-1], oracle2, rtol=2.0e-12, atol=2.0e-12)
    assert np.allclose(j3[1:-1], oracle3, rtol=2.0e-12, atol=2.0e-12)
    assert np.all(np.isnan(j2[[0, -1]]))
    assert np.all(np.isnan(j3[[0, -1]]))
