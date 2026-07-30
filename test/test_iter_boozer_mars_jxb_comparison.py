"""Independent behavioral checks for radial MARS torque conversion."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np


PATH = (
    Path(__file__).parents[1]
    / "python"
    / "plot_iter_boozer_mars_jxb_comparison.py"
)
SPEC = spec_from_file_location("plot_iter_boozer_mars_jxb_comparison", PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_native_contraction_matches_explicit_phase_average():
    phases = np.linspace(0.0, 2.0 * np.pi, 2048, endpoint=False)
    phase = np.exp(1j * phases)
    k1 = np.array([[1.2 - 0.7j], [0.5 + 0.9j], [-0.4 + 0.2j]])
    k2 = np.array([[0.3 + 0.1j], [-0.8 + 1.1j], [0.6 - 0.5j], [0.2 + 0.4j]])
    c1 = np.array([[0.7 + 0.2j], [-0.3 + 0.8j], [1.0 - 0.4j], [0.2 - 0.9j]])
    c2 = np.array([[0.5 - 0.6j], [0.9 + 0.3j], [-0.2 + 0.7j]])
    got = MODULE.native_half_mesh_torque(
        (k1, k2, np.zeros_like(k2)), c1, c2
    )
    cell = 1
    k2_center = 0.5 * (
        np.real(k2[cell, 0] * phase)
        + np.real(k2[cell + 1, 0] * phase)
    )
    c1_center = 0.5 * (
        np.real(c1[cell, 0] * phase)
        + np.real(c1[cell + 1, 0] * phase)
    )
    # KCTORQ=2 averages the endpoint products, not the centered fields.
    endpoint_product = 0.5 * (
        np.real(k2[cell, 0] * phase) * np.real(c1[cell, 0] * phase)
        + np.real(k2[cell + 1, 0] * phase)
        * np.real(c1[cell + 1, 0] * phase)
    )
    expected = np.mean(
        np.real(k1[cell, 0] * phase) * np.real(c2[cell, 0] * phase)
        - endpoint_product
    )
    np.testing.assert_allclose(got[cell], expected, rtol=0.0, atol=2.0e-15)
    assert not np.isclose(np.mean(k2_center * c1_center), np.mean(endpoint_product))


def test_radial_coordinate_conversion_preserves_integrated_torque():
    mars_s = np.array([0.1, 0.2, 0.5, 0.9])
    stor = np.array([0.02, 0.08, 0.45, 1.0])
    raw = np.array([2.0, 1.0, 0.25])
    scale = 7.5
    _, _, density = MODULE.mars_density_in_stor(raw, mars_s, stor, scale)
    expected = 4.0 * np.pi**2 * scale * np.sum(raw * np.diff(mars_s))
    np.testing.assert_allclose(
        np.sum(density * np.diff(stor)), expected, rtol=2.0e-15
    )
