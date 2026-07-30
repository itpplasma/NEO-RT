"""Independent Ampère reconstruction on the native staggered MARS grid."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline


def spectral_derivative(values: np.ndarray, chi: np.ndarray) -> np.ndarray:
    """Differentiate periodic samples along their final axis."""
    if values.shape[-1] != len(chi):
        raise ValueError("poloidal grid and field shape differ")
    if len(chi) < 3:
        raise ValueError("at least three poloidal samples are required")
    spacing = np.diff(chi)
    if not np.allclose(spacing, spacing[0], rtol=2.0e-13, atol=0.0):
        raise ValueError("poloidal grid must be uniform")
    period = spacing[0] * len(chi)
    if not np.isclose(period, 2.0 * np.pi, rtol=2.0e-13, atol=0.0):
        raise ValueError("poloidal grid must cover one 2*pi period")
    modes = np.fft.fftfreq(len(chi), d=1.0 / len(chi))
    return np.fft.ifft(
        1j * modes * np.fft.fft(values, axis=-1), axis=-1
    )


def lower_mars_density(
    radius: np.ndarray,
    radius_s: np.ndarray,
    height_s: np.ndarray,
    radius_chi: np.ndarray,
    height_chi: np.ndarray,
    jacobian: np.ndarray,
    c1: np.ndarray,
    c2: np.ndarray,
    c3: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Lower MARS's Jacobian-weighted contravariant vector components.

    MARS stores ``C^i = sqrt(g) B^i`` in normalized ``(s,chi,phi)``
    coordinates. This routine returns the covariant components ``B_i``.
    """
    arrays = (
        radius,
        radius_s,
        height_s,
        radius_chi,
        height_chi,
        jacobian,
        c1,
        c2,
    )
    if any(np.shape(value) != np.shape(radius) for value in arrays):
        raise ValueError("metric and vector arrays must have identical shapes")
    if np.any(np.isclose(jacobian, 0.0)):
        raise ValueError("MARS metric has a zero Jacobian")
    g11 = radius_s**2 + height_s**2
    g12 = radius_s * radius_chi + height_s * height_chi
    g22 = radius_chi**2 + height_chi**2
    b1 = (g11 * c1 + g12 * c2) / jacobian
    b2 = (g12 * c1 + g22 * c2) / jacobian
    if c3 is None:
        return b1, b2, None
    if np.shape(c3) != np.shape(radius):
        raise ValueError("third vector component has an incompatible shape")
    b3 = radius**2 * c3 / jacobian
    return b1, b2, b3


def curl_covariant_staggered(
    radial: np.ndarray,
    chi: np.ndarray,
    toroidal_mode: int,
    b1_full: np.ndarray,
    b2_half: np.ndarray,
    b3_half: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``sqrt(g) curl(B)`` on MARS's half/full component meshes.

    ``J1`` is returned on the half mesh. ``J2`` and ``J3`` are returned on
    the full mesh; their two boundary rows are NaN because centered
    half-mesh data do not determine a boundary radial derivative.
    """
    radial = np.asarray(radial, dtype=float)
    half = 0.5 * (radial[:-1] + radial[1:])
    expected_full = (len(radial), len(chi))
    expected_half = (len(half), len(chi))
    if b1_full.shape != expected_full:
        raise ValueError(f"B_1 must have shape {expected_full}")
    if b2_half.shape != expected_half or b3_half.shape != expected_half:
        raise ValueError(f"B_2 and B_3 must have shape {expected_half}")
    if np.any(np.diff(radial) <= 0.0):
        raise ValueError("radial grid must be strictly increasing")

    j1_half = spectral_derivative(b3_half, chi) - (
        1j * toroidal_mode * b2_half
    )
    j2_full = np.full(expected_full, np.nan + 1j * np.nan)
    j3_full = np.full(expected_full, np.nan + 1j * np.nan)
    interior = radial[1:-1]
    db3_ds = CubicSpline(half, b3_half, axis=0).derivative()(interior)
    db2_ds = CubicSpline(half, b2_half, axis=0).derivative()(interior)
    j2_full[1:-1] = (
        1j * toroidal_mode * b1_full[1:-1] - db3_ds
    )
    j3_full[1:-1] = db2_ds - spectral_derivative(
        b1_full[1:-1], chi
    )
    return j1_half, j2_full, j3_full


def harmonics_to_grid(
    coefficients: np.ndarray, modes: np.ndarray, chi: np.ndarray
) -> np.ndarray:
    """Evaluate complex poloidal harmonics on a uniform angle grid."""
    return coefficients @ np.exp(1j * np.outer(modes, chi))


def grid_to_harmonics(
    values: np.ndarray, modes: np.ndarray
) -> np.ndarray:
    """Project uniform periodic samples onto selected integer harmonics."""
    transform = np.fft.fft(values, axis=-1) / values.shape[-1]
    return transform[..., np.mod(modes, values.shape[-1])]
