"""Generate a manufactured covariant Boozer-vector NetCDF fixture."""

from pathlib import Path
import sys
import warnings

import netCDF4
import numpy as np


if len(sys.argv) != 2:
    raise SystemExit("usage: gen_boozer_vector_fixture.py OUTPUT.nc")

output = Path(sys.argv[1])
# netCDF4 1.7.x reaches NumPy's deprecated internal shape setter on assignment;
# this is upstream adapter noise, not a fixture-schema warning.
warnings.filterwarnings(
    "ignore", message="Setting the shape on a NumPy array", category=DeprecationWarning
)
s = np.array([0.1, 0.3, 0.7, 1.0])
modes = np.array([-2, 3], dtype=np.int32)
n = -1

# Vector potentials A_phi=p(s), A_theta=q(s).  The stored field is curl(A)
# in the identity chart, with exp(i*(n*phi+m*theta)) angular dependence.
p = np.column_stack(
    (
        (1.0 + 0.5j) + (2.0 - 0.25j) * s + (0.3 + 0.2j) * s**2,
        (-0.4 + 0.7j) + (0.8 + 0.1j) * s + (-0.2 + 0.15j) * s**2,
    )
)
q = np.column_stack(
    (
        (-0.2 + 0.3j) + (0.4 + 0.6j) * s + (-0.1 + 0.05j) * s**2,
        (0.9 - 0.2j) + (-0.3 + 0.4j) * s + (0.25 - 0.1j) * s**2,
    )
)
p_prime = np.column_stack(
    (
        (2.0 - 0.25j) + 2.0 * (0.3 + 0.2j) * s,
        (0.8 + 0.1j) + 2.0 * (-0.2 + 0.15j) * s,
    )
)
q_prime = np.column_stack(
    (
        (0.4 + 0.6j) + 2.0 * (-0.1 + 0.05j) * s,
        (-0.3 + 0.4j) + 2.0 * (0.25 - 0.1j) * s,
    )
)
b_s = 1j * (n * q - p * modes[None, :])
b_phi = -q_prime
b_theta = p_prime

with netCDF4.Dataset(output, "w") as dataset:
    dataset.createDimension("s", len(s))
    dataset.createDimension("mode", len(modes))
    dataset.toroidal_mode = n
    dataset.coordinate_order = "s,phi,theta"
    dataset.component_variance = "covariant"
    dataset.radial_coordinate = "s_tor"
    dataset.fourier_convention = "B_hat(s) exp(i*(n*phi+m*theta))"
    dataset.magnetic_component_units = "T m"
    dataset.createVariable("s", "f8", ("s",))[:] = s
    dataset.createVariable("m", "i4", ("mode",))[:] = modes
    for name, values in (
        ("B_s", b_s),
        ("B_phi", b_phi),
        ("B_theta", b_theta),
    ):
        # Python/NetCDF order is reversed relative to the Fortran array.
        real_variable = dataset.createVariable(f"{name}_real", "f8", ("mode", "s"))
        imag_variable = dataset.createVariable(f"{name}_imag", "f8", ("mode", "s"))
        for mode_index in range(len(modes)):
            real_variable[mode_index, :] = values[:, mode_index].real
            imag_variable[mode_index, :] = values[:, mode_index].imag
