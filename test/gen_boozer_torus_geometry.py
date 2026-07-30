#!/usr/bin/env python3
"""Generate an analytic circular-torus Boozer chartmap geometry fixture."""

from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset


if len(sys.argv) != 2:
    raise SystemExit("usage: gen_boozer_torus_geometry.py OUTPUT.nc")

output = Path(sys.argv[1])
rho = np.array([0.2, 0.4, 0.6, 0.8])
theta = np.arange(16) * 2.0 * np.pi / 16.0
zeta = np.arange(8) * 2.0 * np.pi / 8.0
major_radius = 3.0
minor_scale = 0.5

rr, tt, pp = np.meshgrid(rho, theta, zeta, indexing="ij")
minor_radius = minor_scale * rr
cylindrical_radius = major_radius + minor_radius * np.cos(tt)
x = cylindrical_radius * np.cos(pp)
y = cylindrical_radius * np.sin(pp)
z = minor_radius * np.sin(tt)

with Dataset(output, "w") as dataset:
    dataset.createDimension("rho", len(rho))
    dataset.createDimension("s", len(rho))
    dataset.createDimension("theta", len(theta))
    dataset.createDimension("zeta", len(zeta))
    dataset.createVariable("rho", "f8", ("rho",))[:] = rho
    dataset.createVariable("s", "f8", ("s",))[:] = rho**2
    dataset.createVariable("theta", "f8", ("theta",))[:] = theta
    dataset.createVariable("zeta", "f8", ("zeta",))[:] = zeta
    for name, values in (("x", x), ("y", y), ("z", z)):
        variable = dataset.createVariable(
            name, "f8", ("zeta", "theta", "rho")
        )
        variable[:] = np.transpose(100.0 * values, (2, 1, 0))
        variable.units = "cm"
    dataset.createVariable("A_phi", "f8", ("s",))[:] = 0.0
    dataset["A_phi"].radial_abscissa = "s"
    dataset.createVariable("B_theta", "f8", ("rho",))[:] = 0.0
    dataset.createVariable("B_phi", "f8", ("rho",))[:] = 1.0
    dataset.createVariable(
        "Bmod", "f8", ("zeta", "theta", "rho")
    )[:] = 1.0
    dataset.createVariable("num_field_periods", "i4").assignValue(1)
    dataset.rho_convention = "rho_tor"
    dataset.zeta_convention = "boozer"
    dataset.rho_lcfs = 1.0
    dataset.boozer_field = 1
    dataset.torflux = 1.0
    dataset.manufactured_geometry = (
        "R=3+0.5*rho*cos(theta), Z=0.5*rho*sin(theta)"
    )
