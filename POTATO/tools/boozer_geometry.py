#!/usr/bin/env python3
"""Sample axisymmetric BC geometry for libneo's Boozer inverse (CGS).

The chart contains smooth geometry only. Perturbation coefficients remain in
their original BC file and are evaluated by NEO-RT at runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.io import netcdf_file


def read_geometry(path):
    surfaces, blocks = [], []
    for line in Path(path).read_text().splitlines()[6:]:
        try:
            row = [float(v.replace("D", "E")) for v in line.split()]
        except ValueError:
            continue
        if len(row) == 6:
            surfaces.append(row)
            blocks.append([])
        elif len(row) == 10 and blocks:
            blocks[-1].append(row)
    surfaces, blocks = np.asarray(surfaces), np.asarray(blocks)
    if surfaces.shape[0] < 4 or blocks.shape[0] != surfaces.shape[0]:
        raise ValueError("expected an eight-column geometry BC file")
    if np.any(blocks[:, :, 1] != 0):
        raise ValueError("geometry must be axisymmetric")
    if np.any(np.diff(surfaces[:, 0]) <= 0):
        raise ValueError("BC flux labels must increase")
    coeff = blocks[:, :, 2:8:2] - 1j*blocks[:, :, 3:8:2]
    return surfaces[:, 0], blocks[0, :, 0], coeff


def write_chart(path, s, theta, geometry):
    """geometry=(R [cm], Z [cm], phi_geom-phi_B [rad])."""
    zeta = np.arange(4)*2*np.pi/4
    with netcdf_file(path, "w", version=2) as nc:
        for name, grid in (("rho", s), ("theta", theta), ("zeta", zeta)):
            nc.createDimension(name, len(grid))
            nc.createVariable(name, "d", (name,))[:] = grid
        nc.zeta_convention = "boozer"
        nc.rho_convention = "psi_tor_norm"
        nc.num_field_periods = np.int32(1)
        for k, name in enumerate(("x", "y", "z")):
            v = nc.createVariable(name, "d", ("zeta", "theta", "rho"))
            v.units = "cm"
            for j, phi in enumerate(zeta):
                if k == 0:
                    data = geometry[:, :, 0]*np.cos(phi + geometry[:, :, 2])
                elif k == 1:
                    data = geometry[:, :, 0]*np.sin(phi + geometry[:, :, 2])
                else:
                    data = geometry[:, :, 1]
                v[j] = data.T


def convert(source, output, ns=4001, ntheta=768):
    if ns < 8 or ntheta < 8:
        raise ValueError("geometry grids need at least eight points")
    native_s, modes, coeff = read_geometry(source)
    s = np.linspace(0.0, 1.0, ns)
    theta = np.arange(ntheta)*2*np.pi/ntheta
    # Match the original BC geometry evaluator, including its radial convention.
    c = CubicSpline(np.sqrt(native_s), coeff, axis=0)(np.sqrt(s))
    geometry = np.real(np.einsum("smk,mt->stk", c, np.exp(1j*modes[:, None]*theta)))
    geometry *= [100.0, 100.0, 2*np.pi]
    write_chart(output, s, theta, geometry)
    metadata = dict(source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
                    ns=ns, ntheta=ntheta, units="cm", radial_coordinate="s=psi_tor_norm",
                    shift="phi_geom-phi_B = 2*pi*BC_v", perturbation="absent",
                    output_sha256=hashlib.sha256(Path(output).read_bytes()).hexdigest())
    Path(str(output) + ".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ns", type=int, default=4001)
    parser.add_argument("--ntheta", type=int, default=768)
    args = parser.parse_args()
    convert(args.source, args.output, args.ns, args.ntheta)
