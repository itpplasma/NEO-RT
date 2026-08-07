#!/usr/bin/env python3
"""Prepare the circular POTATO torque case with delta-B/B = 1e-3."""

from __future__ import annotations

import argparse
import shutil
import struct
from pathlib import Path

import numpy as np


R0_M = 1.60
A_M = 0.50
B0_T = 2.0
Q0 = 1.5
QA = 4.0
EPSILON = 1.0e-3


def record(payload: bytes) -> bytes:
    size = len(payload)
    return struct.pack("=i", size) + payload + struct.pack("=i", size)


def write_relative_field(path: Path, r_m: np.ndarray, z_m: np.ndarray) -> None:
    rr, zz = np.meshgrid(r_m, z_m, indexing="ij")
    radius = np.sqrt((rr - R0_M) ** 2 + zz**2)
    radius_clipped = np.minimum(radius, A_M)
    safety_factor = Q0 + (QA - Q0) * (radius_clipped / A_M) ** 2
    dpsi_dr = R0_M * B0_T * radius_clipped / safety_factor
    b_toroidal = R0_M * B0_T / rr
    b_poloidal = dpsi_dr / rr
    bmod_gauss = np.sqrt(b_toroidal**2 + b_poloidal**2) * 1.0e4
    real_part = np.asarray(EPSILON * bmod_gauss, dtype=np.float64)
    imaginary_part = np.zeros_like(real_part)
    r_cm = np.asarray(r_m * 100.0, dtype=np.float64)
    z_cm = np.asarray(z_m * 100.0, dtype=np.float64)
    with path.open("wb") as stream:
        stream.write(record(struct.pack("=ii", r_cm.size, z_cm.size)))
        stream.write(record(r_cm.tobytes(order="F") + z_cm.tobytes(order="F")))
        stream.write(
            record(
                real_part.tobytes(order="F")
                + imaginary_part.tobytes(order="F")
            )
        )


def read_grid(eqdsk: Path) -> tuple[np.ndarray, np.ndarray]:
    # The public circular fixture has a fixed 65x65 box. Its exact dimensions
    # are part of the independently checked EQDSK generator contract.
    nr = nz = 65
    r = R0_M - 1.6 * A_M + np.linspace(0.0, 3.2 * A_M, nr)
    z = -1.6 * A_M + np.linspace(0.0, 3.2 * A_M, nz)
    if not eqdsk.is_file():
        raise FileNotFoundError(eqdsk)
    return r, z


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nenergy", type=int, default=24)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for name in ("circ.eqdsk", "convexwall.dat", "field_divB0.inp"):
        shutil.copy2(args.fixture / name, args.output / name)

    r, z = read_grid(args.fixture / "circ.eqdsk")
    write_relative_field(args.output / "bmod_n.dat", r, z)
    (args.output / "profile_poly_torque.in").write_text(
        "% density, dummy, temperature, zero potential; descending powers of s_pol\n"
        "0 0 0 0 0 0 0 0 -2.5e13 5.0e13\n"
        "0 0 0 0 0 0 0 0 0 0\n"
        "0 0 0 0 0 0 0 0 -1.4e3 2.0e3\n"
        "0 0 0 0 0 0 0 0 0 0\n"
    )
    (args.output / "potato.in").write_text(
        f"""&potato_nml
  itest_type = 3
  edge_extension = .false.
  E_alpha = 5d3
  A_alpha = 2d0
  Z_alpha = 1d0
  rho_pol = 0.5d0
  rho_pol_max = 0.8d0
  scalfac_energy = 1d0
  scalfac_efield = 1d0
  Rmax_orbit = 200d0
  ntimstep = 20
  npoicut = 1000
  m_min = -3
  m_max = 3
  n_tor = 2
  nenerg = {args.nenergy}
  thermen_max = 4d0
  nbox = 50
  adaptive_jperp = .true.
  npoi_init = 9
  nlagr_sampling = 3
  eps_sampling = 1.5d-1
  itermax_sampling = 1
  profile_file = 'profile_poly_torque.in'
/
"""
    )


if __name__ == "__main__":
    main()
