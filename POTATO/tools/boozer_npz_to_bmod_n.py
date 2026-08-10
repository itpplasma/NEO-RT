#!/usr/bin/env python3
"""Map a signed Boozer ``Delta|B|`` spectrum to POTATO's R-Z input grid.

The converter deliberately does not smooth or fit the perturbation.  It uses
the Boozer geometry already stored in a NEO-RT chartmap, reconstructs the
complex single-n amplitude on those surfaces, and linearly interpolates that
amplitude in the poloidal R-Z plane.  Values outside the mapped surface are
zero, and the rectangular grid includes an explicit margin so POTATO does not
silently clamp ordinary orbit excursions to a nonzero boundary value.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import CubicSpline, LinearNDInterpolator, RegularGridInterpolator


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record(payload: bytes) -> bytes:
    size = len(payload)
    return struct.pack("=i", size) + payload + struct.pack("=i", size)


def _text_attribute(value: object) -> str:
    """Return a scalar HDF5 text attribute as a normal Python string."""
    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError("chartmap text attributes must be scalar")
    scalar = array.item()
    if isinstance(scalar, bytes):
        return scalar.decode("utf-8")
    return str(scalar)


def write_bmod_n(path: Path, rad_cm: np.ndarray, zet_cm: np.ndarray,
                 amplitude_tesla: np.ndarray) -> None:
    """Write POTATO's three-record, Fortran-sequential ``bmod_n.dat``."""
    if amplitude_tesla.shape != (rad_cm.size, zet_cm.size):
        raise ValueError("amplitude grid does not match the R-Z axes")
    amplitude_gauss = np.asarray(amplitude_tesla, dtype=np.complex128) * 1.0e4
    with path.open("wb") as stream:
        stream.write(_record(struct.pack("=ii", rad_cm.size, zet_cm.size)))
        stream.write(_record(
            np.asarray(rad_cm, dtype=np.float64).tobytes()
            + np.asarray(zet_cm, dtype=np.float64).tobytes()
        ))
        stream.write(_record(
            amplitude_gauss.real.tobytes(order="F")
            + amplitude_gauss.imag.tobytes(order="F")
        ))


def read_bmod_n(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a POTATO perturbation grid; primarily used for round-trip checks."""
    data = path.read_bytes()
    offset = 0

    def record() -> bytes:
        nonlocal offset
        (size,) = struct.unpack_from("=i", data, offset)
        offset += 4
        payload = data[offset:offset + size]
        offset += size
        (trailer,) = struct.unpack_from("=i", data, offset)
        offset += 4
        if trailer != size:
            raise ValueError("Fortran record markers differ")
        return payload

    nrad, nzet = struct.unpack("=ii", record())
    axes = np.frombuffer(record(), dtype=np.float64)
    values = np.frombuffer(record(), dtype=np.float64)
    if offset != len(data) or axes.size != nrad + nzet or values.size != 2*nrad*nzet:
        raise ValueError("invalid bmod_n.dat record sizes")
    rad = axes[:nrad].copy()
    zet = axes[nrad:].copy()
    count = nrad*nzet
    real = values[:count].reshape((nrad, nzet), order="F")
    imag = values[count:].reshape((nrad, nzet), order="F")
    return rad, zet, (real + 1j*imag)/1.0e4


def convert(chartmap: Path, components: Path, output: Path, *, component: str,
            n_tor: int, s_max: float, nrad: int, nzet: int,
            ntheta: int, margin_fraction: float) -> dict:
    if n_tor == 0:
        raise ValueError("the perturbation toroidal mode must be nonzero")
    if not 0.0 < s_max <= 1.0:
        raise ValueError("s_max must be in (0, 1]")
    if nrad < 8 or nzet < 8:
        raise ValueError("the R-Z grid must have at least 8 points per direction")
    if ntheta < 2:
        raise ValueError("ntheta must be at least two")
    if margin_fraction <= 0.0:
        raise ValueError("margin_fraction must be positive")

    with h5py.File(chartmap, "r") as handle:
        s_geometry = np.asarray(handle["s"], dtype=float)
        theta = np.asarray(handle["theta"], dtype=float)
        zeta = np.asarray(handle["zeta"], dtype=float)
        x = np.asarray(handle["x"], dtype=float)
        y = np.asarray(handle["y"], dtype=float)
        z = np.asarray(handle["z"], dtype=float)
        if "zeta_convention" not in handle.attrs:
            raise ValueError("chartmap lacks the required zeta_convention attribute")
        zeta_convention = _text_attribute(handle.attrs["zeta_convention"])
    if zeta_convention != "boozer":
        raise ValueError(
            f"chartmap zeta_convention must be 'boozer', got {zeta_convention!r}"
        )
    expected = (zeta.size, theta.size, s_geometry.size)
    if x.shape != expected or y.shape != expected or z.shape != expected:
        raise ValueError(f"chartmap x/y/z must have shape {expected}")
    if zeta.size < 1:
        raise ValueError("chartmap zeta axis must not be empty")
    if (np.any(np.diff(s_geometry) <= 0.0)
            or np.any(np.diff(theta) <= 0.0)
            or np.any(np.diff(zeta) <= 0.0)):
        raise ValueError("chartmap s, theta, and zeta axes must increase strictly")

    with np.load(components) as data:
        s_spectrum = np.asarray(data["boozer_s"], dtype=float)
        modes = np.asarray(data["boozer_m"], dtype=int)
        key = f"boozer_{component}"
        if key not in data:
            raise KeyError(f"{components} has no {key}")
        coefficients = np.asarray(data[key], dtype=np.complex128)
    if coefficients.shape != (s_spectrum.size, modes.size):
        raise ValueError("Boozer coefficient matrix has inconsistent dimensions")
    if np.any(np.diff(s_spectrum) <= 0.0):
        raise ValueError("Boozer spectrum surfaces must increase strictly")

    selected = s_geometry <= s_max
    if np.count_nonzero(selected) < 3:
        raise ValueError("s_max leaves fewer than three chartmap surfaces")
    s_map = s_geometry[selected]
    r_chart = np.hypot(x[0][:, selected], y[0][:, selected]).T
    z_chart = z[0][:, selected].T
    theta_map = np.linspace(0.0, 2.0*np.pi, ntheta, endpoint=False)
    theta_closed = np.append(theta, 2.0*np.pi)
    r_map = CubicSpline(
        theta_closed, np.column_stack((r_chart, r_chart[:, 0])), axis=1,
        bc_type="periodic",
    )(theta_map)
    z_map = CubicSpline(
        theta_closed, np.column_stack((z_chart, z_chart[:, 0])), axis=1,
        bc_type="periodic",
    )(theta_map)

    # The spectrum amplitude multiplies exp(i*n*phi_B), whereas POTATO's R-Z
    # amplitude multiplies exp(i*n*phi_geom).  The chartmap is sampled on the
    # first constant-Boozer-zeta plane, whose cylindrical points generally do
    # not all have phi_geom=zeta[0].  Transform the coefficient independently
    # of the signed mode/helicity choice:
    #
    #   A_RZ = A_B * exp(i*n*(phi_B - phi_geom)).
    #
    # Interpolate the unit phasor rather than a wrapped angle, then explicitly
    # restore unit modulus so the coordinate map cannot alter the amplitude.
    phi_geom = np.arctan2(y[0][:, selected], x[0][:, selected]).T
    toroidal_shift_phasor = np.exp(1j*n_tor*(zeta[0] - phi_geom))
    toroidal_shift_map = CubicSpline(
        theta_closed,
        np.column_stack((toroidal_shift_phasor, toroidal_shift_phasor[:, 0])),
        axis=1,
        bc_type="periodic",
    )(theta_map)
    shift_modulus = np.abs(toroidal_shift_map)
    if np.any(shift_modulus <= np.finfo(float).tiny):
        raise ValueError("toroidal-angle phase interpolation produced zero modulus")
    toroidal_shift_map /= shift_modulus

    # Piecewise-linear radial interpolation is exact at the input surfaces and
    # introduces no radial smoothing.  The unresolved axis interval is set to
    # zero, consistent with regular nonaxisymmetric harmonics there.
    coeff_map = np.zeros((s_map.size, modes.size), dtype=np.complex128)
    in_spectrum = s_map >= s_spectrum[0]
    for index in range(modes.size):
        coeff_map[in_spectrum, index] = (
            np.interp(s_map[in_spectrum], s_spectrum, coefficients[:, index].real)
            + 1j*np.interp(s_map[in_spectrum], s_spectrum, coefficients[:, index].imag)
        )
    amplitude_map = (
        coeff_map @ np.exp(1j*np.outer(modes, theta_map))
    ) * toroidal_shift_map

    points = np.column_stack((r_map.ravel(), z_map.ravel()))
    values = amplitude_map.ravel()
    r_span = float(np.ptp(r_map))
    z_span = float(np.ptp(z_map))
    rad_cm = np.linspace(
        float(np.min(r_map) - margin_fraction*r_span),
        float(np.max(r_map) + margin_fraction*r_span), nrad,
    )
    zet_cm = np.linspace(
        float(np.min(z_map) - margin_fraction*z_span),
        float(np.max(z_map) + margin_fraction*z_span), nzet,
    )
    rr, zz = np.meshgrid(rad_cm, zet_cm, indexing="ij")
    interpolator = LinearNDInterpolator(points, values, fill_value=0.0j)
    amplitude_grid = np.asarray(interpolator(rr, zz), dtype=np.complex128)
    if not np.all(np.isfinite(amplitude_grid)):
        raise ValueError("non-finite values produced by R-Z interpolation")

    output.parent.mkdir(parents=True, exist_ok=True)
    write_bmod_n(output, rad_cm, zet_cm, amplitude_grid)
    round_rad, round_zet, round_values = read_bmod_n(output)
    roundtrip_atol = np.finfo(float).eps*max(float(np.max(np.abs(amplitude_grid))), 1.0)
    if not (np.array_equal(round_rad, rad_cm) and np.array_equal(round_zet, zet_cm)
            and np.allclose(round_values, amplitude_grid, rtol=2.0e-16,
                            atol=roundtrip_atol)):
        raise RuntimeError("bmod_n.dat round-trip check failed")

    # Quantify only gridding error at mapped input points.  The no-filter
    # radial and Fourier reconstruction happens before this interpolation.
    grid_interp = RegularGridInterpolator(
        (rad_cm, zet_cm), amplitude_grid, bounds_error=False, fill_value=0.0j
    )
    mapped_back = grid_interp(points).reshape(amplitude_map.shape)
    denominator = max(float(np.linalg.norm(amplitude_map)), np.finfo(float).tiny)
    gridding_relative_l2 = float(np.linalg.norm(mapped_back-amplitude_map)/denominator)
    gridding_absolute_l2_by_surface = np.linalg.norm(
        mapped_back-amplitude_map, axis=1
    )
    surface_reference_l2 = np.linalg.norm(amplitude_map, axis=1)
    reference_floor = np.finfo(float).eps*max(
        float(np.max(surface_reference_l2)), 1.0
    )
    gridding_relative_l2_by_surface = [
        None if reference <= reference_floor else float(error/reference)
        for error, reference in zip(
            gridding_absolute_l2_by_surface, surface_reference_l2, strict=True
        )
    ]
    metadata = {
        "format": "POTATO bmod_n.dat, Fortran sequential, complex amplitude in gauss",
        "inputs": {
            "chartmap": {"path": str(chartmap), "sha256": _sha256(chartmap)},
            "components": {"path": str(components), "sha256": _sha256(components)},
        },
        "output": {"path": str(output), "sha256": _sha256(output)},
        "component": component,
        "signed_toroidal_mode": n_tor,
        "source_fourier_convention": "Delta|B|=Re[A_B*exp(i*n*phi_B)]",
        "output_fourier_convention": "Delta|B|=Re[A_RZ(R,Z)*exp(i*n*phi_geom)]",
        "toroidal_angle_transform": "A_RZ=A_B*exp(i*n*(phi_B-phi_geom))",
        "chartmap_zeta_convention": zeta_convention,
        "chartmap_zeta_slice": float(zeta[0]),
        "toroidal_shift_radians_max_abs": float(np.max(np.abs(np.angle(
            np.exp(1j*(zeta[0] - phi_geom))
        )))),
        "toroidal_phase_correction_radians_max_abs": float(np.max(np.abs(
            np.angle(toroidal_shift_map)
        ))),
        "radial_coordinate": "s_tor",
        "s_map_min": float(s_map[0]),
        "s_map_max": float(s_map[-1]),
        "spectrum_s_min": float(s_spectrum[0]),
        "spectrum_s_max": float(s_spectrum[-1]),
        "rad_cm": [float(rad_cm[0]), float(rad_cm[-1]), nrad],
        "zet_cm": [float(zet_cm[0]), float(zet_cm[-1]), nzet],
        "ntheta_reconstruction": ntheta,
        "chartmap_ntheta": int(theta.size),
        "margin_fraction": margin_fraction,
        "outside_mapped_surface": "zero",
        "radial_interpolation": "piecewise linear; no smoothing or fit",
        "rz_interpolation": "piecewise linear Delaunay; no smoothing or fit",
        "gridding_relative_l2": gridding_relative_l2,
        "gridding_relative_l2_by_surface": {
            "s_tor": s_map.tolist(),
            "absolute_l2_tesla": gridding_absolute_l2_by_surface.tolist(),
            "reference_l2_tesla": surface_reference_l2.tolist(),
            "relative_l2": gridding_relative_l2_by_surface,
            "relative_l2_null_policy": (
                "null when reference_l2_tesla <= eps*max(max_reference,1 T)"
            ),
        },
        "grid_zero_fraction": float(np.count_nonzero(amplitude_grid == 0.0j)/amplitude_grid.size),
        "amplitude_tesla_max": float(np.max(np.abs(amplitude_grid))),
    }
    metadata_path = output.with_suffix(output.suffix + ".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("chartmap", type=Path)
    result.add_argument("components", type=Path)
    result.add_argument("output", type=Path)
    result.add_argument("--component", default="total")
    result.add_argument("--n-tor", type=int, required=True,
                        help="signed native toroidal harmonic used by POTATO")
    result.add_argument("--s-max", type=float, default=0.704,
                        help="largest mapped s_tor surface (default: 0.704)")
    result.add_argument("--nrad", type=int, default=801)
    result.add_argument("--nzet", type=int, default=801)
    result.add_argument("--ntheta", type=int, default=1024)
    result.add_argument("--margin-fraction", type=float, default=0.05)
    return result


def main() -> None:
    args = parser().parse_args()
    metadata = convert(
        args.chartmap, args.components, args.output,
        component=args.component, n_tor=args.n_tor, s_max=args.s_max,
        nrad=args.nrad, nzet=args.nzet, ntheta=args.ntheta,
        margin_fraction=args.margin_fraction,
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
