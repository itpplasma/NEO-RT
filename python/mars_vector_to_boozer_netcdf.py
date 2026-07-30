#!/usr/bin/env python3
"""Convert full MARS magnetic vectors to covariant Boozer harmonics.

The output is the vector NetCDF contract consumed by
``neo_rt_boozer_current.x``.  Unlike the historical ``.bc`` perturbation
conversion, no projection to scalar ``delta|B|`` is made.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tomllib
from typing import Any
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline


class ConversionError(ValueError):
    """Input fields or coordinate systems violate the conversion contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def periodic_interp(
    x: np.ndarray, values: np.ndarray, target: np.ndarray
) -> np.ndarray:
    """Periodic cubic interpolation along the final axis."""
    x0 = x[0]
    wrapped = x0 + np.mod(target - x0, 2.0 * np.pi)
    x_closed = np.append(x, x0 + 2.0 * np.pi)
    values_closed = np.concatenate((values, values[..., :1]), axis=-1)
    real = CubicSpline(
        x_closed, values_closed.real, axis=-1, bc_type="periodic"
    )(wrapped)
    imag = CubicSpline(
        x_closed, values_closed.imag, axis=-1, bc_type="periodic"
    )(wrapped)
    return real + 1j * imag


def evaluate_axis_geometry(
    equilibrium: Any, s: np.ndarray, theta: np.ndarray
) -> dict[str, np.ndarray]:
    """Evaluate R, Z, phi_B-phi_geom and their Boozer derivatives."""
    if np.any(equilibrium.n != 0):
        raise ConversionError("axisymmetric Boozer equilibrium required")
    phase = np.exp(1j * np.outer(equilibrium.m, theta))
    dphase = 1j * equilibrium.m[:, None] * phase

    def evaluate(cosine: np.ndarray, sine: np.ndarray):
        coefficients = cosine - 1j * sine
        interpolator = CubicSpline(
            equilibrium.s, coefficients, axis=0, extrapolate=False
        )
        at_s = interpolator(s)
        radial = interpolator.derivative()(s)
        return (
            np.real(at_s @ phase),
            np.real(radial @ phase),
            np.real(at_s @ dphase),
        )

    radius, radius_s, radius_theta = evaluate(
        equilibrium.rmnc, equilibrium.rmns
    )
    height, height_s, height_theta = evaluate(
        equilibrium.zmnc, equilibrium.zmns
    )
    shift, shift_s, shift_theta = evaluate(
        equilibrium.vmnc, equilibrium.vmns
    )
    shift_factor = 2.0 * np.pi / equilibrium.nper
    return {
        "radius": radius,
        "radius_s": radius_s,
        "radius_theta": radius_theta,
        "height": height,
        "height_s": height_s,
        "height_theta": height_theta,
        "shift": shift_factor * shift,
        "shift_s": shift_factor * shift_s,
        "shift_theta": shift_factor * shift_theta,
    }


def toroidal_flux_coordinate(
    geometry: Any, dpsi_ds: np.ndarray, f_profile: np.ndarray, chi: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized toroidal flux and q from the delivered MARS equilibrium."""
    metric = geometry.evaluate(chi)
    ns = geometry.ns_plasma
    if len(dpsi_ds) != ns or len(f_profile) != ns:
        raise ConversionError("XPLASMA equilibrium profiles do not match geometry")
    if np.any(np.isclose(dpsi_ds[1:], 0.0)):
        raise ConversionError("MARS dpsi/ds vanishes away from the axis")
    dphi_dchi = np.empty_like(metric.r)
    dphi_dchi[1:] = (
        f_profile[1:, None]
        * metric.jacobian[1:]
        / (metric.r[1:] ** 2 * dpsi_ds[1:, None])
    )
    dphi_dchi[0] = dphi_dchi[1]
    q = np.mean(dphi_dchi, axis=1)
    toroidal_flux = cumulative_trapezoid(q, metric.s**2, initial=0.0)
    if np.isclose(toroidal_flux[-1], 0.0):
        raise ConversionError("MARS toroidal flux normalization vanishes")
    return toroidal_flux / toroidal_flux[-1], q


def combine_vectors(
    vectors: list[Any], weights: list[complex]
) -> Any:
    from rmp_torque.mars import MarsVectorHarmonics

    if not vectors:
        raise ValueError("at least one MARS vector is required")
    reference = vectors[0]
    for vector in vectors[1:]:
        if (
            vector.s_count != reference.s_count
            or vector.toroidal_mode != reference.toroidal_mode
            or not np.array_equal(vector.modes, reference.modes)
        ):
            raise ConversionError("MARS vectors have incompatible grids or modes")
    components = []
    for name in ("c1", "c2", "c3"):
        components.append(
            sum(
                weight * np.asarray(getattr(vector, name))
                for vector, weight in zip(vectors, weights, strict=True)
            )
        )
    return MarsVectorHarmonics(
        reference.modes,
        reference.s_count,
        reference.toroidal_mode,
        *components,
    )


def physical_mars_vector(
    geometry: Any,
    vector: Any,
    chi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct the physical cylindrical amplitude in tesla."""
    metric = geometry.evaluate(chi)
    c1, c2, c3 = vector.reconstruct(chi, radial_count=geometry.ns_plasma)
    jacobian = metric.jacobian
    b_r = (c1 * metric.dr_ds + c2 * metric.dr_dchi) / jacobian
    b_z = (c1 * metric.dz_ds + c2 * metric.dz_dchi) / jacobian
    b_phi = c3 * metric.r / jacobian
    return tuple(geometry.b0_t * value for value in (b_r, b_phi, b_z))


def map_chi_to_boozer(
    geometry: Any,
    stor: np.ndarray,
    mars_indices: np.ndarray,
    chi: np.ndarray,
    axis: dict[str, np.ndarray],
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    metric = geometry.evaluate(chi)
    mapped = np.empty((len(stor), len(theta)))
    errors = np.empty(len(stor))
    r_axis = metric.r[0, 0] * geometry.r0_m
    z_axis = metric.z[0, 0] * geometry.r0_m

    for output_index, mars_index in enumerate(mars_indices):
        r_mars = metric.r[mars_index] * geometry.r0_m
        z_mars = metric.z[mars_index] * geometry.r0_m
        angle_mars = np.unwrap(np.arctan2(z_mars - z_axis, r_mars - r_axis))
        angle_boozer = np.unwrap(
            np.arctan2(
                axis["height"][output_index] - z_axis,
                axis["radius"][output_index] - r_axis,
            )
        )
        if angle_mars[-1] < angle_mars[0]:
            angle_mars = angle_mars[::-1]
            chi_branch = chi[::-1]
        else:
            chi_branch = chi
        angle_closed = np.append(angle_mars, angle_mars[0] + 2.0 * np.pi)
        chi_closed = np.append(
            chi_branch,
            chi_branch[0]
            + np.sign(chi_branch[-1] - chi_branch[0]) * 2.0 * np.pi,
        )
        aligned = angle_boozer + 2.0 * np.pi * np.round(
            (np.mean(angle_mars) - np.mean(angle_boozer)) / (2.0 * np.pi)
        )
        aligned = angle_closed[0] + np.mod(
            aligned - angle_closed[0], 2.0 * np.pi
        )
        mapped[output_index] = np.interp(aligned, angle_closed, chi_closed)

        radius_mapped = periodic_interp(
            chi, r_mars[None, :], mapped[output_index]
        )[0].real
        height_mapped = periodic_interp(
            chi, z_mars[None, :], mapped[output_index]
        )[0].real
        coordinate_error = np.hypot(
            radius_mapped - axis["radius"][output_index],
            height_mapped - axis["height"][output_index],
        )
        surface_radius = np.hypot(
            axis["radius"][output_index] - r_axis,
            axis["height"][output_index] - z_axis,
        )
        errors[output_index] = np.linalg.norm(coordinate_error) / np.linalg.norm(
            surface_radius
        )
    return mapped, errors


def transform_vector(
    geometry: Any,
    vector: Any,
    equilibrium: Any,
    stor_all: np.ndarray,
    *,
    nchi: int,
    ntheta: int,
    m_max: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Return stor, modes, and covariant Boozer coefficients (3,s,m)."""
    chi = np.linspace(-np.pi, np.pi, nchi, endpoint=False)
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    valid = (stor_all >= equilibrium.s[0]) & (stor_all <= equilibrium.s[-1])
    valid[0] = False
    mars_indices = np.flatnonzero(valid)
    if len(mars_indices) < 3:
        raise ConversionError("MARS and Boozer grids have insufficient overlap")
    stor = stor_all[mars_indices]
    axis = evaluate_axis_geometry(equilibrium, stor, theta)
    mapped_chi, coordinate_error = map_chi_to_boozer(
        geometry, stor, mars_indices, chi, axis, theta
    )

    mars_components = physical_mars_vector(geometry, vector, chi)
    mapped_components = []
    for component in mars_components:
        mapped = np.empty((len(stor), ntheta), dtype=complex)
        for output_index, mars_index in enumerate(mars_indices):
            mapped[output_index] = periodic_interp(
                chi,
                component[mars_index : mars_index + 1],
                mapped_chi[output_index],
            )[0]
        mapped_components.append(mapped)

    # Re-express the coefficient relative to exp(i*n*phi_B).
    phase = np.exp(-1j * vector.toroidal_mode * axis["shift"])
    b_r, b_phi, b_z = [component * phase for component in mapped_components]
    b_covariant = np.empty((3, len(stor), ntheta), dtype=complex)
    b_covariant[0] = (
        b_r * axis["radius_s"]
        + b_z * axis["height_s"]
        - b_phi * axis["radius"] * axis["shift_s"]
    )
    b_covariant[1] = b_phi * axis["radius"]
    b_covariant[2] = (
        b_r * axis["radius_theta"]
        + b_z * axis["height_theta"]
        - b_phi * axis["radius"] * axis["shift_theta"]
    )

    full_modes = np.fft.fftshift(
        np.fft.fftfreq(ntheta, d=1.0 / ntheta)
    ).astype(int)
    selected = (full_modes >= -m_max) & (full_modes <= m_max)
    modes = full_modes[selected]
    full_coefficients = np.fft.fftshift(
        np.fft.fft(b_covariant, axis=2) / ntheta, axes=2
    )
    coefficients = full_coefficients[:, :, selected]
    total_power = np.sum(np.abs(full_coefficients) ** 2, axis=(0, 2))
    retained_power = np.sum(np.abs(coefficients) ** 2, axis=(0, 2))
    discarded = 1.0 - retained_power / np.maximum(
        total_power, np.finfo(float).tiny
    )
    diagnostics = {
        "max_coordinate_relative_l2": float(np.max(coordinate_error)),
        "max_discarded_vector_power_fraction": float(np.max(discarded)),
    }
    return stor, modes, coefficients, diagnostics


def projection_coefficients(
    equilibrium: Any,
    stor: np.ndarray,
    modes: np.ndarray,
    coefficients: np.ndarray,
    ntheta: int,
) -> np.ndarray:
    """Project a covariant perturbation onto the equilibrium unit field."""
    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    axis = evaluate_axis_geometry(equilibrium, stor, theta)
    phase = np.exp(1j * np.outer(modes, theta))
    b_covariant = np.einsum("csm,mt->cst", coefficients, phase)
    iota = CubicSpline(equilibrium.s, equilibrium.iota)(stor)[:, None]

    # e_phi + iota*e_theta is parallel to the Boozer equilibrium field.
    tangent_r = iota * axis["radius_theta"]
    tangent_z = iota * axis["height_theta"]
    tangent_phi = axis["radius"] * (1.0 - iota * axis["shift_theta"])
    tangent_norm = np.sqrt(tangent_r**2 + tangent_z**2 + tangent_phi**2)
    parallel = (b_covariant[1] + iota * b_covariant[2]) / tangent_norm
    return np.fft.fftshift(
        np.fft.fft(parallel, axis=1) / ntheta, axes=1
    )[:, np.isin(
        np.fft.fftshift(np.fft.fftfreq(ntheta, d=1.0 / ntheta)).astype(int),
        modes,
    )]


def write_vector_netcdf(
    path: Path,
    stor: np.ndarray,
    modes: np.ndarray,
    toroidal_mode: int,
    coefficients: np.ndarray,
    provenance: dict[str, object],
) -> None:
    import netCDF4

    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("s", len(stor))
        dataset.createDimension("mode", len(modes))
        dataset.toroidal_mode = toroidal_mode
        dataset.coordinate_order = "s,phi,theta"
        dataset.component_variance = "covariant"
        dataset.radial_coordinate = "s_tor"
        dataset.fourier_convention = "B_hat(s) exp(i*(n*phi+m*theta))"
        dataset.magnetic_component_units = "T m"
        dataset.source_format = "MARS BPLASMA"
        dataset.provenance_json = json.dumps(provenance, sort_keys=True)
        dataset.createVariable("s", "f8", ("s",))[:] = stor
        dataset.createVariable("m", "i4", ("mode",))[:] = modes
        for component_index, name in enumerate(("B_s", "B_phi", "B_theta")):
            real_variable = dataset.createVariable(
                f"{name}_real", "f8", ("mode", "s")
            )
            imag_variable = dataset.createVariable(
                f"{name}_imag", "f8", ("mode", "s")
            )
            for mode_index in range(len(modes)):
                with warnings.catch_warnings():
                    # netCDF4 1.7.x uses NumPy's deprecated internal shape
                    # setter while copying a strided slice.
                    warnings.filterwarnings(
                        "ignore",
                        message="Setting the shape on a NumPy array",
                        category=DeprecationWarning,
                    )
                    real_variable[mode_index, :] = coefficients[
                        component_index, :, mode_index
                    ].real
                    imag_variable[mode_index, :] = coefficients[
                        component_index, :, mode_index
                    ].imag


def load_manifest(
    path: Path, component: str
) -> tuple[Any, Any, Any, np.ndarray, np.ndarray, dict]:
    from rmp_torque.boozer import read_boozer_file
    from rmp_torque.mars import (
        MarsVectorHarmonics,
        read_bplasma,
        read_rmzm,
        read_xplasma,
    )

    manifest = tomllib.loads(path.read_text())
    root = path.parent
    geometry_path = (root / manifest["geometry"]).resolve()
    equilibrium_path = (root / manifest["equilibrium_boozer"]).resolve()
    geometry = read_rmzm(geometry_path)
    equilibrium = read_boozer_file(equilibrium_path)
    vectors = []
    weights = []
    profile_dpsi = None
    profile_f = None
    inputs: dict[str, object] = {
        "manifest": {"path": str(path.resolve()), "sha256": sha256(path)},
        "geometry": {
            "path": str(geometry_path),
            "sha256": sha256(geometry_path),
        },
        "equilibrium": {
            "path": str(equilibrium_path),
            "sha256": sha256(equilibrium_path),
        },
        "component": component,
        "rows": {},
    }
    for name, row in manifest["rows"].items():
        weight = (row["current_at"] / row["run_current_at"]) * np.exp(
            1j * np.deg2rad(row["phase_deg_mars"])
        )
        plasma_path = (root / row["with_plasma"]).resolve()
        vacuum_path = (root / row["vacuum"]).resolve()
        displacement_path = (root / row["displacement"]).resolve()
        plasma = read_bplasma(plasma_path, geometry.s, center_half_mesh=True)
        vacuum = read_bplasma(vacuum_path, geometry.s, center_half_mesh=True)
        if component == "total":
            vector = plasma
        elif component == "vacuum":
            vector = vacuum
        else:
            vector = MarsVectorHarmonics(
                plasma.modes,
                plasma.s_count,
                plasma.toroidal_mode,
                plasma.c1 - vacuum.c1,
                plasma.c2 - vacuum.c2,
                plasma.c3 - vacuum.c3,
            )
        vectors.append(vector)
        weights.append(weight)
        displacement = read_xplasma(
            displacement_path, geometry.s, center_half_mesh=True
        )
        if profile_dpsi is None:
            profile_dpsi = displacement.dpsi_ds
            profile_f = displacement.f
        elif not (
            np.allclose(profile_dpsi, displacement.dpsi_ds)
            and np.allclose(profile_f, displacement.f)
        ):
            raise ConversionError("row XPLASMA equilibrium profiles differ")
        inputs["rows"][name] = {
            "weight_real": float(weight.real),
            "weight_imag": float(weight.imag),
            "with_plasma_sha256": sha256(plasma_path),
            "vacuum_sha256": sha256(vacuum_path),
            "displacement_sha256": sha256(displacement_path),
        }
    assert profile_dpsi is not None and profile_f is not None
    return (
        geometry,
        equilibrium,
        combine_vectors(vectors, weights),
        np.asarray(profile_dpsi),
        np.asarray(profile_f),
        inputs,
    )


def compare_scalar_reference(
    reference_path: Path,
    key: str,
    stor: np.ndarray,
    modes: np.ndarray,
    projection: np.ndarray,
) -> dict[str, float]:
    with np.load(reference_path) as reference:
        reference_s = reference["boozer_s"]
        reference_m = reference["boozer_m"]
        reference_values = reference[key]
    if not np.array_equal(reference_m, modes):
        raise ValueError("scalar reference mode grid differs")
    if not np.allclose(reference_s, stor, rtol=0.0, atol=2.0e-12):
        raise ValueError("scalar reference radial grid differs")
    residual = projection - reference_values
    return {
        "scalar_projection_relative_l2": float(
            np.linalg.norm(residual) / np.linalg.norm(reference_values)
        ),
        "scalar_projection_max_abs_t": float(np.max(np.abs(residual))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--component",
        choices=("total", "vacuum", "plasma_response"),
        default="total",
    )
    parser.add_argument("--m-max", type=int)
    parser.add_argument("--nchi", type=int)
    parser.add_argument("--ntheta", type=int)
    parser.add_argument("--scalar-reference", type=Path)
    parser.add_argument("--scalar-key", default="boozer_eulerian")
    args = parser.parse_args()

    manifest = tomllib.loads(args.manifest.read_text())
    m_max = args.m_max if args.m_max is not None else int(manifest["m_max"])
    nchi = args.nchi if args.nchi is not None else int(manifest["nchi"])
    ntheta = args.ntheta if args.ntheta is not None else int(manifest["ntheta"])
    geometry, equilibrium, vector, dpsi_ds, f_profile, provenance = load_manifest(
        args.manifest, args.component
    )
    chi = np.linspace(-np.pi, np.pi, nchi, endpoint=False)
    stor_all, _ = toroidal_flux_coordinate(geometry, dpsi_ds, f_profile, chi)
    stor, modes, coefficients, diagnostics = transform_vector(
        geometry,
        vector,
        equilibrium,
        stor_all,
        nchi=nchi,
        ntheta=ntheta,
        m_max=m_max,
    )
    projection = projection_coefficients(
        equilibrium, stor, modes, coefficients, ntheta
    )
    if args.scalar_reference is not None:
        diagnostics.update(
            compare_scalar_reference(
                args.scalar_reference,
                args.scalar_key,
                stor,
                modes,
                projection,
            )
        )
        provenance["scalar_reference"] = {
            "path": str(args.scalar_reference.resolve()),
            "sha256": sha256(args.scalar_reference),
            "key": args.scalar_key,
        }
    diagnostics.update(
        {
            "surfaces": len(stor),
            "modes": len(modes),
            "toroidal_mode": vector.toroidal_mode,
        }
    )
    provenance["diagnostics"] = diagnostics
    write_vector_netcdf(
        args.output,
        stor,
        modes,
        vector.toroidal_mode,
        coefficients,
        provenance,
    )
    args.output.with_suffix(".json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(diagnostics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
