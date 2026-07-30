#!/usr/bin/env python3
"""Plot a full-radius MARS JPLASMA versus curl(BPLASMA) benchmark."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mars_ampere import (  # noqa: E402
    curl_covariant_staggered,
    grid_to_harmonics,
    harmonics_to_grid,
    lower_mars_density,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evaluate_metric(geometry, radial: np.ndarray, chi: np.ndarray):
    """Evaluate the normalized MARS axisymmetric chart at arbitrary radii."""
    source_s = geometry.s[: geometry.ns_plasma]
    radius_modes = CubicSpline(source_s, geometry.rm[: geometry.ns_plasma], axis=0)
    height_modes = CubicSpline(source_s, geometry.zm[: geometry.ns_plasma], axis=0)
    mode = geometry.modes
    phase = np.exp(1j * np.outer(mode, chi))
    phase_derivative = 1j * mode[:, None] * phase
    radius_coefficients = radius_modes(radial)
    height_coefficients = height_modes(radial)
    radius = np.real(radius_coefficients @ phase)
    height = np.real(height_coefficients @ phase)
    radius_s = np.real(radius_modes.derivative()(radial) @ phase)
    height_s = np.real(height_modes.derivative()(radial) @ phase)
    radius_chi = np.real(radius_coefficients @ phase_derivative)
    height_chi = np.real(height_coefficients @ phase_derivative)
    jacobian = radius * (
        radius_s * height_chi - radius_chi * height_s
    )
    return (
        radius,
        radius_s,
        height_s,
        radius_chi,
        height_chi,
        jacobian,
    )


def reconstruct_current(geometry, magnetic, nchi: int):
    """Derive native MARS current harmonics without reading JPLASMA."""
    radial = geometry.s[: geometry.ns_plasma]
    half = 0.5 * (radial[:-1] + radial[1:])
    chi = np.linspace(0.0, 2.0 * np.pi, nchi, endpoint=False)
    full_metric = evaluate_metric(geometry, radial, chi)
    half_metric = evaluate_metric(geometry, half, chi)
    c1_full = magnetic.c1[: len(radial)]
    c1_half = 0.5 * (c1_full[:-1] + c1_full[1:])
    c2_half = magnetic.c2[: len(half)]
    c3_half = magnetic.c3[: len(half)]
    c2_full = CubicSpline(half, c2_half, axis=0)(radial)

    c1_full_grid = harmonics_to_grid(c1_full, magnetic.modes, chi)
    c2_full_grid = harmonics_to_grid(c2_full, magnetic.modes, chi)
    c1_half_grid = harmonics_to_grid(c1_half, magnetic.modes, chi)
    c2_half_grid = harmonics_to_grid(c2_half, magnetic.modes, chi)
    c3_half_grid = harmonics_to_grid(c3_half, magnetic.modes, chi)
    b1_full = np.full_like(c1_full_grid, np.nan + 1j * np.nan)
    b1_full[1:], _, _ = lower_mars_density(
        *(value[1:] for value in full_metric),
        c1_full_grid[1:],
        c2_full_grid[1:],
    )
    _, b2_half, b3_half = lower_mars_density(
        *half_metric, c1_half_grid, c2_half_grid, c3_half_grid
    )
    assert b3_half is not None
    current_grid = curl_covariant_staggered(
        radial,
        chi,
        magnetic.toroidal_mode,
        b1_full,
        b2_half,
        b3_half,
    )
    return radial, tuple(
        grid_to_harmonics(component, magnetic.modes)
        for component in current_grid
    )


def radial_profiles(actual: np.ndarray, predicted: np.ndarray):
    """Return spectral L2 amplitudes and pointwise relative residual."""
    actual_norm = np.linalg.norm(actual, axis=1)
    predicted_norm = np.linalg.norm(predicted, axis=1)
    scale = max(float(np.max(actual_norm)), np.finfo(float).tiny)
    if np.max(actual_norm) == 0.0:
        raise ValueError("MARS current component is identically zero")
    denominator = np.maximum(actual_norm, 1.0e-12 * scale)
    residual = np.linalg.norm(predicted - actual, axis=1) / denominator
    return actual_norm, predicted_norm, residual


def relative_l2(
    actual: np.ndarray, predicted: np.ndarray, selected: np.ndarray
) -> float:
    denominator = np.linalg.norm(actual[selected])
    if denominator == 0.0:
        raise ValueError("selected MARS current reference is identically zero")
    return float(np.linalg.norm(predicted[selected] - actual[selected]) / denominator)


def write_csv(path: Path, profiles: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("component", "rho_pol", "mars_rms", "neort_rms", "relative_residual")
        )
        for profile in profiles:
            for row in zip(
                profile["radial"],
                profile["actual_norm"],
                profile["predicted_norm"],
                profile["residual"],
                strict=True,
            ):
                writer.writerow((profile["component"], *row))


def render(profiles: list[dict[str, object]], output: Path, case_name: str) -> None:
    colors = ("#0072B2", "#E69F00", "#CC79A7")
    fig, (amplitude_axis, residual_axis) = plt.subplots(
        2, 1, figsize=(7.2, 7.0), sharex=True, constrained_layout=True
    )
    for profile, color in zip(profiles, colors, strict=True):
        label = str(profile["component"])
        radial = np.asarray(profile["radial"])
        amplitude_axis.semilogy(
            radial,
            profile["actual_norm"],
            color=color,
            linestyle="-",
            linewidth=1.7,
            label=f"MARS {label}",
        )
        amplitude_axis.semilogy(
            radial,
            profile["predicted_norm"],
            color=color,
            linestyle="--",
            linewidth=1.5,
            label=f"NEO-RT curl(B) {label}",
        )
        residual_axis.semilogy(
            radial,
            profile["residual"],
            color=color,
            linestyle=profile["residual_style"],
            linewidth=1.6,
            label=label,
        )

    amplitude_axis.set_ylabel(
        r"spectral norm of $\sqrt{g}J^i$ (MARS units)"
    )
    amplitude_axis.set_title(
        f"Independent Ampère reconstruction of MARS current — {case_name}"
    )
    amplitude_axis.grid(True, which="both", alpha=0.25)
    amplitude_axis.legend(
        ncol=2, fontsize=8, loc="upper center", frameon=False
    )
    residual_axis.set_xlabel(r"MARS $\rho_{\rm pol}$")
    residual_axis.set_ylabel("relative spectral residual")
    residual_axis.set_xlim(0.0, 1.0)
    residual_axis.grid(True, which="both", alpha=0.25)
    residual_axis.legend(ncol=3, fontsize=8, loc="upper center", frameon=False)
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("geometry", type=Path)
    parser.add_argument("bplasma", type=Path)
    parser.add_argument("jplasma", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--nchi", type=int, default=1024)
    parser.add_argument("--case-name", default="MARS case")
    parser.add_argument(
        "--interior-cells",
        type=int,
        default=10,
        help="cells excluded only from the separately reported interior metric",
    )
    args = parser.parse_args()

    from rmp_torque.mars import read_bplasma, read_rmzm

    geometry = read_rmzm(args.geometry)
    magnetic = read_bplasma(
        args.bplasma, geometry.s, center_half_mesh=False
    )
    current = read_bplasma(
        args.jplasma, geometry.s, center_half_mesh=False
    )
    if (
        magnetic.toroidal_mode != current.toroidal_mode
        or not np.array_equal(magnetic.modes, current.modes)
    ):
        raise ValueError("BPLASMA and JPLASMA harmonic grids differ")
    radial, predicted = reconstruct_current(geometry, magnetic, args.nchi)
    half = 0.5 * (radial[:-1] + radial[1:])
    actual = (
        current.c1[: len(half)],
        current.c2[: len(radial)],
        current.c3[: len(radial)],
    )
    radial_grids = (half, radial, radial)
    residual_styles = ("-", "--", ":")
    profiles: list[dict[str, object]] = []
    diagnostics: dict[str, object] = {
        "case_name": args.case_name,
        "inputs": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in (
                ("geometry", args.geometry),
                ("bplasma", args.bplasma),
                ("jplasma_reference", args.jplasma),
            )
        },
        "nchi": args.nchi,
        "toroidal_mode": magnetic.toroidal_mode,
        "poloidal_mode_min": int(magnetic.modes[0]),
        "poloidal_mode_max": int(magnetic.modes[-1]),
        "components": {},
        "boundary_radial_derivatives": "not determined from half-mesh B2/B3",
    }
    for index, name in enumerate(("J1", "J2", "J3")):
        valid = np.all(np.isfinite(predicted[index]), axis=1)
        actual_valid = actual[index][valid]
        predicted_valid = predicted[index][valid]
        radial_valid = radial_grids[index][valid]
        actual_norm, predicted_norm, residual = radial_profiles(
            actual_valid, predicted_valid
        )
        interior = np.ones(len(radial_valid), dtype=bool)
        cells = args.interior_cells
        if 2 * cells >= len(interior):
            raise ValueError("interior exclusion removes the entire profile")
        interior[:cells] = False
        interior[-cells:] = False
        component_diagnostics = {
            "resolved_relative_l2": relative_l2(
                actual_valid, predicted_valid, np.ones(len(interior), dtype=bool)
            ),
            "interior_relative_l2": relative_l2(
                actual_valid, predicted_valid, interior
            ),
            "interior_rho_min": float(radial_valid[interior][0]),
            "interior_rho_max": float(radial_valid[interior][-1]),
            "resolved_rho_min": float(radial_valid[0]),
            "resolved_rho_max": float(radial_valid[-1]),
        }
        diagnostics["components"][name] = component_diagnostics
        profiles.append(
            {
                "component": name,
                "radial": radial_valid,
                "actual_norm": actual_norm,
                "predicted_norm": predicted_norm,
                "residual": residual,
                "residual_style": residual_styles[index],
            }
        )

    render(profiles, args.output, args.case_name)
    write_csv(args.output.with_suffix(".csv"), profiles)
    args.output.with_suffix(".json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(diagnostics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
