#!/usr/bin/env python3
"""Compare paired-Boozer and native-staggered MARS j×b profiles."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mars_vector_to_boozer_netcdf import (  # noqa: E402
    load_manifest,
    toroidal_flux_coordinate,
)
from plot_mars_ampere_comparison import reconstruct_current  # noqa: E402


COMPONENTS = ("total", "vacuum", "plasma_response")
LABELS = {
    "total": "total response",
    "vacuum": "applied vacuum",
    "plasma_response": "plasma response",
}
SOURCE_STYLES = {
    "mars_native_staggered": ("#000000", "-"),
    "neo_rt_boozer_collocated": ("#D55E00", "--"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_outrmar_scale(path: Path) -> tuple[float, float, float]:
    """Return MARS aspect ratio, R0 [m], and B0 [T] from OUTRMAR."""
    with path.open(encoding="ascii") as stream:
        next(stream)
        values = [float(value) for value in next(stream).split()]
    if len(values) < 3:
        raise ValueError("OUTRMAR second record needs aspect ratio, R0, and B0")
    aspect, major_radius, magnetic_field = values[:3]
    if min(aspect, major_radius, magnetic_field) <= 0.0:
        raise ValueError("OUTRMAR normalization values must be positive")
    return aspect, major_radius, magnetic_field


def native_half_mesh_torque(
    weighted_current: tuple[np.ndarray, ...],
    c1_full: np.ndarray,
    c2_half: np.ndarray,
) -> np.ndarray:
    """MARS KCTORQ=2 contraction where exported derivatives are determined."""
    k1_half, k2_full, _ = weighted_current
    cell_count = len(c1_full) - 1
    if k1_half.shape[0] != cell_count or c2_half.shape[0] < cell_count:
        raise ValueError("native MARS half-mesh arrays are incompatible")
    if k2_full.shape[0] != len(c1_full):
        raise ValueError("native MARS full-mesh arrays are incompatible")
    result = np.full(cell_count, np.nan)
    for cell in range(1, cell_count - 1):
        result[cell] = 0.5 * np.sum(
            np.real(
                np.conj(k1_half[cell]) * c2_half[cell]
                - 0.5
                * (
                    np.conj(k2_full[cell]) * c1_full[cell]
                    + np.conj(k2_full[cell + 1]) * c1_full[cell + 1]
                )
            )
        )
    return result


def mars_density_in_stor(
    raw: np.ndarray,
    mars_s: np.ndarray,
    stor: np.ndarray,
    torque_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert MARS dT/dsqrt(psi_pol) to dT/ds_tor without changing cells."""
    if len(raw) + 1 != len(mars_s) or len(stor) != len(mars_s):
        raise ValueError("MARS torque grids are incompatible")
    delta_stor = np.diff(stor)
    if np.any(delta_stor <= 0.0):
        raise ValueError("toroidal flux grid must increase")
    density = (
        4.0
        * np.pi**2
        * torque_scale
        * raw
        * np.diff(mars_s)
        / delta_stor
    )
    half_stor = 0.5 * (stor[:-1] + stor[1:])
    return half_stor, np.sqrt(half_stor), density


def cumulative_cells(
    stor: np.ndarray, density: np.ndarray
) -> np.ndarray:
    """Integrate cell-average density while retaining unresolved boundary NaNs."""
    result = np.full_like(density, np.nan)
    valid = np.isfinite(density)
    running = 0.0
    for index in range(len(density)):
        if valid[index]:
            running += density[index] * (stor[index + 1] - stor[index])
            result[index] = running
    return result


def load_native_profiles(
    manifest: Path, outrmar: Path, nchi: int
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, object]]:
    aspect, major_radius, magnetic_field = read_outrmar_scale(outrmar)
    mu0 = 4.0e-7 * np.pi
    torque_scale = major_radius**3 * magnetic_field**2 / (mu0 * aspect**2)
    profiles: dict[str, dict[str, np.ndarray]] = {}
    reference = None
    provenance: dict[str, object] = {
        "manifest": {"path": str(manifest.resolve()), "sha256": sha256(manifest)},
        "outrmar": {"path": str(outrmar.resolve()), "sha256": sha256(outrmar)},
        "mars_aspect_ratio": aspect,
        "mars_r0_m": major_radius,
        "mars_b0_t": magnetic_field,
        "mars_torqfac_N_m": torque_scale,
    }
    for component in COMPONENTS:
        geometry, _, vector, dpsi, field, _ = load_manifest(
            manifest, component, center_half_mesh=False
        )
        mars_s, current = reconstruct_current(geometry, vector, nchi)
        chi = np.linspace(-np.pi, np.pi, nchi, endpoint=False)
        stor, _ = toroidal_flux_coordinate(geometry, dpsi, field, chi)
        if reference is None:
            reference = (mars_s, stor)
            if not np.isclose(geometry.r0_m, major_radius):
                raise ValueError("OUTRMAR and RMZM R0 differ")
            if not np.isclose(geometry.b0_t, magnetic_field):
                raise ValueError("OUTRMAR and RMZM B0 differ")
        elif not (
            np.array_equal(mars_s, reference[0])
            and np.array_equal(stor, reference[1])
        ):
            raise ValueError("MARS components use different radial grids")
        c1 = vector.c1[: len(mars_s)]
        c2 = vector.c2[: len(mars_s) - 1]
        raw = native_half_mesh_torque(current, c1, c2)
        half_stor, rho, density = mars_density_in_stor(
            raw, mars_s, stor, torque_scale
        )
        profiles[component] = {
            "s_tor": half_stor,
            "rho_tor": rho,
            "density": density,
            "cumulative": cumulative_cells(stor, density),
        }
    return profiles, provenance


def load_boozer_profiles(
    directory: Path,
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, object]]:
    profiles = {}
    provenance: dict[str, object] = {}
    for component in COMPONENTS:
        path = directory / f"{component}_jxb_profile.out"
        values = np.loadtxt(path)
        if values.ndim != 2 or values.shape[1] != 3:
            raise ValueError(f"invalid paired-Boozer profile: {path}")
        profiles[component] = {
            "s_tor": values[:, 0],
            "rho_tor": values[:, 1],
            "density": values[:, 2],
            "cumulative": np.concatenate(
                (
                    [0.0],
                    np.cumsum(
                        0.5
                        * (values[:-1, 2] + values[1:, 2])
                        * np.diff(values[:, 0])
                    ),
                )
            ),
        }
        provenance[component] = {
            "path": str(path.resolve()),
            "sha256": sha256(path),
        }
    return profiles, provenance


def sparse_symlog_ticks(axis: plt.Axes, curves: list[np.ndarray]) -> None:
    finite = np.concatenate(
        [curve[np.isfinite(curve)] for curve in curves]
    )
    magnitude = np.abs(finite[np.nonzero(finite)])
    if not len(magnitude):
        return
    linear_threshold = float(np.percentile(magnitude, 15))
    axis.set_yscale("symlog", linthresh=linear_threshold, linscale=1.5)
    first = int(np.floor(np.log10(linear_threshold)))
    last = int(np.ceil(np.log10(np.max(magnitude))))
    exponents = np.unique(np.rint(np.linspace(first, last, 4)).astype(int))
    positive = 10.0**exponents
    axis.set_yticks(np.concatenate((-positive[::-1], positive)))
    axis.tick_params(axis="y", labelsize=8)


def render(
    mars: dict[str, dict[str, np.ndarray]],
    boozer: dict[str, dict[str, np.ndarray]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(
        3, 2, figsize=(11.2, 9.0), sharex="col", constrained_layout=True
    )
    for row, component in enumerate(COMPONENTS):
        density_curves = []
        for source, profiles in (
            ("mars_native_staggered", mars),
            ("neo_rt_boozer_collocated", boozer),
        ):
            color, style = SOURCE_STYLES[source]
            label = {
                "mars_native_staggered": "native-staggered MARS B→J→j×b",
                "neo_rt_boozer_collocated": "paired Boozer NetCDF (collocated)",
            }[source]
            curve = profiles[component]
            axes[row, 0].plot(
                curve["rho_tor"],
                curve["density"],
                color=color,
                linestyle=style,
                linewidth=1.8,
                label=label,
            )
            density_curves.append(curve["density"])
            axes[row, 1].plot(
                curve["rho_tor"],
                curve["cumulative"],
                color=color,
                linestyle=style,
                linewidth=1.8,
                label=label,
            )
        sparse_symlog_ticks(axes[row, 0], density_curves)
        axes[row, 0].set_ylabel(
            rf"{LABELS[component]}" + "\n" + r"$dT_\phi/ds_{\rm tor}$ (N m)"
        )
        axes[row, 0].grid(True, alpha=0.25, which="both")
        axes[row, 1].grid(True, alpha=0.25, which="both")
    axes[0, 0].set_title("signed local torque density")
    axes[0, 1].set_title("cumulative toroidal torque (N m)")
    axes[0, 0].legend(loc="best", fontsize=8)
    axes[0, 1].legend(loc="best", fontsize=8)
    axes[-1, 0].set_xlabel(r"$\rho_{\rm tor}=\sqrt{s_{\rm tor}}$")
    axes[-1, 1].set_xlabel(r"$\rho_{\rm tor}=\sqrt{s_{\rm tor}}$")
    fig.suptitle(
        "ITER TC24 phiI000 electromagnetic torque: "
        "native MARS staggering vs paired Boozer NetCDF"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def write_csv(
    path: Path,
    sources: dict[str, dict[str, dict[str, np.ndarray]]],
) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            [
                "source",
                "component",
                "s_tor",
                "rho_tor",
                "dT_phi_ds_N_m",
                "cumulative_T_phi_N_m",
            ]
        )
        for source, profiles in sources.items():
            for component, curve in profiles.items():
                for values in zip(
                    curve["s_tor"],
                    curve["rho_tor"],
                    curve["density"],
                    curve["cumulative"],
                    strict=True,
                ):
                    writer.writerow((source, component, *values))


def comparison_metrics(
    mars: dict[str, dict[str, np.ndarray]],
    boozer: dict[str, dict[str, np.ndarray]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for component in COMPONENTS:
        native = mars[component]
        paired = boozer[component]
        valid = np.isfinite(native["density"])
        interpolated = np.interp(
            native["s_tor"][valid],
            paired["s_tor"],
            paired["density"],
        )
        denominator = np.linalg.norm(native["density"][valid])
        result[component] = {
            "profile_relative_l2": float(
                np.linalg.norm(interpolated - native["density"][valid])
                / denominator
            ),
            "native_integrated_N_m": float(
                native["cumulative"][np.flatnonzero(valid)[-1]]
            ),
            "boozer_integrated_N_m": float(paired["cumulative"][-1]),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("outrmar", type=Path)
    parser.add_argument("boozer_profile_directory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--nchi", type=int, default=1024)
    args = parser.parse_args()

    mars, mars_provenance = load_native_profiles(
        args.manifest, args.outrmar, args.nchi
    )
    boozer, boozer_provenance = load_boozer_profiles(
        args.boozer_profile_directory
    )
    render(mars, boozer, args.output)
    write_csv(
        args.output.with_suffix(".csv"),
        {
            "mars_native_staggered": mars,
            "neo_rt_boozer_collocated": boozer,
        },
    )
    report = {
        "case_name": "ITER TC24 phiI000",
        "published_radial_coordinate": "rho_tor=sqrt(s_tor)",
        "comparison": comparison_metrics(mars, boozer),
        "mars_provenance": mars_provenance,
        "boozer_provenance": boozer_provenance,
        "interpretation": (
            "The collocated Boozer contract is a smooth-input path. "
            "Its disagreement with native MARS quantifies loss of exported "
            "half-mesh shielding-layer derivatives."
        ),
    }
    args.output.with_suffix(".json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report["comparison"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
