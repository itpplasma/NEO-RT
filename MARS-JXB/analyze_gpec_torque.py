#!/usr/bin/env python3
"""Extract and plot native GPEC ideal and kinetic torque evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

CONTROL_TORQUE_PATTERN = re.compile(
    r"toroidal torque\s*=\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)
RATIONAL_SURFACE_PATTERN = re.compile(
    r"psi\s*=\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
    r",\s*q\s*=\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)


def parse_control_torque(path: Path) -> float:
    """Read the scalar torque printed by gpec_control_n3.out."""
    match = CONTROL_TORQUE_PATTERN.search(path.read_text())
    if match is None:
        raise ValueError(f"no toroidal torque found in {path}")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise ValueError(f"non-finite toroidal torque in {path}")
    return value


def parse_dw_profile(path: Path) -> list[tuple[float, ...]]:
    """Read the six numeric columns in gpec_dw_profile_n3.out."""
    rows: list[tuple[float, ...]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 6:
            continue
        try:
            row = tuple(float(field) for field in fields)
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"non-finite profile row in {path}: {line}")
        rows.append(row)

    if len(rows) < 2:
        raise ValueError(f"fewer than two profile rows found in {path}")
    if any(right[0] <= left[0] for left, right in zip(rows, rows[1:])):
        raise ValueError(f"psi_N is not strictly increasing in {path}")
    return rows


def parse_rational_surfaces(path: Path, toroidal_mode: int) -> list[tuple[float, int]]:
    """Read interior q=m/n surfaces reported by ideal DCON."""
    surfaces: list[tuple[float, int]] = []
    for psi_text, q_text in RATIONAL_SURFACE_PATTERN.findall(path.read_text()):
        psi = float(psi_text)
        q = float(q_text)
        poloidal_mode = round(q * toroidal_mode)
        if abs(q - poloidal_mode / toroidal_mode) < 5.0e-4:
            surfaces.append((psi, poloidal_mode))
    return surfaces


def write_profile_csv(rows: list[tuple[float, ...]], path: Path) -> None:
    """Write the native profile without interpolation or smoothing."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "psi_n",
                "dT_dpsi_native_N_m",
                "two_n_delta_W_native",
                "T_cumulative_native_N_m",
                "two_n_delta_W_cumulative_native",
                "dV_dpsi_n_native",
            )
        )
        writer.writerows(rows)


def plot_profiles(
    rows: list[tuple[float, ...]],
    rational_surfaces: list[tuple[float, int]],
    ideal_torque: float,
    kinetic_torque: float,
    output_dir: Path,
) -> None:
    """Create the two diagnostic figures."""
    import matplotlib.pyplot as plt

    psi = [row[0] for row in rows]
    differential = [row[1] for row in rows]
    cumulative = [row[3] for row in rows]

    plt.rcParams.update(
        {
            "axes.grid": True,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "font.size": 10,
            "grid.alpha": 0.25,
            "savefig.bbox": "tight",
        }
    )
    blue = "#0072B2"
    orange = "#D55E00"
    green = "#009E73"
    gray = "#5C5C5C"
    light_gray = "#B8B8B8"

    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    axes[0].plot(
        psi,
        differential,
        color=blue,
        linewidth=1.6,
        label="kinetic GPEC",
    )
    axes[0].axhline(
        0.0,
        color=gray,
        linewidth=1.2,
        linestyle="--",
        label="ideal torque null",
    )
    for index, (surface_psi, _) in enumerate(rational_surfaces):
        axes[0].axvline(
            surface_psi,
            color=light_gray,
            linewidth=0.8,
            linestyle=":",
            label=r"ideal $q=m/3$ surfaces" if index == 0 else None,
        )
    axes[0].set_xlabel(r"native $\psi_N$")
    axes[0].set_ylabel(r"$dT_\phi/d\psi_N$ [N m]")
    axes[0].set_yscale("symlog", linthresh=1.0e2)
    axes[0].set_title("Differential torque, signed log scale")
    axes[0].legend(frameon=False)

    axes[1].plot(
        psi,
        cumulative,
        color=orange,
        linewidth=1.8,
        label="kinetic profile integral",
    )
    axes[1].axhline(
        kinetic_torque,
        color=green,
        linewidth=1.3,
        linestyle=":",
        label="kinetic global scalar",
    )
    axes[1].axhline(
        ideal_torque,
        color=gray,
        linewidth=1.2,
        linestyle="--",
        label="ideal global residual",
    )
    for surface_psi, _ in rational_surfaces:
        axes[1].axvline(
            surface_psi,
            color=light_gray,
            linewidth=0.8,
            linestyle=":",
        )
    axes[1].set_xlabel(r"native $\psi_N$")
    axes[1].set_ylabel(r"$T_\phi(\leq\psi_N)$ [N m]")
    axes[1].set_yscale("symlog", linthresh=1.0e1)
    axes[1].set_title("Cumulative torque, signed log scale")
    axes[1].legend(frameon=False)

    figure.suptitle("TC24 GPEC kinetic torque, native GPEC sign convention")
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(output_dir / f"gpec_kinetic_torque_profile.{suffix}", dpi=220)
    plt.close(figure)

    labels = (
        "ideal global\nresidual",
        "kinetic global\nscalar",
        "kinetic profile\nendpoint",
    )
    values = (ideal_torque, kinetic_torque, cumulative[-1])
    figure, axis = plt.subplots(figsize=(6.5, 4.2))
    bars = axis.bar(
        labels,
        values,
        color=(gray, green, orange),
        edgecolor="black",
        linewidth=0.5,
    )
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_ylabel(r"$T_\phi$ [N m]")
    axis.set_title("TC24 GPEC torque closure, native GPEC sign convention")
    axis.bar_label(bars, fmt="%.5g", padding=3)
    relative_closure = abs(cumulative[-1] - kinetic_torque) / abs(kinetic_torque)
    axis.text(
        0.02,
        0.93,
        f"profile closure: {relative_closure:.3g} relative",
        ha="left",
        va="top",
        transform=axis.transAxes,
    )
    figure.tight_layout()
    for suffix in ("png", "pdf"):
        figure.savefig(output_dir / f"gpec_ideal_kinetic_torque_closure.{suffix}", dpi=220)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    arguments = parser.parse_args()

    ideal_control = arguments.campaign_root / "run-ideal" / "gpec_control_n3.out"
    kinetic_control = (
        arguments.campaign_root / "run-kinetic" / "gpec_control_n3.out"
    )
    kinetic_profile = (
        arguments.campaign_root / "run-kinetic" / "gpec_dw_profile_n3.out"
    )
    ideal_profile = arguments.campaign_root / "run-ideal" / "gpec_dw_profile_n3.out"
    if ideal_profile.exists():
        raise ValueError("ideal case unexpectedly produced a kinetic torque profile")

    rows = parse_dw_profile(kinetic_profile)
    rational_surfaces = parse_rational_surfaces(
        arguments.campaign_root / "run-ideal" / "dcon.log", toroidal_mode=3
    )
    ideal_torque = parse_control_torque(ideal_control)
    kinetic_torque = parse_control_torque(kinetic_control)
    endpoint = rows[-1][3]
    absolute_closure = abs(endpoint - kinetic_torque)
    relative_closure = absolute_closure / max(abs(kinetic_torque), 1.0e-300)

    arguments.output_dir.mkdir(parents=True, exist_ok=False)
    write_profile_csv(rows, arguments.output_dir / "gpec_kinetic_torque_profile.csv")
    metrics = {
        "ideal_control_torque_native_N_m": ideal_torque,
        "ideal_dw_profile_present": False,
        "kinetic_control_torque_native_N_m": kinetic_torque,
        "kinetic_profile_endpoint_native_N_m": endpoint,
        "profile_control_absolute_closure_N_m": absolute_closure,
        "profile_control_relative_closure": relative_closure,
        "profile_rows": len(rows),
        "psi_n_min": rows[0][0],
        "psi_n_max": rows[-1][0],
        "rational_surfaces": [
            {"psi_n": psi, "m": poloidal_mode, "n": 3}
            for psi, poloidal_mode in rational_surfaces
        ],
    }
    (arguments.output_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    )
    plot_profiles(
        rows,
        rational_surfaces,
        ideal_torque,
        kinetic_torque,
        arguments.output_dir,
    )


if __name__ == "__main__":
    main()
