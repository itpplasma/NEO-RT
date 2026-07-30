#!/usr/bin/env python3
"""Plot a full-radius MARS jxb profile and its NEO-RT integration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FOUR_PI_SQUARED = 4.0 * np.pi**2
MARS_COLOR = "#0072B2"
NEORT_COLOR = "#D55E00"
RAW_COLOR = "#4D4D4D"


def load_profiles(edge_path: Path, torque_path: Path) -> tuple[np.ndarray, np.ndarray]:
    edges = np.loadtxt(edge_path, usecols=0, ndmin=1)
    mars = np.loadtxt(torque_path, usecols=(0, 1), ndmin=2)
    if edges.size != mars.shape[0] + 1:
        raise ValueError("PROFEQ must have one more row than TORQUEJXB")
    if not np.all(np.diff(edges) > 0.0):
        raise ValueError("PROFEQ radial mesh must be strictly increasing")
    return edges, mars


def cumulative_torque(edges: np.ndarray, density: np.ndarray) -> np.ndarray:
    cell_torque = FOUR_PI_SQUARED * density * np.diff(edges)
    return np.concatenate(([0.0], np.cumsum(cell_torque)))


def plot_comparison(
    edges: np.ndarray,
    mars: np.ndarray,
    reconstruction: np.ndarray,
    output: Path,
    title: str,
) -> dict[str, float | int | str]:
    radius, density = mars.T
    if reconstruction.shape != (radius.size, 3):
        raise ValueError("reconstruction must have radius, raw, and processed columns")
    if not np.allclose(reconstruction[:, 0], radius, rtol=0.0, atol=1.0e-4):
        raise ValueError("MARS and reconstructed radial meshes differ")
    raw = reconstruction[:, 1]
    reconstructed = reconstruction[:, 2]
    mars_cumulative = cumulative_torque(edges, density)
    raw_cumulative = cumulative_torque(edges, raw)
    neort_cumulative = cumulative_torque(edges, reconstructed)
    figure, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)

    axes[0].plot(
        radius,
        raw,
        color=RAW_COLOR,
        lw=1.0,
        alpha=0.75,
        label="NEO-RT raw J/B contraction",
    )
    axes[0].plot(radius, density, color=MARS_COLOR, lw=2.0, label="MARS output")
    axes[0].plot(
        radius,
        reconstructed,
        color=NEORT_COLOR,
        lw=1.5,
        ls="--",
        label="NEO-RT + MARS postprocessing",
    )
    axes[0].set_ylabel(r"native $T_{j\times b}$")
    axes[0].legend(frameon=False, loc="best")
    axes[0].axhline(0.0, color="0.5", lw=0.7)

    axes[1].plot(
        edges, mars_cumulative, color=MARS_COLOR, lw=2.0, label="MARS midpoint sum"
    )
    axes[1].plot(
        edges,
        raw_cumulative,
        color=RAW_COLOR,
        lw=1.0,
        alpha=0.75,
        label="NEO-RT raw cumulative",
    )
    axes[1].plot(
        edges,
        neort_cumulative,
        color=NEORT_COLOR,
        lw=1.5,
        ls="--",
        label="NEO-RT integration",
    )
    axes[1].set(xlabel=r"MARS $\rho_\mathrm{pol}=\sqrt{\psi_\mathrm{pol}}$",
                ylabel="cumulative native torque")
    axes[1].legend(frameon=False, loc="best")
    axes[1].axhline(0.0, color="0.5", lw=0.7)
    axes[1].text(
        0.02,
        0.04,
        f"processed endpoint = {neort_cumulative[-1]:.8e}\n"
        f"max profile residual = {np.max(np.abs(reconstructed - density)):.1e}",
        transform=axes[1].transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.margins(x=0.0)
    figure.suptitle(title)
    figure.text(
        0.5,
        0.005,
        "Raw: unsmoothed. Processed: exact MARS stencil. No fit, sign change, "
        "normalization, or radial remapping.",
        ha="center",
        fontsize=8,
    )
    figure.tight_layout(rect=(0.0, 0.025, 1.0, 0.96))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=200)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)
    return {
        "points": int(radius.size),
        "radial_coordinate": "MARS rho_pol=sqrt(psi_pol)",
        "native_total": float(neort_cumulative[-1]),
        "raw_native_total": float(raw_cumulative[-1]),
        "max_profile_residual": float(np.max(np.abs(reconstructed - density))),
        "relative_l2_residual": float(
            np.linalg.norm(reconstructed - density) / np.linalg.norm(density)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("profeq", type=Path)
    parser.add_argument("torquejxb", type=Path)
    parser.add_argument("reconstruction", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--title", default="MARS–NEO-RT electromagnetic torque profile")
    args = parser.parse_args()
    edges, mars = load_profiles(args.profeq, args.torquejxb)
    reconstruction = np.loadtxt(args.reconstruction, ndmin=2)
    metrics = plot_comparison(edges, mars, reconstruction, args.output, args.title)
    args.output.with_suffix(".json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
