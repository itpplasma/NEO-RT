#!/usr/bin/env python3
"""Plot S3 CSV evidence from diagnose_line_drive.x without changing its units.

The Okabe-Ito blue/orange palette is safe for common color-vision deficiencies.
Solid line versus dashed line and circle markers also distinguish the routes
in grayscale. Sources remain beside the PNG/PDF figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_table(path: Path) -> np.ndarray:
    table = np.genfromtxt(path, names=True, delimiter=",", dtype=None,
                          encoding="utf-8")
    table = np.atleast_1d(table)
    for name in table.dtype.names:
        if table[name].dtype.kind in "if":
            if not np.all(np.isfinite(table[name])):
                raise ValueError(f"{path.name}: non-finite column {name}")
    return table


def save_figure(fig, directory: Path, name: str) -> None:
    for suffix in ("png", "pdf"):
        fig.savefig(directory / f"{name}.{suffix}", dpi=180,
                    bbox_inches="tight")
    plt.close(fig)


def plot_integrands(directory: Path, colors: tuple[str, str], suffix: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 5.8), sharex=True)
    for row, (orbit, harmonic) in enumerate((("passing", -4), ("trapped", 0))):
        table = read_table(directory / f"integrand_{orbit}.csv")
        t = table["time_over_period"]
        if t[0] != 0.0 or abs(t[-1] - 1.0) > 1e-12:
            raise ValueError(f"{orbit}: incomplete sampled orbit")
        if not np.all(np.diff(t) > 0):
            raise ValueError(f"{orbit}: non-monotone sampling times")
        for col, (part, label) in enumerate((("re", "Re"), ("im", "Im"))):
            ax = axes[row, col]
            ax.plot(t, table[f"line_{part}"], color=colors[0],
                    linewidth=1.9, label="Line integral")
            ax.plot(t, table[f"boozer_{part}"], color=colors[1],
                    linewidth=1.5, linestyle="--", label="Boozer scalar")
            ax.set_title(f"{orbit.capitalize()}, $m_b={harmonic}$: {label}")
            ax.set_ylabel(r"Fourier integrand / $E$ [1]")
            ax.grid(alpha=0.2)
            ax.set_xlim(0, 1)
            if row == 1:
                ax.set_xlabel(r"Time / bounce or transit period $t/T$ [1]")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95),
               ncol=2, frameon=False)
    fig.suptitle("Ideal-helical perturbation on examples/base, on resonance",
                 fontsize=11, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save_figure(fig, directory, f"integrands{suffix}")


def plot_harmonics(directory: Path, colors: tuple[str, str], suffix: str) -> None:
    table = read_table(directory / "harmonics.csv")
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0), sharex="col")
    for col, orbit in enumerate(("passing", "trapped")):
        data = table[table["orbit"] == orbit]
        if len(data) == 0:
            raise ValueError(f"Missing {orbit} harmonics")
        for key in ("line_abs2", "boozer_abs2"):
            if np.any(data[key] <= 0):
                raise ValueError(f"{orbit}: zero harmonic power, cannot use log scale")
        if len(np.unique(data["m_b"])) != 1:
            raise ValueError(f"{orbit}: mixed harmonic indices")
        mb = data["m_b"][0]
        pitch = data["pitch_fraction"]
        ax = axes[0, col]
        ax.plot(pitch, data["line_abs2"], color=colors[0], linewidth=1.9,
                label="Line integral")
        ax.plot(pitch, data["boozer_abs2"], color=colors[1], linewidth=1.5,
                linestyle="--", marker="o", markersize=3.4,
                markerfacecolor="none", label="Boozer scalar")
        ax.set_yscale("log")
        ax.set_title(f"{orbit.capitalize()}, $m_b={mb}$")
        ax.set_ylabel(r"$|H_m/E|^2$ [1]")
        ax.grid(alpha=0.2)
        error_ax = axes[1, col]
        error_ax.plot(pitch, data["relative_error"], color=colors[0],
                      marker="o", markersize=3.4, linewidth=1.4)
        error_ax.set_yscale("symlog", linthresh=1e-12)
        error_ax.set_ylim(bottom=0)
        error_ax.set_ylabel(r"$|H_m^{line}-H_m^{Boozer}|/|H_m^{Boozer}|$ [1]")
        error_ax.set_xlabel(r"Pitch fraction [1]")
        error_ax.grid(alpha=0.2)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.95),
               ncol=2, frameon=False)
    fig.suptitle("Ideal-helical harmonic scan; resonant rotation at each pitch",
                 fontsize=11, y=0.99)
    fig.text(0.5, 0.015, "Pitch 0.5 is the trapped-passing boundary; "
             "error axes are linear below $10^{-12}$.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.89))
    save_figure(fig, directory, f"harmonics{suffix}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing S3 CSV files")
    parser.add_argument("--grayscale", action="store_true",
                        help="Render additional grayscale figures for inspection")
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    colors = ("#0072B2", "#E69F00")
    suffix = ""
    if args.grayscale:
        colors, suffix = ("0.15", "0.55"), "_grayscale"
    plot_integrands(args.directory, colors, suffix)
    plot_harmonics(args.directory, colors, suffix)
    print(f"Figures written to {args.directory.resolve()}")


if __name__ == "__main__":
    main()
