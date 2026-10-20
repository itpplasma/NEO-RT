#!/usr/bin/env python3
"""Plot the circular Kasilov offset-rotation behavioral scan."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


HERE = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=HERE / "circular_offset_rotation.csv")
    parser.add_argument("--output", type=Path, default=HERE / "circular_offset_rotation.png")
    parser.add_argument("--pdf", type=Path, default=HERE / "circular_offset_rotation.pdf")
    parser.add_argument("--grayscale", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = np.genfromtxt(args.input, delimiter=",", names=True)
    delta = data["delta_omega_s1"]
    torque = data["torque_native"]
    residual = data["residual_s1"]

    blue = "0.20" if args.grayscale else "#0072B2"
    orange = "0.45" if args.grayscale else "#D55E00"

    fig, (ax_torque, ax_residual) = plt.subplots(
        2, 1, figsize=(6.4, 6.4), sharex=True, constrained_layout=True
    )

    ax_torque.axhline(0.0, color="0.25", linewidth=0.9)
    ax_torque.axvline(0.0, color="0.55", linewidth=0.9, linestyle=":")
    ax_torque.plot(
        delta,
        torque,
        color=blue,
        linewidth=1.7,
        marker="o",
        markersize=4.5,
        label="NEO-RT scan",
    )
    ax_torque.scatter(
        [delta[0], delta[len(delta) // 2], delta[-1]],
        [torque[0], torque[len(torque) // 2], torque[-1]],
        color=blue,
        edgecolor="white",
        linewidth=0.8,
        marker="D",
        s=42,
        zorder=3,
        label="below / offset / above",
    )
    ax_torque.set_ylabel("NTV torque [native units]")
    ax_torque.set_title(r"Circular Kasilov offset: $k=1.17$, $s=0.075$")
    ax_torque.legend(frameon=False, loc="best")
    ax_torque.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax_torque.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))

    ax_residual.axhline(0.0, color="0.25", linewidth=0.9)
    ax_residual.axvline(0.0, color="0.55", linewidth=0.9, linestyle=":")
    ax_residual.plot(
        delta,
        residual,
        color=orange,
        linewidth=1.7,
        linestyle="--",
        marker="s",
        markersize=4.2,
        label=r"$V^\phi-V^\phi_{\rm in}$",
    )
    ax_residual.plot(
        delta,
        delta,
        color="0.25",
        linewidth=1.0,
        linestyle=":",
        label="unit-slope reference",
    )
    ax_residual.set_xlabel(
        r"imposed shift from the offset, $\Delta\Omega_{tE}=\Delta V^\phi$ [$s^{-1}$]"
    )
    ax_residual.set_ylabel(r"fixed-point residual [$s^{-1}$]")
    ax_residual.legend(frameon=False, loc="best")

    for axis in (ax_torque, ax_residual):
        axis.grid(True, color="0.88", linewidth=0.7)
        axis.set_axisbelow(True)

    fig.savefig(args.output, dpi=200)
    fig.savefig(args.pdf)


if __name__ == "__main__":
    main()
