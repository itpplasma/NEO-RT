#!/usr/bin/env python3
"""Plot legacy versus real-space GC-thin comparison artifacts.

The input directory is an accepted run directory containing the NEO-RT
transport files and the frequency-report tables.  The script does not alter
the numerical data; all signs and component labels are retained as written by
the executable.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MODEL_LABEL = {"legacy": "Legacy", "gc_thin": "GC thin", "gc": "GC thin"}
MODEL_COLOR = {"legacy": "#355c7d", "gc_thin": "#c06c84", "gc": "#c06c84"}


def read_frequency(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        rows.append(
            (fields[0], fields[1], *map(float, fields[2:]))
        )
    return rows


def read_resonances(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        rows.append((fields[0], fields[1], int(fields[2]), *map(float, fields[3:])))
    return rows


def read_numeric_row(path: Path):
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    return np.fromstring(lines[-1].replace("D", "E"), sep=" ")


def read_numeric_table(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        rows.append(np.fromstring(line.replace("D", "E"), sep=" "))
    return np.asarray(rows)


def read_timing(path: Path):
    text = path.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time .*: ([^\n]+)", text)
    user = re.search(r"User time \(seconds\): ([^\n]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): ([^\n]+)", text)
    elapsed = wall.group(1).strip() if wall else "unknown"
    parts = elapsed.split(":")
    if len(parts) == 3:
        wall_seconds = 3600 * float(parts[0]) + 60 * float(parts[1]) + float(parts[2])
    elif len(parts) == 2:
        wall_seconds = 60 * float(parts[0]) + float(parts[1])
    else:
        wall_seconds = float(parts[0])
    return {
        "wall_seconds": wall_seconds,
        "user_seconds": float(user.group(1)) if user else float("nan"),
        "rss_kib": int(rss.group(1)) if rss else -1,
    }


def configure_axes(ax):
    ax.grid(True, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_frequencies(data, output):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    columns = [(3, r"$\omega_b$ [s$^{-1}$]"), (4, r"$\omega_\phi$ [s$^{-1}$]"),
               (5, r"$\omega_B$ [s$^{-1}$]"), (6, r"$\omega_E$ [s$^{-1}$]")]
    for ax, (column, ylabel) in zip(axes.flat, columns):
        for model in ("legacy", "gc_thin"):
            for orbit_class in ("passing", "trapped"):
                rows = [row for row in data if row[0] == model and row[1] == orbit_class]
                if not rows:
                    continue
                x = np.array([row[2] for row in rows])
                y = np.array([row[column] for row in rows])
                style = "-" if orbit_class == "trapped" else "--"
                ax.plot(x, y, style, color=MODEL_COLOR[model],
                        label=f"{MODEL_LABEL[model]} / {orbit_class}")
        ax.set_ylabel(ylabel)
        configure_axes(ax)
    axes[1, 0].set_xlabel(r"$\eta/\eta_{tp}$")
    axes[1, 1].set_xlabel(r"$\eta/\eta_{tp}$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.985),
               ncol=2, frameon=False)
    fig.suptitle("Direct real-space frequency comparison", y=1.06)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_resonances(data, output):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    markers = {"passing": "o", "trapped": "s"}
    for model in ("legacy", "gc_thin"):
        for orbit_class in ("passing", "trapped"):
            rows = [row for row in data if row[0] == model and row[1] == orbit_class]
            if not rows:
                continue
            ax.scatter([row[3] for row in rows], [row[2] for row in rows],
                       marker=markers[orbit_class], color=MODEL_COLOR[model],
                       edgecolors="black", linewidths=0.35, s=50,
                       label=f"{MODEL_LABEL[model]} / {orbit_class}")
        
    ax.axvline(1.0, color="black", linewidth=0.7, alpha=0.5)
    ax.set_xlabel(r"resonant $\eta/\eta_{tp}$")
    ax.set_ylabel(r"poloidal harmonic $m_{th}$")
    ax.set_title(r"Roots of $m_{th}\omega_b+3\omega_\phi=0$")
    configure_axes(ax)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_transport(run_dir, output):
    rows = {}
    for model in ("legacy", "gc"):
        values = read_numeric_row(run_dir / f"{model}.out")
        rows[model] = values
    labels = ["D11 co", "D11 ctr", "D11 trapped", "D11 total",
              "D12 co", "D12 ctr", "D12 trapped", "D12 total"]
    positions = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for offset, model in [(-width / 2, "legacy"), (width / 2, "gc")]:
        ax.bar(positions + offset, rows[model][1:], width,
               color=MODEL_COLOR[model], label=MODEL_LABEL[model])
    ax.set_xticks(positions, labels, rotation=35, ha="right")
    ax.set_ylabel("coefficient")
    ax.set_title("Transport coefficients at s=0.25")
    ax.set_yscale("log")
    configure_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_torque(run_dir, output):
    values = {}
    for model in ("legacy", "gc"):
        values[model] = read_numeric_row(run_dir / f"{model}_torque.out")[3:]
    labels = ["co", "counter", "trapped"]
    positions = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for offset, model in [(-width / 2, "legacy"), (width / 2, "gc")]:
        ax.bar(positions + offset, values[model], width,
               color=MODEL_COLOR[model], label=MODEL_LABEL[model])
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(positions, labels)
    ax.set_ylabel("torque density [native executable units]")
    ax.set_title("Torque density at s=0.25; native signs retained")
    ax.set_yscale("symlog", linthresh=1.0)
    configure_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_torque_spectrum(run_dir, output):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for model in ("legacy", "gc"):
        table = read_numeric_table(run_dir / f"{model}_torque_integral.out")
        mth = table[:, 0]
        total = table[:, 1:].sum(axis=1)
        ax.plot(mth, total, marker="o", color=MODEL_COLOR[model],
                label=MODEL_LABEL[model])
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel(r"$m_{th}$")
    ax.set_ylabel("torque integral by harmonic [native units]")
    ax.set_title("Torque spectrum")
    ax.set_yscale("symlog", linthresh=100.0)
    configure_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_timing(run_dir, output, summary):
    entries = []
    for model in ("legacy", "gc"):
        for omp in (1, 4):
            timing = read_timing(run_dir / f"{model}-omp{omp}.time")
            entries.append((f"{MODEL_LABEL[model]}\nOMP={omp}", timing["wall_seconds"],
                            timing["rss_kib"] / 1024.0))
            summary.setdefault("timing", {})[f"{model}_omp{omp}"] = timing
    positions = np.arange(len(entries))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = [MODEL_COLOR["legacy"] if "Legacy" in entry[0]
              else MODEL_COLOR["gc_thin"] for entry in entries]
    bars = ax.bar(positions, [entry[1] for entry in entries], color=colors)
    ax.set_xticks(positions, [entry[0] for entry in entries])
    ax.set_ylabel("wall time [s]")
    ax.set_title("End-to-end timing benchmark on faepkub4")
    configure_axes(ax)
    for bar, entry in zip(bars, entries):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{entry[1]:.2f}s\n{entry[2]:.1f} MiB", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frequency-dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = {"run_dir": str(args.run_dir)}
    frequency_dir = args.frequency_dir or args.run_dir
    frequency = read_frequency(frequency_dir / "frequency_scan.out")
    resonances = read_resonances(frequency_dir / "frequency_resonances.out")
    plot_frequencies(frequency, args.output_dir / "frequencies.png")
    plot_resonances(resonances, args.output_dir / "resonances.png")
    plot_transport(args.run_dir, args.output_dir / "transport.png")
    plot_torque(args.run_dir, args.output_dir / "torque_density.png")
    plot_torque_spectrum(args.run_dir, args.output_dir / "torque_spectrum.png")
    plot_timing(args.run_dir, args.output_dir / "timing.png", summary)

    for model in ("legacy", "gc"):
        transport = read_numeric_row(args.run_dir / f"{model}.out")
        torque = read_numeric_row(args.run_dir / f"{model}_torque.out")
        summary.setdefault("transport", {})[model] = transport.tolist()
        summary.setdefault("torque_density", {})[model] = torque.tolist()
    summary["resonance_rows"] = len(resonances)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
