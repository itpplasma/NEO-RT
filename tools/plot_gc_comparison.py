#!/usr/bin/env python3
"""Plot Boozer-thin, full-FOW, and POTATO comparison artifacts.

The input directory is an accepted run directory containing the NEO-RT
transport files and the frequency-report tables.  The script does not alter
the numerical data; all signs and component labels are retained as written by
the executable.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class ModelPresentation:
    label: str
    color: str
    torque_magnitude_only: bool = False


# These are presentation contracts, not numerical model definitions.  The
# aliases reflect identifiers emitted by the current and accepted older
# diagnostics.  In particular, ``gc`` is model 2 (full real-space FOW), not
# the removed model-1 real-space-thin path.
MODEL_PRESENTATION = {
    "boozer_thin": ModelPresentation("Boozer thin", "#0072B2"),
    "full_fow": ModelPresentation("Full real-space FOW", "#D55E00"),
    "potato": ModelPresentation("POTATO", "#009E73", torque_magnitude_only=True),
}
MODEL_ORDER = tuple(MODEL_PRESENTATION)
MODEL_ALIASES = {
    "boozer_thin": "boozer_thin",
    "legacy": "boozer_thin",
    "legacy_boozer": "boozer_thin",
    "model0": "boozer_thin",
    "0": "boozer_thin",
    "full_fow": "full_fow",
    "gc": "full_fow",
    "gc_full": "full_fow",
    "neort_full": "full_fow",
    "model2": "full_fow",
    "2": "full_fow",
    "potato": "potato",
    "potato_full": "potato",
}
DEPRECATED_REAL_SPACE_THIN = {
    "gc_thin", "gc_thin_full", "real_space_thin", "realspace_thin", "model1", "1"
}
MODEL_FILE_STEMS = {
    "boozer_thin": ("boozer_thin", "legacy"),
    "full_fow": ("full_fow", "gc", "neort_full"),
    "potato": ("potato",),
}
SIGNED_TORQUE_MODELS = ("boozer_thin", "full_fow")
OMEGA_B_PROVENANCE = (
    r"EQDSK full-FOW: $\omega_b=2\pi/\tau_b>0$; signed $\omega_\phi$ retained."
)
POTATO_TORQUE_POLICY = (
    "POTATO torque weight is magnitude-only conjugate-pair mode power; "
    "it is not a signed torque component without an explicit sign gate."
)
POTATO_SIGN_GATE = (
    "sign gate OPEN for magnitude display only; CLOSED for a signed toroidal component"
)
POTATO_TORQUE_FILENAMES = (
    "integral_torque.dat", "potato_torque.out", "potato_torque.dat"
)


def canonical_model(identifier: str) -> str:
    """Return the presentation key, rejecting the removed real-space-thin ID."""
    normalized = identifier.strip().lower().replace("-", "_")
    if normalized in DEPRECATED_REAL_SPACE_THIN:
        raise ValueError(
            f"deprecated ambiguous model identifier {identifier!r}: "
            "real-space thin (frequency_model=1) was removed; use "
            "'boozer_thin' or 'full_fow'"
        )
    try:
        return MODEL_ALIASES[normalized]
    except KeyError as error:
        raise ValueError(
            f"unknown comparison model identifier {identifier!r}; expected "
            "boozer_thin, full_fow, or potato"
        ) from error


def model_presentation(model: str) -> ModelPresentation:
    return MODEL_PRESENTATION[canonical_model(model)]


def model_label(model: str, torque: bool = False) -> str:
    presentation = model_presentation(model)
    if torque and presentation.torque_magnitude_only:
        return "POTATO (magnitude-only; conjugate-pair mode power)"
    return presentation.label


def ordered_models(rows):
    present = {row[0] for row in rows}
    return [model for model in MODEL_ORDER if model in present]


def reject_deprecated_model_files(run_dir: Path) -> None:
    """Fail before plotting an artifact whose filename still means model 1."""
    for path in run_dir.iterdir():
        if not (path.name.lower().endswith((".out", ".time"))):
            continue
        stem = path.name.lower().replace("-", "_").split(".", 1)[0]
        if any(stem == candidate or stem.startswith(candidate + "_")
               for candidate in DEPRECATED_REAL_SPACE_THIN):
            raise ValueError(
                f"deprecated ambiguous real-space-thin artifact {path.name!r}; "
                "refusing to label it as full real-space FOW"
            )


def find_model_file(run_dir: Path, model: str, suffix: str) -> Path | None:
    """Find one unambiguous file for a canonical presentation model."""
    canonical = canonical_model(model)
    candidates = [run_dir / f"{stem}{suffix}"
                  for stem in MODEL_FILE_STEMS[canonical]]
    existing = [path for path in candidates if path.is_file()]
    if len(existing) > 1:
        names = ", ".join(path.name for path in existing)
        raise ValueError(f"ambiguous {canonical} artifacts: {names}")
    return existing[0] if existing else None


def discover_report(run_dir: Path, current_glob: str, legacy_name: str) -> Path:
    """Resolve one current or legacy report, rejecting ambiguous inventories."""
    current = sorted(run_dir.glob(current_glob))
    legacy = run_dir / legacy_name
    candidates = current + ([legacy] if legacy.is_file() else [])
    if not candidates:
        raise FileNotFoundError(
            f"no report in {run_dir}: expected {current_glob!r} or {legacy_name!r}"
        )
    if len(candidates) > 1:
        names = ", ".join(path.name for path in candidates)
        raise ValueError(f"ambiguous report in {run_dir}: {names}")
    return candidates[0]


def find_potato_torque_file(run_dir: Path) -> Path | None:
    candidates = [run_dir / name for name in POTATO_TORQUE_FILENAMES]
    existing = [path for path in candidates if path.is_file()]
    if len(existing) > 1:
        names = ", ".join(path.name for path in existing)
        raise ValueError(f"ambiguous POTATO torque artifacts: {names}")
    return existing[0] if existing else None


def read_potato_torque_scalar(path: Path) -> float:
    """Read the scalar POTATO integral schema; reject combined artifacts."""
    values = np.fromstring(path.read_text().replace("D", "E"), sep=" ")
    if values.size != 1 or not np.isfinite(values[0]):
        raise ValueError(
            f"incompatible POTATO torque artifact {path.name!r}: expected one "
            "finite scalar from integral_torque.dat"
        )
    return float(values[0])


def configure_signed_axis(ax, values) -> None:
    """Use a zero-containing signed scale without changing plotted values."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("cannot configure signed axis from an empty data set")
    largest = max(float(np.max(finite)), float(-np.min(finite)), 0.0)
    linthresh = largest * 1.0e-3 if largest > 0.0 else 1.0
    ax.set_yscale("symlog", linthresh=linthresh)
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.7)
    lower = min(float(np.min(finite)), 0.0)
    upper = max(float(np.max(finite)), 0.0)
    span = upper - lower
    if span == 0.0:
        span = max(linthresh, 1.0)
    padding = 0.05 * span
    ax.set_ylim(lower - padding, upper + padding)


def add_provenance(fig, include_potato_torque_policy: bool = False) -> None:
    note = OMEGA_B_PROVENANCE
    if include_potato_torque_policy:
        note += "  " + POTATO_TORQUE_POLICY
    fig.text(0.01, 0.006, note, ha="left", va="bottom", fontsize=7)


def read_frequency(path: Path):
    rows = []
    lines = path.read_text().splitlines()
    inventory = any(line.startswith("# frequency_inventory_v1") for line in lines)
    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if inventory:
            if len(fields) < 12:
                raise ValueError(f"malformed frequency inventory row in {path}: {line!r}")
            rows.append((canonical_model(fields[0]), fields[1], float(fields[5]),
                         float(fields[6]), float(fields[7]), float(fields[8]),
                         float(fields[9])))
            continue
        if len(fields) < 7:
            raise ValueError(f"malformed frequency row in {path}: {line!r}")
        rows.append((canonical_model(fields[0]), fields[1], *map(float, fields[2:])))
    return rows


def read_resonances(path: Path):
    rows = []
    lines = path.read_text().splitlines()
    inventory = any(line.startswith("# resonance_inventory_v1") for line in lines)
    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if inventory:
            if len(fields) < 12:
                raise ValueError(f"malformed resonance inventory row in {path}: {line!r}")
            rows.append((canonical_model(fields[0]), fields[1], int(fields[4]),
                         float(fields[7]), float(fields[8]), float(fields[9])))
            continue
        if len(fields) < 6:
            raise ValueError(f"malformed resonance row in {path}: {line!r}")
        rows.append((canonical_model(fields[0]), fields[1], int(fields[2]),
                     *map(float, fields[3:])))
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
        for model in ordered_models(data):
            for orbit_class in ("passing", "trapped"):
                rows = [row for row in data if row[0] == model and row[1] == orbit_class]
                if not rows:
                    continue
                x = np.array([row[2] for row in rows])
                y = np.array([row[column] for row in rows])
                style = "-" if orbit_class == "trapped" else "--"
                presentation = MODEL_PRESENTATION[model]
                ax.plot(x, y, style, color=presentation.color,
                        label=f"{model_label(model)} / {orbit_class}")
        ax.set_ylabel(ylabel)
        configure_axes(ax)
    axes[1, 0].set_xlabel(r"$\eta/\eta_{tp}$")
    axes[1, 1].set_xlabel(r"$\eta/\eta_{tp}$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.985),
               ncol=2, frameon=False)
    fig.suptitle("Boozer-thin, full real-space FOW, and POTATO frequencies", y=1.06)
    add_provenance(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 0.88))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_resonances(data, output):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    markers = {"passing": "o", "trapped": "s"}
    for model in ordered_models(data):
        for orbit_class in ("passing", "trapped"):
            rows = [row for row in data if row[0] == model and row[1] == orbit_class]
            if not rows:
                continue
            presentation = MODEL_PRESENTATION[model]
            ax.scatter([row[3] for row in rows], [row[2] for row in rows],
                       marker=markers[orbit_class], color=presentation.color,
                       edgecolors="black", linewidths=0.35, s=50,
                       label=f"{model_label(model)} / {orbit_class}")
        
    ax.axvline(1.0, color="black", linewidth=0.7, alpha=0.5)
    ax.set_xlabel(r"resonant $\eta/\eta_{tp}$")
    ax.set_ylabel(r"poloidal harmonic $m_{th}$")
    ax.set_title(r"Roots of $m_{th}\omega_b+3\omega_\phi=0$")
    configure_axes(ax)
    ax.legend(frameon=False, ncol=2)
    add_provenance(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_transport(run_dir, output):
    rows = {}
    for model in SIGNED_TORQUE_MODELS:
        path = find_model_file(run_dir, model, ".out")
        if path is not None:
            rows[model] = read_numeric_row(path)
    if not rows:
        raise FileNotFoundError("no Boozer-thin or full-FOW transport artifact found")
    labels = ["D11 co", "D11 ctr", "D11 trapped", "D11 total",
              "D12 co", "D12 ctr", "D12 trapped", "D12 total"]
    positions = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for offset, model in zip(
        np.linspace(-width * (len(rows) - 1) / 2,
                    width * (len(rows) - 1) / 2, len(rows)), rows
    ):
        ax.bar(positions + offset, rows[model][1:], width,
               color=MODEL_PRESENTATION[model].color,
               label=MODEL_PRESENTATION[model].label)
    ax.set_xticks(positions, labels, rotation=35, ha="right")
    ax.set_ylabel("coefficient")
    ax.set_title("Transport coefficients at s=0.25")
    configure_signed_axis(ax, np.concatenate([values[1:] for values in rows.values()]))
    configure_axes(ax)
    ax.legend(frameon=False)
    add_provenance(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def read_signed_torque_values(run_dir: Path):
    values = {}
    for model in SIGNED_TORQUE_MODELS:
        path = find_model_file(run_dir, model, "_torque.out")
        if path is None:
            raise FileNotFoundError(
                f"missing {model} torque artifact; all signed comparison models "
                "are required"
            )
        row = read_numeric_row(path)
        if row.size != 6:
            raise ValueError(
                f"incompatible {path.name!r}: expected six native torque columns, "
                f"got {row.size}"
            )
        values[model] = row[3:]
    return values


def plot_potato_magnitude_axis(ax, raw_value: float, title: str) -> None:
    # This is an explicitly labelled display magnitude for POTATO's
    # conjugate-pair weight.  The raw native value remains in provenance.
    magnitude = max(raw_value, -raw_value)
    ax.bar(["POTATO"], [magnitude], color=MODEL_PRESENTATION["potato"].color,
           edgecolor="black", hatch="//", label="POTATO magnitude-only")
    ax.set_ylim(0.0, max(1.0, 1.15 * magnitude))
    ax.set_ylabel("magnitude [native POTATO units]")
    ax.set_title(title)
    ax.text(0.02, 0.98, POTATO_SIGN_GATE, transform=ax.transAxes,
            ha="left", va="top", fontsize=8, wrap=True)
    configure_axes(ax)
    ax.legend(frameon=False, loc="upper right")


def make_torque_figure(run_dir: Path):
    values = read_signed_torque_values(run_dir)
    potato_path = find_potato_torque_file(run_dir)
    if potato_path is None:
        raise FileNotFoundError(
            "missing POTATO torque artifact; expected integral_torque.dat "
            "or an explicitly named potato torque scalar"
        )
    potato_raw = read_potato_torque_scalar(potato_path)
    labels = ["co", "counter", "trapped"]
    positions = np.arange(len(labels))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True,
                             gridspec_kw={"width_ratios": (2.0, 1.0)})
    ax = axes[0]
    for offset, model in zip(
        np.linspace(-width * (len(values) - 1) / 2,
                    width * (len(values) - 1) / 2, len(values)), values
    ):
        ax.bar(positions + offset, values[model], width,
               color=MODEL_PRESENTATION[model].color,
               label=MODEL_PRESENTATION[model].label)
    ax.set_xticks(positions, labels)
    ax.set_ylabel("torque density [native executable units]")
    ax.set_title("Torque density at s=0.25; native signs retained")
    configure_signed_axis(ax, np.concatenate(list(values.values())))
    configure_axes(ax)
    ax.legend(frameon=False)
    plot_potato_magnitude_axis(
        axes[1], potato_raw,
        "POTATO integrated torque magnitude; not torque density"
    )
    add_provenance(fig, include_potato_torque_policy=True)
    return fig


def plot_torque(run_dir, output):
    fig = make_torque_figure(run_dir)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_torque_spectrum(run_dir, output):
    potato_path = find_potato_torque_file(run_dir)
    if potato_path is None:
        raise FileNotFoundError(
            "missing POTATO torque artifact; expected integral_torque.dat "
            "or an explicitly named potato torque scalar"
        )
    potato_raw = read_potato_torque_scalar(potato_path)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), constrained_layout=True,
                             gridspec_kw={"width_ratios": (2.0, 1.0)})
    ax = axes[0]
    plotted = []
    for model in SIGNED_TORQUE_MODELS:
        path = find_model_file(run_dir, model, "_torque_integral.out")
        if path is None:
            raise FileNotFoundError(
                f"missing {model} torque spectrum; all signed comparison models "
                "are required"
            )
        table = read_numeric_table(path)
        if table.ndim != 2 or table.shape[1] < 2:
            raise ValueError(
                f"incompatible {path.name!r}: expected harmonic plus torque "
                f"columns, got shape {table.shape}"
            )
        mth = table[:, 0]
        total = table[:, 1:].sum(axis=1)
        plotted.append(model)
        ax.plot(mth, total, marker="o", color=MODEL_PRESENTATION[model].color,
                label=MODEL_PRESENTATION[model].label)
    if not plotted:
        raise FileNotFoundError("no Boozer-thin or full-FOW torque spectrum found")
    ax.set_xlabel(r"$m_{th}$")
    ax.set_ylabel("torque integral by harmonic [native units]")
    ax.set_title("Torque spectrum")
    totals = []
    for model in plotted:
        path = find_model_file(run_dir, model, "_torque_integral.out")
        totals.append(read_numeric_table(path)[:, 1:].sum(axis=1))
    configure_signed_axis(ax, np.concatenate(totals))
    configure_axes(ax)
    ax.legend(frameon=False)
    plot_potato_magnitude_axis(
        axes[1], potato_raw,
        "POTATO integrated torque magnitude; no harmonic decomposition"
    )
    add_provenance(fig, include_potato_torque_policy=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_timing(run_dir, output, summary):
    entries = []
    for model in MODEL_ORDER:
        for omp in (1, 4):
            path = find_model_file(run_dir, model, f"-omp{omp}.time")
            if path is None:
                continue
            timing = read_timing(path)
            entries.append((model, f"{MODEL_PRESENTATION[model].label}\nOMP={omp}",
                            timing["wall_seconds"],
                            timing["rss_kib"] / 1024.0))
            summary.setdefault("timing", {})[f"{model}_omp{omp}"] = timing
    if not entries:
        raise FileNotFoundError("no timing artifacts found")
    positions = np.arange(len(entries))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = [MODEL_PRESENTATION[entry[0]].color for entry in entries]
    bars = ax.bar(positions, [entry[2] for entry in entries], color=colors)
    ax.set_xticks(positions, [entry[1] for entry in entries])
    ax.set_ylabel("wall time [s]")
    ax.set_title("End-to-end timing benchmark")
    configure_axes(ax)
    for bar, entry in zip(bars, entries):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{entry[2]:.2f}s\n{entry[3]:.1f} MiB", ha="center", va="bottom", fontsize=8)
    add_provenance(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frequency-dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reject_deprecated_model_files(args.run_dir)
    summary = {"run_dir": str(args.run_dir)}
    frequency_dir = args.frequency_dir or args.run_dir
    if frequency_dir != args.run_dir:
        reject_deprecated_model_files(frequency_dir)
    frequency_path = discover_report(
        frequency_dir, "*_frequency_inventory.out", "frequency_scan.out"
    )
    resonance_path = discover_report(
        frequency_dir, "*_resonance_inventory.out", "frequency_resonances.out"
    )
    frequency = read_frequency(frequency_path)
    resonances = read_resonances(resonance_path)
    plot_frequencies(frequency, args.output_dir / "frequencies.png")
    plot_resonances(resonances, args.output_dir / "resonances.png")
    plot_transport(args.run_dir, args.output_dir / "transport.png")
    plot_torque(args.run_dir, args.output_dir / "torque_density.png")
    plot_torque_spectrum(args.run_dir, args.output_dir / "torque_spectrum.png")
    plot_timing(args.run_dir, args.output_dir / "timing.png", summary)

    summary["provenance"] = {
        "omega_b_convention": (
            "EQDSK full-FOW provider uses omega_b=2*pi/tau_b > 0; "
            "signed omega_phi is retained"
        ),
        "potato_torque_policy": POTATO_TORQUE_POLICY,
        "signed_torque_models": list(SIGNED_TORQUE_MODELS),
    }
    summary["model_labels"] = {
        model: model_label(model) for model in MODEL_ORDER
    }
    for model in MODEL_ORDER:
        transport_path = find_model_file(args.run_dir, model, ".out")
        if transport_path is not None:
            transport_target = "transport" if model in SIGNED_TORQUE_MODELS \
                else "potato_transport"
            summary.setdefault(transport_target, {})[model] = \
                read_numeric_row(transport_path).tolist()
        if model in SIGNED_TORQUE_MODELS:
            torque_path = find_model_file(args.run_dir, model, "_torque.out")
            if torque_path is not None:
                summary.setdefault("torque_density", {})[model] = \
                    read_numeric_row(torque_path).tolist()
    potato_torque_path = find_potato_torque_file(args.run_dir)
    if potato_torque_path is not None:
        summary["potato_torque"] = {
            "raw_native_value": read_potato_torque_scalar(potato_torque_path),
            "displayed_quantity": "magnitude-only",
            "sign_gate": POTATO_SIGN_GATE,
        }
    summary["resonance_rows"] = len(resonances)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
