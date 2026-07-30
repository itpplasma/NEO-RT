#!/usr/bin/env python3
"""Plot ITER shielding current independently reconstructed from MARS B."""

from __future__ import annotations

import argparse
import csv
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


FIELD_STYLES = {
    "total": ("#000000", "-"),
    "vacuum": ("#0072B2", ":"),
    "plasma_response": ("#E69F00", "--"),
}
MODE_COLORS = (
    "#0072B2",
    "#E69F00",
    "#CC79A7",
    "#009E73",
    "#D55E00",
    "#56B4E9",
    "#6F4E7C",
)


def spectral_norm(coefficients: np.ndarray) -> np.ndarray:
    """Poloidal-harmonic L2 amplitude on every radial surface."""
    return np.linalg.norm(coefficients, axis=1)


def rational_surfaces(
    stor: np.ndarray, q: np.ndarray, toroidal_mode: int
) -> dict[int, float]:
    """Locate all integer-m resonances inside the supplied q interval."""
    if toroidal_mode == 0:
        return {}
    start = int(np.argmin(q))
    q_branch = q[start:]
    stor_branch = stor[start:]
    if np.any(np.diff(q_branch) <= 0.0):
        raise ValueError("q is not monotone outside its axis minimum")
    n_abs = abs(toroidal_mode)
    first = int(np.ceil(n_abs * q_branch[0]))
    last = int(np.floor(n_abs * q_branch[-1]))
    result = {}
    for mode in range(first, last + 1):
        target = mode / n_abs
        result[mode] = float(np.interp(target, q_branch, stor_branch))
    return result


def load_currents(manifest: Path, nchi: int):
    """Return grids, q, modes, and native curl(B) for three field pieces."""
    fields = {}
    provenance = {}
    reference = None
    for component in ("total", "vacuum", "plasma_response"):
        geometry, _, vector, dpsi, field, inputs = load_manifest(
            manifest, component, center_half_mesh=False
        )
        radial, current = reconstruct_current(geometry, vector, nchi)
        chi = np.linspace(-np.pi, np.pi, nchi, endpoint=False)
        stor, q = toroidal_flux_coordinate(geometry, dpsi, field, chi)
        signature = (radial, stor, q, vector.modes, vector.toroidal_mode)
        if reference is None:
            reference = signature
        else:
            for actual, expected in zip(signature[:-1], reference[:-1], strict=True):
                if not np.array_equal(actual, expected):
                    raise ValueError("MARS field pieces use different grids")
            if signature[-1] != reference[-1]:
                raise ValueError("MARS field pieces use different toroidal modes")
        fields[component] = current
        provenance[component] = inputs
    assert reference is not None
    return (*reference, fields, provenance)


def profile_rows(
    stor: np.ndarray,
    modes: np.ndarray,
    fields: dict[str, tuple[np.ndarray, ...]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    half_stor = 0.5 * (stor[:-1] + stor[1:])
    grids = (half_stor, stor, stor)
    for component_index, component_name in enumerate(("J1", "J2", "J3")):
        for field_name, current in fields.items():
            coefficients = current[component_index]
            valid = np.all(np.isfinite(coefficients), axis=1)
            rows.append(
                {
                    "component": component_name,
                    "field": field_name,
                    "stor": grids[component_index][valid],
                    "amplitude": spectral_norm(coefficients[valid]),
                    "coefficients": coefficients[valid],
                    "modes": modes,
                }
            )
    return rows


def render_profiles(
    rows: list[dict[str, object]], output: Path, case_name: str
) -> None:
    fig, axes = plt.subplots(
        3, 1, figsize=(7.2, 8.2), sharex=True, constrained_layout=True
    )
    coordinates = {"J1": "s", "J2": r"\chi", "J3": r"\phi"}
    for axis, component in zip(axes, coordinates, strict=True):
        for row in rows:
            if row["component"] != component:
                continue
            field = str(row["field"])
            color, style = FIELD_STYLES[field]
            axis.semilogy(
                row["stor"],
                row["amplitude"],
                color=color,
                linestyle=style,
                linewidth=1.6,
                label=field.replace("_", " "),
            )
        axis.set_ylabel(
            rf"$\|\sqrt{{g}}J^{{{coordinates[component]}}}\|_m$"
            + "\n(MARS units)"
        )
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(loc="upper center", ncol=3, frameon=False, fontsize=8)
    axes[0].set_title(
        f"Native-staggered current from MARS magnetic field — {case_name}"
    )
    axes[-1].set_xlabel(r"normalized toroidal flux $s_{\rm tor}$")
    axes[-1].set_xlim(0.0, 1.0)
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def render_resonant(
    stor: np.ndarray,
    modes: np.ndarray,
    response: tuple[np.ndarray, ...],
    toroidal_mode: int,
    resonances: dict[int, float],
    output: Path,
    case_name: str,
) -> list[tuple[int, np.ndarray, np.ndarray]]:
    valid = np.all(np.isfinite(response[1]), axis=1)
    valid &= np.all(np.isfinite(response[2]), axis=1)
    radial = stor[valid]
    curves = []
    fig, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for index, (mode, surface) in enumerate(resonances.items()):
        mode_index = np.flatnonzero(modes == mode)
        if len(mode_index) != 1:
            raise ValueError(f"MARS spectrum does not contain resonant m={mode}")
        k = int(mode_index[0])
        amplitude = np.sqrt(
            np.abs(response[1][valid, k]) ** 2
            + np.abs(response[2][valid, k]) ** 2
        )
        color = MODE_COLORS[index % len(MODE_COLORS)]
        style = ("-", "--", ":", "-.")[index % 4]
        axis.semilogy(
            radial,
            amplitude,
            color=color,
            linestyle=style,
            linewidth=1.6,
            label=rf"$m={mode}$",
        )
        axis.axvline(
            surface, color=color, linestyle=":", linewidth=0.8, alpha=0.55
        )
        curves.append((mode, radial, amplitude))
    axis.set_title(
        f"Plasma-response current at n={toroidal_mode} "
        f"rational harmonics — {case_name}"
    )
    axis.set_xlabel(r"normalized toroidal flux $s_{\rm tor}$")
    axis.set_ylabel(
        r"$\sqrt{|\,\sqrt{g}J^\chi_m|^2+|\,\sqrt{g}J^\phi_m|^2}$"
        "\n(MARS units)"
    )
    axis.set_xlim(0.0, 1.0)
    axis.grid(True, which="both", alpha=0.25)
    axis.legend(ncol=4, loc="upper center", frameon=False, fontsize=8)
    axis.text(
        0.99,
        0.92,
        r"dotted vertical: $q=m/|n|$",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
    )
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)
    return curves


def write_profile_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("component", "field", "s_tor", "spectral_norm"))
        for row in rows:
            for radial, amplitude in zip(
                row["stor"], row["amplitude"], strict=True
            ):
                writer.writerow(
                    (row["component"], row["field"], radial, amplitude)
                )


def write_resonant_csv(
    path: Path,
    curves: list[tuple[int, np.ndarray, np.ndarray]],
    resonances: dict[int, float],
) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("m", "s_tor", "tangential_current", "rational_s_tor"))
        for mode, radial, amplitude in curves:
            for s_value, value in zip(radial, amplitude, strict=True):
                writer.writerow((mode, s_value, value, resonances[mode]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("profile_output", type=Path)
    parser.add_argument("resonant_output", type=Path)
    parser.add_argument("--nchi", type=int, default=1024)
    parser.add_argument("--case-name", default="ITER TC24 phiI000")
    args = parser.parse_args()

    radial, stor, q, modes, toroidal_mode, fields, provenance = load_currents(
        args.manifest, args.nchi
    )
    rows = profile_rows(stor, modes, fields)
    resonances = rational_surfaces(stor, q, toroidal_mode)
    curves = render_resonant(
        stor,
        modes,
        fields["plasma_response"],
        toroidal_mode,
        resonances,
        args.resonant_output,
        args.case_name,
    )
    render_profiles(rows, args.profile_output, args.case_name)
    write_profile_csv(args.profile_output.with_suffix(".csv"), rows)
    write_resonant_csv(
        args.resonant_output.with_suffix(".csv"), curves, resonances
    )

    flat = {
        name: np.concatenate(
            [component[np.isfinite(component).all(axis=1)].ravel() for component in value]
        )
        for name, value in fields.items()
    }
    diagnostics = {
        "case_name": args.case_name,
        "nchi": args.nchi,
        "toroidal_mode": toroidal_mode,
        "surfaces": len(radial),
        "poloidal_mode_min": int(modes[0]),
        "poloidal_mode_max": int(modes[-1]),
        "vacuum_over_total_current_l2": float(
            np.linalg.norm(flat["vacuum"]) / np.linalg.norm(flat["total"])
        ),
        "linearity_closure_relative_l2": float(
            np.linalg.norm(
                flat["total"] - flat["vacuum"] - flat["plasma_response"]
            )
            / np.linalg.norm(flat["total"])
        ),
        "rational_surfaces_s_tor": {
            str(mode): value for mode, value in resonances.items()
        },
        "provenance": provenance["total"],
        "radial_derivative": (
            "cubic derivative of natural MARS B2/B3 half-mesh values; "
            "no pre-curl collocation"
        ),
    }
    args.profile_output.with_suffix(".json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n"
    )
    args.resonant_output.with_suffix(".json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(diagnostics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
