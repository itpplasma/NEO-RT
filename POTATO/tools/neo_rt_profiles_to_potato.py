#!/usr/bin/env python3
"""Convert NEO-RT s_tor profiles to POTATO polynomials in explicit s_pol."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

C_CGS = 2.99792458e10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "%")):
            continue
        try:
            rows.append([float(value) for value in stripped.split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"{path} contains no numeric rows")
    return rows


def fit_descending(x: np.ndarray, y: np.ndarray, degree: int) -> tuple[np.ndarray, float]:
    polynomial = Polynomial.fit(x, y, degree, domain=[0.0, 1.0]).convert()
    coefficients = np.pad(polynomial.coef, (0, degree + 1-polynomial.coef.size))
    scale = max(float(np.linalg.norm(y)), np.finfo(float).tiny)
    return coefficients[::-1], float(np.linalg.norm(polynomial(x)-y)/scale)


def convert(profile: Path, plasma: Path, geometry: Path, output: Path, *,
            r0_cm: float, psi_span_tm2: float, relation_sign: int,
            degree: int = 9) -> dict:
    if r0_cm <= 0.0 or psi_span_tm2 <= 0.0:
        raise ValueError("r0_cm and psi_span_tm2 must be positive")
    if relation_sign not in (-1, 0, 1):
        raise ValueError("relation_sign must be -1, 0, or +1")

    rotation = np.asarray(numeric_rows(profile), dtype=float)
    plasma_rows = numeric_rows(plasma)
    if rotation.shape[1] < 3 or len(plasma_rows) < 2:
        raise ValueError("truncated NEO-RT profile input")
    species = np.asarray(plasma_rows[0], dtype=float)
    if int(round(species[0])) != 50 or species.size < 5:
        raise ValueError("expected the NEO-RT two-ion header with leading grid count")
    charges = species[3:5]
    data = np.asarray(plasma_rows[1:], dtype=float)
    candidates = [index for index, charge in enumerate(charges)
                  if charge > 0.0 and np.any(data[:, 1+index] > 0.0)]
    if len(candidates) != 1:
        raise ValueError(f"cannot select one physical ion from charges {charges}")
    ion = candidates[0]

    with np.load(geometry) as mapping:
        s_pol_map = np.asarray(mapping["s_sqrt_poloidal"], dtype=float)**2
        s_tor_map = np.asarray(mapping["s_toroidal"], dtype=float)
    if np.any(np.diff(s_pol_map) <= 0.0) or np.any(np.diff(s_tor_map) <= 0.0):
        raise ValueError("geometry flux map must increase strictly")

    s_pol = np.linspace(0.0, 1.0, 1001)
    s_tor = np.interp(s_pol, s_pol_map, s_tor_map)
    density = np.interp(s_tor, data[:, 0], data[:, 1+ion])
    temperature_ev = np.interp(s_tor, data[:, 0], data[:, 3+ion])
    mach = np.interp(s_tor, rotation[:, 0], rotation[:, 1])
    vth_cm_s = np.interp(s_tor, rotation[:, 0], rotation[:, 2])
    omega_e_s = mach*vth_cm_s/r0_cm

    # Phi is statvolt and psi_pol is gauss cm^2. The sign is explicit because
    # NEO-RT and POTATO encode their poloidal/toroidal orientation differently.
    dphi_ds_pol = relation_sign*omega_e_s*(psi_span_tm2*1.0e8)/C_CGS
    potential_statv = np.zeros_like(s_pol)
    potential_statv[1:] = np.cumsum(
        0.5*(dphi_ds_pol[1:]+dphi_ds_pol[:-1])*np.diff(s_pol)
    )

    density_coef, density_error = fit_descending(s_pol, density, degree)
    temperature_coef, temperature_error = fit_descending(s_pol, temperature_ev, degree)
    potential_coef, potential_error = fit_descending(s_pol, potential_statv, degree)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as stream:
        stream.write("% NEO-RT profiles mapped from s_tor to s_pol\n")
        stream.write("% density, dummy, ion temperature, potential; descending powers\n")
        for values in (density_coef, np.zeros(degree+1),
                       temperature_coef, potential_coef):
            stream.write(" ".join(f"{value:.16e}" for value in values) + "\n")

    metadata = {
        "inputs": {
            "profile": {"path": str(profile), "sha256": sha256(profile)},
            "plasma": {"path": str(plasma), "sha256": sha256(plasma)},
            "geometry": {"path": str(geometry), "sha256": sha256(geometry)},
        },
        "output": {"path": str(output), "sha256": sha256(output)},
        "input_abscissa": "s_tor",
        "output_abscissa": "s_pol=rho_pol^2",
        "selected_ion_index_one_based": ion + 1,
        "selected_ion_charge": float(charges[ion]),
        "r0_cm": r0_cm,
        "psi_pol_span_tm2": psi_span_tm2,
        "potential_relation": "dPhi/dpsi_pol = relation_sign*Omega_E/c",
        "potential_relation_sign": relation_sign,
        "potential_units": "statvolt",
        "omega_e_s_range": [float(np.min(omega_e_s)), float(np.max(omega_e_s))],
        "potential_statv_edge": float(potential_statv[-1]),
        "polynomial_degree": degree,
        "relative_l2_fit": {
            "density": density_error,
            "temperature": temperature_error,
            "potential": potential_error,
        },
    }
    output.with_suffix(output.suffix + ".json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path)
    parser.add_argument("plasma", type=Path)
    parser.add_argument("geometry", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--r0-cm", type=float, required=True)
    parser.add_argument("--psi-span-tm2", type=float, required=True)
    parser.add_argument("--relation-sign", type=int, choices=(-1, 0, 1), required=True)
    parser.add_argument("--degree", type=int, default=9)
    args = parser.parse_args()
    print(json.dumps(convert(
        args.profile, args.plasma, args.geometry, args.output,
        r0_cm=args.r0_cm, psi_span_tm2=args.psi_span_tm2,
        relation_sign=args.relation_sign, degree=args.degree,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
