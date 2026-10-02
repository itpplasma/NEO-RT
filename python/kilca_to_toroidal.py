#!/usr/bin/env python3
"""Export KiLCA magnetic or paired laboratory harmonics through an explicit chart."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from periodic_cylinder import (NativeMode, RadialGauge, AxialGauge, laboratory_omega,
                               read_geometry, transform_mode, transform_electromagnetic_mode)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("eb", type=Path)
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--R0-cm", type=float, required=True)
    parser.add_argument("--r-min-cm", type=float)
    parser.add_argument("--r-max-cm", type=float)
    parser.add_argument("--geometry", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-relative-br-residual", type=float)
    parser.add_argument("--gauge", choices=("radial-magnetic", "axial-electromagnetic"),
                        default="radial-magnetic")
    parser.add_argument("--mode-data", type=Path)
    parser.add_argument("--max-relative-B-residual", type=float)
    parser.add_argument("--max-relative-E-residual", type=float)
    parser.add_argument("--max-relative-Faraday-residual", type=float)
    parser.add_argument("--max-relative-div-B", type=float)
    parser.add_argument("--diagnostics-only", action="store_true")
    args = parser.parse_args()
    mode = NativeMode.from_eb(args.eb, m=args.m, n=args.n, R0=args.R0_cm,
                             r_min=args.r_min_cm, r_max=args.r_max_cm)
    residual = RadialGauge(mode).relative_radial_residual()
    metadata = {"schema": "periodic-cylinder-toroidal-samples-v1",
                "input_EB_sha256": digest(args.eb), "m": mode.m, "n": mode.n,
                "R0_cm": mode.R0, "source_r_range_cm": [float(mode.r[0]), float(mode.r[-1])],
                "max_relative_Br_residual": residual,
                "native_phase": "exp(i*(m*theta+n*z/R0))",
                "target_phase": "exp(i*n*phi)",
                "component_order": ["R", "phi", "Z"],
                "length_units": "cm", "A_units": "G cm", "B_units": "G",
                "electric_field_role": "finite-frequency diagnostic; no deltaPhi inferred",
                "potential_gauge": "A_r=0, A_z(r_min)=0, A_theta(r_min)=i*B_r(r_min)/(n/R0)"}
    omega = None
    if args.gauge == "axial-electromagnetic":
        mode_data = args.mode_data or args.eb.with_name("mode_data.dat")
        omega = laboratory_omega(mode_data, m=mode.m, n=mode.n)
        metadata["mode_data_sha256"] = digest(mode_data)
        metadata["omega_lab_rad_per_s"] = [float(np.real(omega)), float(np.imag(omega))]
        metadata["potential_gauge"] = "A_z=0; A_r=-i*B_theta/kz, A_theta=i*B_r/kz, Phi=i*E_z/kz"
        metadata["electric_field_role"] = "paired laboratory A/Phi; exp(-i*omega_lab*t)"
        metadata["native_EB_diagnostics"] = AxialGauge(mode, omega).diagnostics()
    if args.diagnostics_only:
        print(json.dumps(metadata, indent=2))
        return 0
    if args.geometry is None or args.output is None:
        parser.error("export requires --geometry and --output")
    geometry = read_geometry(args.geometry)
    if args.gauge == "radial-magnetic":
        if args.max_relative_br_residual is None:
            parser.error("magnetic export requires --max-relative-br-residual")
        fields = transform_mode(mode, geometry,
                                max_relative_br_residual=args.max_relative_br_residual)
    else:
        limits = {key: getattr(args, key) for key in ("max_relative_B_residual",
                  "max_relative_E_residual", "max_relative_Faraday_residual", "max_relative_div_B")}
        if any(value is None for value in limits.values()):
            parser.error("electromagnetic export requires all four E/B/Faraday/divergence limits")
        fields = transform_electromagnetic_mode(mode, geometry, omega_lab=omega, **limits)
    metadata["geometry_sha256"] = digest(args.geometry)
    with np.load(args.geometry, allow_pickle=False) as source:
        metadata["geometry_provenance"] = str(source["provenance"].item())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as output:
        np.savez_compressed(output, **fields, metadata_json=json.dumps(metadata, sort_keys=True))
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
