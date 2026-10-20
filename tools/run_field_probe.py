#!/usr/bin/env python3
"""Run the native field-only probe and record exact file identities."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--bc", type=Path, required=True)
    parser.add_argument("--points", type=Path, required=True,
                        help="CSV with s_tor,phi_B_rad,theta_B_rad or whitespace triples")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inp-swi", type=int, choices=(8, 9), default=9)
    args = parser.parse_args()
    args.probe, args.bc, args.points = (
        path.resolve(strict=True) for path in (args.probe, args.bc, args.points)
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    query = output.with_suffix(output.suffix + ".points.txt")
    metadata = output.with_suffix(output.suffix + ".json")
    if any(path in (args.bc, args.points) for path in (output, query, metadata)):
        raise SystemExit("output must not overwrite the carrier or input points")
    hashes = {"bc": sha256(args.bc), "points": sha256(args.points),
              "probe": sha256(args.probe)}
    text = args.points.read_text()
    if "," in text.splitlines()[0]:
        rows = csv.DictReader(text.splitlines())
        triples = [(row["s_tor"], row["phi_B_rad"], row["theta_B_rad"])
                   for row in rows]
        query.write_text("".join(" ".join(row) + "\n" for row in triples))
    else:
        query.write_text(text)
    command = [str(args.probe), str(args.bc), str(args.inp_swi),
               str(query), str(output)]
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    subprocess.run(command, check=True, env=environment)
    if hashes["bc"] != sha256(args.bc) or hashes["points"] != sha256(args.points):
        raise SystemExit("carrier or input points changed during field readback")
    metadata.write_text(json.dumps({
        "command": command, "bc_path": str(args.bc), "bc_sha256": hashes["bc"],
        "points_path": str(args.points), "points_sha256": hashes["points"],
        "probe_sha256": hashes["probe"], "query_sha256": sha256(query),
        "output_sha256": sha256(output), "inp_swi": args.inp_swi, "bfac": 1,
        "coordinates": ["s_tor", "phi_B_rad", "theta_B_rad"],
        "units": {"Bmod": "G", "Bphcov": "G cm", "Bthcov": "G cm",
                  "iota": "dimensionless", "sqrtg": "cm^3", "psi_pr": "G cm^2"},
        "scope": "axisymmetric field only; no profile, orbit or torque initialization",
    }, indent=2) + "\n")
    print(metadata)


if __name__ == "__main__":
    main()
