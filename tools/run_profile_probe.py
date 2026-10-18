#!/usr/bin/env python3
"""Run the native profile initializer with frozen carrier and profile identities."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('probe', 'bc', 'plasma', 'mach', 'points', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--inp-swi', type=int, choices=(8, 9), default=9)
    parser.add_argument('--efac', type=float, required=True)
    parser.add_argument('--bfac', type=float, default=1.)
    args = parser.parse_args()
    inputs = {name: getattr(args, name).resolve(strict=True)
              for name in ('probe', 'bc', 'plasma', 'mach', 'points')}
    output = args.output.resolve()
    query = output.with_suffix(output.suffix + '.points.txt')
    sidecar = output.with_suffix(output.suffix + '.json')
    for path in (output, query, sidecar):
        for source in inputs.values():
            if path == source or (path.exists() and path.samefile(source)):
                raise SystemExit('output must not overwrite an input or the probe')
    output.parent.mkdir(parents=True, exist_ok=True)
    sidecar.unlink(missing_ok=True)
    hashes = {name: sha256(path) for name, path in inputs.items()}
    text = inputs['points'].read_text()
    if ',' in text.splitlines()[0]:
        rows = csv.DictReader(text.splitlines())
        query.write_text(''.join(row['s_tor'] + '\n' for row in rows))
    else:
        query.write_text(text)
    command = [str(inputs['probe']), str(inputs['bc']), str(args.inp_swi),
               str(inputs['plasma']), str(inputs['mach']), str(args.efac),
               str(args.bfac), str(query), str(output)]
    subprocess.run(command, check=True,
                   env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    if any(sha256(path) != hashes[name] for name, path in inputs.items()):
        raise SystemExit('an input or the probe changed during readback')
    metadata = {name + '_path': str(path) for name, path in inputs.items()}
    metadata.update({name + '_sha256': value for name, value in hashes.items()})
    metadata.update(command=command, query_sha256=sha256(query),
                    output_sha256=sha256(output), inp_swi=args.inp_swi,
                    efac=args.efac, bfac=args.bfac, coordinates=['s_tor'],
                    scope='native profile and collision initialization only; no orbits or torque',
                    scientific_admission=False)
    sidecar.write_text(json.dumps(metadata, indent=2) + '\n')
    print(sidecar)


if __name__ == '__main__':
    main()
