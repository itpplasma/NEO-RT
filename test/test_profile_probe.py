#!/usr/bin/env python3
"""Manufactured knot/midpoint oracle for the executed profile initializer."""
import os
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline
from test_field_probe import write_fixture

MU, EV, QE = 1.660538e-24, 1.602176e-12, 4.803204e-10
PMASS, EMASS, CEV, CQE = 1.6726e-24, 9.1094e-28, 1.6022e-12, 4.8032e-10
M1, M2, Z1, Z2 = 2.014539865995238, 12., 1., 6.
EFAC, BFAC, R0 = -2.5, .8, 500.


def run(probe, directory, x, plasma, rotation, query, invalid=False):
    bc, p, m, q, out = [directory / n for n in ('field.bc', 'plasma.in', 'profile.in', 'points.in', 'profiles.csv')]
    write_fixture(bc, -1, 1, 1)
    with p.open('w') as stream:
        stream.write(f'# manufactured background\n{len(x)} {M1:.17e} {M2} {Z1} {Z2}\n# s n1 n2 T1 T2 Te\n')
        np.savetxt(stream, np.column_stack((x, plasma)), fmt='%.17e')
    np.savetxt(m, rotation, fmt='%.17e')
    np.savetxt(q, query, fmt='%.17e')
    command = [str(probe), str(bc), '9', str(p), str(m), str(EFAC), str(BFAC), str(q), str(out)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=30,
                            env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    if invalid:
        assert result.returncode != 0, 'unsupported or invalid input accepted'
        return
    assert result.returncode == 0, result.stdout + result.stderr
    return np.genfromtxt(out, delimiter=',', names=True, comments='#', skip_header=sum(
        line.startswith('#') for line in out.read_text().splitlines()))


def reference(x, plasma, rotation, query):
    values = []
    for table, col in [(np.column_stack((x, plasma)), k) for k in range(1, 6)] + [(rotation, 1)]:
        spline = CubicSpline(table[:, 0], table[:, col], bc_type='natural')
        values.extend((spline(query), spline(query, 1)))
    n1, dn1, n2, dn2, t1, dt1, t2, dt2, te, dte, mach, dmach = values
    mach, dmach = mach * EFAC / BFAC, dmach * EFAC / BFAC
    vth = np.sqrt(2 * t1 * EV / (M1 * MU))
    dvth = vth * dt1 / (2 * t1)
    ne = Z1 * n1 + Z2 * n2
    zeff = (Z1**2 * n1 + Z2**2 * n2) / ne
    pressure = (ne * te + n1 * t1 + n2 * t2) * EV
    base = [query, n1, dn1, n2, dn2, t1, dt1, t2, dt2, te, dte, mach, dmach,
            vth, dvth, vth * mach / R0, (dvth * mach + vth * dmach) / R0]
    constant = lambda value: np.full_like(query, value)
    base.extend(map(constant, (M1 * MU, Z1 * QE)))
    base.extend((ne, zeff, pressure))
    base.extend(map(constant, (M1 * PMASS, M2 * PMASS, EMASS)))
    # Independently assemble the legacy collision normalization, keeping its
    # own mass/energy constants separate from orbit normalization.
    energy = M1 * PMASS * vth**2 / (2 * EV)
    v0 = np.sqrt(2 * energy * CEV / (M1 * PMASS))
    vel = [v0 / np.sqrt(2 * t * CEV / mass) for t, mass in
           ((t1, M1 * PMASS), (t2, M2 * PMASS), (te, EMASS))]
    lambdas = [23 - np.log(np.maximum(np.finfo(float).eps,
               np.sqrt(n * z**2 / t) * Z1 * z * (M1 + mass) / (M1 * t + mass * energy)))
               for n, z, t, mass in ((n1, Z1, t1, M1), (n2, Z2, t2, M2))]
    lambdas.append(24 - np.log(np.sqrt(ne) / te))
    factor = 2 * 3.14159265358979 * ne * CQE**4 * Z1**2 / ((M1 * PMASS)**2 * v0**4)
    base.extend((factor * Z1**2 * lambdas[0] * n1 / ne * vel[0],
                 factor * Z2**2 * lambdas[1] * n2 / ne * vel[1], factor * lambdas[2] * vel[2]))
    base.extend(vel)
    base.extend(energy / t for t in (t1, t2, te))
    base.append(constant(R0))
    return np.column_stack(base)


def main():
    probe = Path(sys.argv[1]).resolve()
    x = np.array([.1, .13, .24, .43, .68, .9])
    xm = np.array([.1, .16, .29, .54, .79, .9])
    query = np.unique(np.r_[x, xm, (x[1:] + x[:-1]) / 2, (xm[1:] + xm[:-1]) / 2])[::-1]
    rotation = np.column_stack((xm, .08 - .03 * xm + .02 * xm**2))
    with tempfile.TemporaryDirectory(prefix='neort-native-profile-') as tmp:
        directory = Path(tmp)
        for carbon in (False, True):
            plasma = np.column_stack((2e13 * (1 + .2 * x),
                1e12 * (1 + .4 * x + .7 * x**2) if carbon else np.zeros_like(x),
                900 + 800 * (1 - x**2), 1100 + 200 * np.sin(x), 700 + 300 * np.exp(-x)))
            output = run(probe, directory, x, plasma, rotation, query)
            actual = np.column_stack([output[name] for name in output.dtype.names])
            expected = reference(x, plasma, rotation, query)
            scale = np.maximum(np.max(abs(expected), axis=0), np.finfo(float).tiny)
            np.testing.assert_allclose((actual - expected) / scale, 0., rtol=0., atol=3e-11)
            assert np.any(actual[:, 17] != actual[:, 22]), 'orbit/kernel mass distinction lost'
        if len(sys.argv) == 3:
            check_wrapper(probe, Path(sys.argv[2]).resolve(), directory, query, expected)
        for q in (np.array([.099]), np.array([.901]), np.array([np.nan])):
            run(probe, directory, x, plasma, rotation, q, invalid=True)
        for bad_x in (np.r_[x[:2], x[1], x[3:]], x[::-1]):
            run(probe, directory, bad_x, plasma, rotation, query, invalid=True)
        run(probe, directory, x, plasma, rotation[1:-1], query, invalid=True)
    print('Native eD/eDC profiles, derivatives, mass conventions and collision arrays: PASS')


def check_wrapper(probe, wrapper, directory, query, expected):
    points, output = directory / 'points.csv', directory / 'wrapped.csv'
    points.write_text('s_tor,phi_B_rad,theta_B_rad\n' + ''.join(
        f'{s:.17e},0,0\n' for s in query))
    paths = dict(probe=probe, bc=directory / 'field.bc', plasma=directory / 'plasma.in',
                 mach=directory / 'profile.in', points=points)
    command = [sys.executable, str(wrapper), '--output', str(output),
               '--efac', str(EFAC), '--bfac', str(BFAC)]
    command.extend(value for key, path in paths.items() for value in ('--' + key, str(path)))
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    sidecar = json.loads(output.with_suffix('.csv.json').read_text())
    paths.update(output=output, query=output.with_suffix('.csv.points.txt'))
    for name, path in paths.items():
        assert sidecar[name + '_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert sidecar['efac'] == EFAC and sidecar['bfac'] == BFAC
    assert sidecar['scientific_admission'] is False
    values = np.genfromtxt(output, delimiter=',', names=True, comments='#',
                          skip_header=sum(line.startswith('#') for line in output.read_text().splitlines()))
    actual = np.column_stack([values[name] for name in values.dtype.names])
    scale = np.maximum(np.max(abs(expected), axis=0), np.finfo(float).tiny)
    np.testing.assert_allclose((actual - expected) / scale, 0., rtol=0., atol=3e-11)
    # Sidecar aliases must be caught before the carrier can be overwritten.
    protected = directory / 'protected.csv'
    protected.with_suffix('.csv.json').hardlink_to(paths['bc'])
    old_hash = hashlib.sha256(paths['bc'].read_bytes()).hexdigest()
    command[command.index('--output') + 1] = str(protected)
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert hashlib.sha256(paths['bc'].read_bytes()).hexdigest() == old_hash


if __name__ == '__main__':
    main()
