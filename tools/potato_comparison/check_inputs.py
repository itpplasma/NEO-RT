"""Check native field exports against Cartesian oracle samples and profiles.

The oracle holds physical samples of the analytic frozen field, not solver
outputs. Complex amplitudes are peak phasors: B = Re[A exp(i n phi)].
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline


def read_records(path):
    records = []
    with path.open("rb") as f:
        while marker := f.read(4):
            size = struct.unpack("=i", marker)[0]
            data = f.read(size)
            if len(data) != size or f.read(4) != marker:
                raise ValueError("invalid Fortran record")
            records.append(data)
    return records


def read_bmod(path):
    dims, coords, fields = read_records(path)
    nr, nz = struct.unpack("=ii", dims)
    coords, fields = np.frombuffer(coords, "=f8"), np.frombuffer(fields, "=f8")
    amplitude = (fields[:nr*nz] + 1j * fields[nr*nz:]).reshape((nr, nz), order="F")
    return coords[:nr] / 100, coords[nr:] / 100, amplitude / 1e4


def read_bc(path):
    lines = path.read_text().splitlines()
    nm, nn, ns = map(int, lines[5].split()[:3])
    count, position, surfaces = (nm + 1) * (nn + 1), 6, []
    for _ in range(ns):
        position += 2
        header = np.array(list(map(float, lines[position].split())))
        position += 2
        rows = np.array([list(map(float, line.split()))
                         for line in lines[position:position + count]])
        position += count
        surfaces.append((header, rows))
    return surfaces


def fourier_real(rows, theta, phi, column):
    phase = theta[:, None] * rows[:, 0] + phi[:, None] * rows[:, 1]
    return np.sum(rows[:, column] * np.cos(phase) +
                  rows[:, column + 1] * np.sin(phase), axis=1)


def error_metric(actual, expected):
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise ValueError("nonfinite physical samples")
    return float(np.max(np.abs(actual - expected)) / np.max(np.abs(expected)))


def read_equilibrium(path):
    lines = path.read_text().splitlines()
    nr, nz = map(int, lines[0].split()[-2:])
    count = 20 + 5 * nr + nr * nz
    values = []
    for line in lines[1:]:
        values.extend(float(line[i:i + 16]) for i in range(0, len(line), 16))
        if len(values) >= count:
            break
    values = np.array(values)
    rd, zd, _, rl, zm = values[:5]
    r, z = np.linspace(rl, rl + rd, nr), np.linspace(zm - zd / 2, zm + zd / 2, nz)
    start = 20 + 4 * nr
    psi = values[start:start + nr * nz].reshape((nr, nz), order="F")
    return (RectBivariateSpline(r, z, psi, kx=5, ky=5),
            values[20:20 + nr], values[start + nr * nz:count], values[7:9])


def audit(root):
    root = Path(root)
    common = root / "inputs/common"
    oracle = np.load(common / "field_oracle.npz")
    eq, pert = read_bc(common / "eq.bc"), read_bc(common / "pert.bc")
    r, z, amp = read_bmod(common / "bmod_n.dat")
    re = RectBivariateSpline(r, z, amp.real, kx=5, ky=5)
    im = RectBivariateSpline(r, z, amp.imag, kx=5, ky=5)
    psi, fpol, qpol, psi_limits = read_equilibrium(common / "efit.eqdsk")
    errors = {}
    for index, s in enumerate(oracle["s"]):
        eqrow = min(eq, key=lambda surface: abs(surface[0][0] - s))
        prow = min(pert, key=lambda surface: abs(surface[0][0] - s))
        if abs(eqrow[0][0] - s) > 1e-9 or abs(prow[0][0] - s) > 1e-9:
            raise ValueError("oracle sample must coincide with an exported surface")
        theta, phi = oracle["theta_b"][index], oracle["phi_b"][index]
        expected = oracle["deltaB"][index]
        bn = fourier_real(prow[1], theta, phi, 8)
        rp, zp = oracle["R"][index], oracle["Z"][index]
        ap = re.ev(rp, zp) + 1j * im.ev(rp, zp)
        bp = np.real(ap * np.exp(-3j * oracle["phi"][index]))
        name = f"s={s:.3f}"
        errors[name] = {
            "NEORT_deltaB": error_metric(bn, expected),
            "POTATO_deltaB": error_metric(bp, expected),
            "NEORT_B0": error_metric(fourier_real(eqrow[1], theta, phi, 8), oracle["B0"][index]),
            "NEORT_R": error_metric(fourier_real(eqrow[1], theta, phi, 2), rp),
            "NEORT_Z": error_metric(fourier_real(eqrow[1], theta, phi, 4), zp),
        }
        br = -psi.ev(rp, zp, dy=1) / rp
        bz = psi.ev(rp, zp, dx=1) / rp
        b0p = np.sqrt(br**2 + bz**2 + (fpol[0] / rp)**2)
        errors[name]["POTATO_B0"] = error_metric(b0p, oracle["B0"][index])
        spol = (psi.ev(rp, zp) - psi_limits[0]) / (psi_limits[1] - psi_limits[0])
        errors[name]["s_pol_equals_s_tor"] = float(np.max(abs(spol - s)))
        errors[name]["q_relative_error"] = float(np.max(abs(qpol / (13 / 6) - 1)))
        phasor = np.sum((prow[1][:, 8] - 1j * prow[1][:, 9]) *
                       np.exp(1j * theta[:, None] * prow[1][:, 0]), axis=1)
        errors[name]["relative_phasor"] = error_metric(phasor, np.conj(ap) *
                                                 np.exp(-3j * oracle["nu"][index]))
    # Direct phi quadrature: the two-sided coefficient is half the peak phasor.
    ph = np.arange(128) * 2 * np.pi / 128
    a = 2.3 + 0.7j
    physical = np.real(a * np.exp(-3j * ph))
    fourier = np.mean(physical * np.exp(3j * ph))
    errors["two_sided_coefficient_over_peak"] = float(abs(fourier / a))
    if abs(fourier - a / 2) > 1e-13:
        raise ValueError("independent toroidal Fourier quadrature fails")
    # POTATO normalized H/Eref divided by E/Eref is NEO-RT's H/E.
    b, eta, e, bn = 4.8, .13, .73, 2.4e-4 + 1.2e-4j
    hp = (2 * e / b - eta * e) * bn
    hn = (2 - eta * b) * bn / b
    errors["normalized_hamiltonian_relative_error"] = float(abs(hp / e - hn) / abs(hn))
    for values in errors.values():
        if isinstance(values, dict) and max(values.values()) > 1e-3:
            raise ValueError(f"field export mismatch: {values}")
    profile_error = 0.
    for c in json.loads((root / "inputs/scan.json").read_text())["cases"]:
        d = root / "inputs" / c["id"]
        if c["kind"] == "half_amplitude":
            _, _, half = read_bmod(d / "bmod_n.dat")
            errors[c["id"] + "_POTATO_phasor"] = error_metric(half, amp / 2)
            control = read_bc(d / "pert.bc")
            for base, half_surface in zip(pert, control):
                theta = np.linspace(0, 2 * np.pi, 64, endpoint=False)
                phi = np.full_like(theta, .37)
                expected = fourier_real(base[1], theta, phi, 8) / 2
                actual = fourier_real(half_surface[1], theta, phi, 8)
                if error_metric(actual, expected) > 1e-12:
                    raise ValueError("half-amplitude Boozer control mismatch")
            if errors[c["id"] + "_POTATO_phasor"] > 1e-12:
                raise ValueError("half-amplitude cylindrical control mismatch")
        poly = np.loadtxt(d / "profile_poly.in", skiprows=2)
        table = np.loadtxt(d / "plasma.in", skiprows=3)
        velocity = np.loadtxt(d / "profile.in")
        mass = float((d / "plasma.in").read_text().splitlines()[1].split()[1]) * 1.660538e-24
        expected = c["T0_eV"] * (1 - .5 * (table[:, 0] - c["s"]))
        errors_here = [error_metric(np.polyval(poly[2], table[:, 0]), expected),
                       error_metric(table[:, 3], expected),
                       error_metric(velocity[:, 2]**2 * mass / (2 * 1.602176e-12), expected),
                       error_metric(np.polyval(poly[0], table[:, 0]), np.full(len(table), 1e10)),
                       error_metric(table[:, 1], np.full(len(table), 1e10))]
        profile_error = max(profile_error, *errors_here)
    errors["profiles_relative_error"] = profile_error
    if profile_error > 1e-12:
        raise ValueError("profile mismatch")
    print(json.dumps(errors, indent=2))
    return errors


if __name__ == "__main__":
    audit(sys.argv[1])
