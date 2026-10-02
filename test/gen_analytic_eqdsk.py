#!/usr/bin/env python3
"""Write the analytic GEQDSK fixtures used by the direct-EQDSK tests.

    gen_analytic_eqdsk.py circ|solovev|solovev_bt_reversed OUTPUT.geqdsk

Both fields have a constant F = R0*B0 and a purely toroidal current, so they
admit Boozer coordinates (libneo's efit_to_boozer.x turns them into the
Boozer files, see gen_analytic_boozer.py).

circ: concentric circles about (R0, 0) with

    psi(rho) = B0/(2*c) * log(1 + c*rho**2/Q0),  c = (QA - Q0)/a**2,

so psi'(rho) = B0*rho/qm(rho), qm = Q0 + c*rho**2.  Closed forms, none of
which goes through a GEQDSK reader:

    |B| = sqrt(F**2 + psi'(rho)**2)/R,  R = R0 + rho*cos(theta)
    local pitch  B^phi/B^theta = qm(rho)*R0/R
    safety factor q(rho) = qm(rho)/sqrt(1 - (rho/R0)**2)
    toroidal flux per radian  F*(R0 - sqrt(R0**2 - rho**2))

On this field |B| is constant times 1/R on every surface, so the Boozer
toroidal angle equals the cylindrical one.

solovev: psi = cs*((R**2 - R0**2)**2/(4*R0**2) + R**2*Z**2/KAPPA**2), whose
Grad-Shafranov source is proportional to R**2 (constant p', F' = 0).  The
surfaces are shaped and |B| is not proportional to 1/R on them, so the
Boozer stream function is nonzero; this is the cross-backend fixture.
"""
import sys

import numpy as np

R0 = 1.65  # magnetic axis [m]
A = 0.50  # outboard midplane minor radius of the boundary [m]
B0 = 2.0  # toroidal field on axis [T]
Q0 = 1.2  # circ: qm on axis
QA = 3.5  # circ: qm at rho = a
ELONG = 1.3  # solovev: elongation of the axis ellipse
KAPPA = ELONG * R0  # solovev: Z scale of psi [m]
NR = 65
NZ = 65
C = (QA - Q0) / A**2


def psi_circ(R, Z):
    rho = np.hypot(R - R0, Z)
    return B0 / (2.0 * C) * np.log1p(C * rho**2 / Q0)


def q_circ(psi):
    rho = np.sqrt(Q0 * np.expm1(psi * 2.0 * C / B0) / C)
    return (Q0 + C * rho**2) / np.sqrt(1.0 - (rho / R0) ** 2)


def solovev_scale():
    """Scale cs so that q on the axis is Q0."""
    # Near the axis psi ~ cs*(x**2 + Z**2/ELONG**2), x = R - R0, which gives
    # toroidal flux pi*B0*ELONG*psi/cs per unit psi and q(axis) = B0*ELONG/(2*cs).
    return B0 * ELONG / (2.0 * Q0)


def psi_solovev(R, Z):
    cs = solovev_scale()
    return cs * ((R**2 - R0**2) ** 2 / (4.0 * R0**2) + R**2 * Z**2 / KAPPA**2)


def boundary_contour(psi_fn, psi_b, n=129):
    """Points of psi = psi_b on rays from the axis (bisection)."""
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    lo = np.zeros(n)
    hi = np.full(n, 1.5 * A)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        inside = psi_fn(R0 + mid * np.cos(theta), mid * np.sin(theta)) < psi_b
        lo = np.where(inside, mid, lo)
        hi = np.where(inside, hi, mid)
    return np.column_stack([R0 + lo * np.cos(theta), lo * np.sin(theta)])


def records(values):
    """Format a flat sequence in the fixed GEQDSK 5e16.9 layout."""
    values = np.asarray(values, dtype=float).ravel()
    lines = []
    for start in range(0, values.size, 5):
        lines.append("".join(f"{v:16.9e}" for v in values[start:start + 5]))
    return "\n".join(lines) + "\n"


def main(kind, path):
    toroidal_sign = -1.0 if kind.endswith("_bt_reversed") else 1.0
    psi_fn = {"circ": psi_circ, "solovev": psi_solovev}[
        kind.removesuffix("_bt_reversed")]
    rdim = 3.2 * A
    # The Solovev LCFS reaches |Z|=0.808 m; ±0.8 m clips edge-flux contours.
    zdim = (3.6 if kind.startswith("solovev") else 3.2) * A
    rleft = R0 - 0.5 * rdim
    R = rleft + np.linspace(0.0, rdim, NR)
    Z = -0.5 * zdim + np.linspace(0.0, zdim, NZ)
    RR, ZZ = np.meshgrid(R, Z)  # (NZ, NR): psirz is written R fastest
    psirz = psi_fn(RR, ZZ)

    simag = 0.0
    sibry = float(psi_fn(R0 + A, 0.0))
    psi_grid = np.linspace(simag, sibry, NR)
    # qpsi is informational only: libneo recomputes q from the field.
    qpsi = q_circ(psi_grid) if kind == "circ" else np.full(NR, Q0)
    fpol = np.full(NR, toroidal_sign * R0 * B0)
    zeros = np.zeros(NR)
    current = 1.0e6

    lcfs = boundary_contour(psi_fn, sibry)
    theta = np.linspace(0.0, 2.0 * np.pi, 129)
    lim = np.column_stack([R0 + 1.5 * A * np.cos(theta), 1.5 * A * np.sin(theta)])

    with open(path, "w") as out:
        out.write(f"{'NEO-RT analytic ' + kind:48s}{0:4d}{NR:4d}{NZ:4d}\n")
        out.write(records([rdim, zdim, R0, rleft, 0.0]))
        out.write(records([R0, 0.0, simag, sibry, toroidal_sign * B0]))
        out.write(records([current, simag, 0.0, R0, 0.0]))
        out.write(records([0.0, 0.0, sibry, 0.0, 0.0]))
        out.write(records(fpol))
        out.write(records(zeros))
        out.write(records(zeros))
        out.write(records(zeros))
        out.write(records(psirz))
        out.write(records(qpsi))
        out.write(f"{lcfs.shape[0]:5d}{lim.shape[0]:5d}\n")
        out.write(records(lcfs))
        out.write(records(lim))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
