#!/usr/bin/env python3
"""Generate the public-safe circular tokamak inputs for the POTATO resonance gate.

The resonance regression that motivates this test (a dummy-procedure trampoline
corrupting the threadprivate get_rescond by-products under OpenMP) is
equilibrium-INDEPENDENT: it corrupts root weights regardless of the field. So a
synthetic concentric-circle equilibrium reproduces it without shipping any
experimental AUG #30835 g-file or kinetic profiles, which must not enter public
NEO-RT.

The equilibrium is a low-beta circular tokamak with concentric circular flux
surfaces

    psi(R,Z) = psi_edge * ((R-R0)^2 + Z^2) / a^2 ,   psi(axis)=0

so the poloidal field is B_pol = |grad psi| / R, growing linearly with minor
radius like a uniform-current screw pinch. The toroidal field is a constant
f = R0*B0 (vacuum, low beta), pressure is a small parabola. Everything is in
SI (m, T, Vs); POTATO's libneo reader converts to its internal cgs grid.

Outputs (written next to this script):
  circ.eqdsk       g-file consumed as the gfile by field_divB0.inp
  bmod_n.dat       unformatted n-harmonic |B| perturbation on the eqdsk R-Z grid
                   (the resonance driver; a smooth nonzero pattern so resonant
                    points carry weight)
  convexwall.dat   circular stretch-coords wall (cm) enclosing the boundary
  profile_poly.in  monotonic density/temperature/potential polynomials
The static potato.in and field_divB0.inp are committed alongside.
"""
import os
import struct
import sys

import numpy as np

# Round-trip helper lives in the libneo python package fetched by the build.
THIS = os.path.dirname(os.path.abspath(__file__))
LIBNEO_PY = os.path.join(
    THIS, "..", "..", "build", "_deps", "libneo-src", "python"
)
if os.path.isdir(LIBNEO_PY):
    sys.path.insert(0, LIBNEO_PY)
from libneo.eqdsk_base import read_eqdsk, write_eqdsk  # noqa: E402

# --- circular equilibrium parameters (public-safe synthetic) -----------------
R0 = 1.60        # magnetic axis major radius [m]
A = 0.50         # minor radius of the boundary flux surface [m]
B0 = 2.0         # toroidal field on axis [T]
P0 = 5.0e3       # on-axis pressure [Pa]; low beta
# Monotonically rising (sheared) safety factor q(r) = Q0 + (QA-Q0)*(r/a)^2.
# Magnetic shear is essential: a constant-q (degenerate) equilibrium makes every
# poloidal harmonic resonate at the same frequency ratio across all J_perp, so
# the resonance root finder hits its oscillatory-root cap on nearly every mode
# and the run neither finishes fast nor produces isolated resonance lines. A
# rising q with q=1..few crossed over the radius gives well-separated lines.
Q0 = 1.5         # safety factor on axis
QA = 4.0         # safety factor at the boundary r=a
NR = 65
NZ = 65


def q_of_r(r):
    """Sheared safety factor as a function of minor radius r in [0, A]."""
    return Q0 + (QA - Q0) * (r / A) ** 2


def psi_of_r():
    """Poloidal flux psi(r) from the prescribed q, psi(0)=0.

    On the outboard midplane B_tor = R0*B0/R and B_pol = psi'(r)/R, so
    q = r*B_tor/(R*B_pol) = r*R0*B0/psi'(r), giving psi'(r) = R0*B0*r/q(r).
    Integrate outward; the same psi(rho) is then laid on concentric circles
    rho = sqrt((R-R0)^2 + Z^2) (geometric-axis model equilibrium).
    """
    n = 2048
    r = np.linspace(0.0, A, n)
    dpsidr = np.zeros(n)
    dpsidr[1:] = R0 * B0 * r[1:] / q_of_r(r[1:])
    psi = np.concatenate([[0.0], np.cumsum(0.5 * (dpsidr[1:] + dpsidr[:-1])
                                          * np.diff(r))])
    return r, psi


def build_eqdsk():
    # Box generously around the boundary so orbits stay inside the grid.
    rboxleft = R0 - 1.6 * A
    rboxlength = 3.2 * A
    zboxmid = 0.0
    zboxlength = 3.2 * A

    R = rboxleft + np.linspace(0.0, rboxlength, NR)
    Z = zboxmid - 0.5 * zboxlength + np.linspace(0.0, zboxlength, NZ)
    RR, ZZ = np.meshgrid(R, Z)  # shape (NZ, NR), matching read/write layout

    # psi(r) from the sheared q; psi_edge = psi(a). Concentric circular flux
    # surfaces: interpolate psi(rho) on rho = sqrt((R-R0)^2 + Z^2), clamped to
    # the boundary value outside r=a so the SOL stays monotone.
    r_tab, psi_tab = psi_of_r()
    psi_edge = psi_tab[-1]
    rho = np.sqrt((RR - R0) ** 2 + ZZ ** 2)
    PsiVs = np.interp(rho, r_tab, psi_tab, right=psi_edge)

    # Toroidal field: constant f = R0*B0 (vacuum, low beta), so B_tor = f/R.
    fprof = np.full(NR, R0 * B0)
    fdfdpsiprof = np.zeros(NR)  # f constant -> f f' = 0

    # Small parabolic pressure in normalized flux s=psi/psi_edge, zero at edge.
    s = np.linspace(0.0, 1.0, NR)
    ptotprof = P0 * (1.0 - s)
    dpressdpsiprof = np.full(NR, -P0 / psi_edge)

    # q profile on the uniform s grid (s = psi/psi_edge -> r via the inverse).
    r_on_s = np.interp(s * psi_edge, psi_tab, r_tab)
    qprof = q_of_r(r_on_s)

    # Plasma current from Ampere on the boundary loop: Ip = B_pol(a)*2*pi*a/mu0,
    # with B_pol(a) = psi'(a)/(R0+a) = R0*B0*a/(q(a)*(R0+a)).
    mu0 = 4.0e-7 * np.pi
    Bpol_edge = R0 * B0 * A / (QA * (R0 + A))
    Ip = Bpol_edge * 2.0 * np.pi * A / mu0

    # Circular boundary and a slightly larger circular limiter.
    theta = np.linspace(0.0, 2.0 * np.pi, 129)
    lcfs = np.column_stack([R0 + A * np.cos(theta), A * np.sin(theta)])
    lim = np.column_stack(
        [R0 + 1.1 * A * np.cos(theta), 1.1 * A * np.sin(theta)]
    )

    eqdata = {
        "header": "POTATO circular synthetic gate  CIRC 00000        ",
        "nrgr": NR,
        "nzgr": NZ,
        "rboxlength": rboxlength,
        "zboxlength": zboxlength,
        "R0": R0,
        "rboxleft": rboxleft,
        "zboxmid": zboxmid,
        "Rpsi0": R0,
        "Zpsi0": 0.0,
        "PsiaxisVs": 0.0,
        "PsiedgeVs": psi_edge,
        "Btor_at_R0": B0,
        "Ip": Ip,
        "fprof": fprof,
        "ptotprof": ptotprof,
        "fdfdpsiprof": fdfdpsiprof,
        "dpressdpsiprof": dpressdpsiprof,
        "PsiVs": PsiVs,
        "qprof": qprof,
        "npbound": lcfs.shape[0],
        "nplimiter": lim.shape[0],
        "Lcfs": lcfs,
        "Limiter": lim,
    }
    return eqdata, R, Z


def write_bmod_n(path, R, Z):
    """Write the synthetic n-harmonic |B| perturbation POTATO's pertham reads.

    Format (Fortran unformatted sequential, default record markers): nrad,nzet;
    then rad,zet (cm, matching the field grid); then the real and imaginary
    amplitude arrays on the (nrad,nzet) grid. The pattern is a smooth poloidal
    m-like ripple peaked near the edge so resonant points carry nonzero weight;
    the magnitude is arbitrary (the gate compares resonance-line counts, not
    torque values).
    """
    rad_cm = R * 100.0
    zet_cm = Z * 100.0
    nrad = rad_cm.size
    nzet = zet_cm.size

    RR, ZZ = np.meshgrid(rad_cm, zet_cm, indexing="ij")  # (nrad,nzet)
    r0_cm = R0 * 100.0
    a_cm = A * 100.0
    rho = np.sqrt((RR - r0_cm) ** 2 + ZZ ** 2)
    tht = np.arctan2(ZZ, RR - r0_cm)
    env = (rho / a_cm) ** 2  # grows toward the edge
    bre = (1.0e3 * env * np.cos(3.0 * tht)).astype(np.float64)
    bim = (1.0e3 * env * np.sin(3.0 * tht)).astype(np.float64)

    def rec(payload):
        n = len(payload)
        return struct.pack("=i", n) + payload + struct.pack("=i", n)

    with open(path, "wb") as f:
        f.write(rec(struct.pack("=ii", nrad, nzet)))
        f.write(rec(rad_cm.tobytes() + zet_cm.tobytes()))
        # Fortran column-major: write arrays transposed so on-disk order is
        # (i fastest). numpy default is C-order, so flatten in Fortran order.
        f.write(rec(bre.tobytes(order="F") + bim.tobytes(order="F")))


def write_convexwall(path):
    """Circular stretch-coords wall (cm) enclosing the boundary with margin.

    stretch_coords reads (R,Z) pairs in the field grid units (cm) and pulls
    points outside this convex wall back in. A circle of radius 1.12*a around
    the axis sits just outside the a=A boundary so traced orbits are not
    distorted but far-SOL excursions are bounded.
    """
    rw = 1.12 * A * 100.0
    r0_cm = R0 * 100.0
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    with open(path, "w") as f:
        for t in theta:
            f.write(f"{r0_cm + rw*np.cos(t):24.16e}{rw*np.sin(t):24.16e}\n")


def write_profile_poly(path):
    """Monotonic species profiles as descending-power polynomials in s_pol.

    POTATO reads ten coefficients per array (index 0..9, descending powers of
    the normalized poloidal flux s_pol in [0,1]); evaluation is
    sum_k poly(k) * s_pol**(9-k). Linear profiles padded with leading zeros:
      density     n = 5e13 * (1 - 0.5 s_pol)  [cm^-3]
      temperature T = 2000  * (1 - 0.7 s_pol)  [eV]
      potential   Phi = 300 * (1 - s_pol)      [V]
    The reader takes line 1 = density, line 2 = dummy, line 3 = temperature,
    line 4 = potential; lines 1-2 of the file are comment/dummy header.
    """
    def lin(slope, intercept):
        return [0.0] * 8 + [slope, intercept]

    rows = [
        lin(-2.5e13, 5.0e13),   # density
        [0.0] * 10,             # dummy (reader discards)
        lin(-1.4e3, 2.0e3),     # temperature
        lin(-3.0e2, 3.0e2),     # potential
    ]
    with open(path, "w") as f:
        f.write("% Public-safe synthetic circular-case profiles "
                "(NOT AUG #30835)\n")
        f.write("% ten coefficients per line, descending powers of s_pol; "
                "order: density, dummy, temperature, potential\n")
        for row in rows:
            f.write(" ".join(f"{c:.16e}" for c in row) + "\n")


def main():
    eqdata, R, Z = build_eqdsk()
    eqdsk_path = os.path.join(THIS, "circ.eqdsk")
    write_eqdsk(eqdsk_path, eqdata)
    # libneo's Fortran GEQDSK reader follows the conventional fixed-width
    # ``2i5`` boundary-count record.  The Python writer emits those two integers
    # list-directed, so normalize that one record for cross-reader portability.
    with open(eqdsk_path) as stream:
        lines = stream.readlines()
    count_record = f"{eqdata['npbound']} {eqdata['nplimiter']}"
    matches = [i for i, line in enumerate(lines) if line.strip() == count_record]
    assert len(matches) == 1
    lines[matches[0]] = f"{eqdata['npbound']:5d}{eqdata['nplimiter']:5d}\n"
    with open(eqdsk_path, "w") as stream:
        stream.writelines(lines)

    psi_edge = eqdata["PsiedgeVs"]

    # Round-trip check: read back and confirm the structural fields survive.
    back = read_eqdsk(eqdsk_path)
    assert back["nrgr"] == NR and back["nzgr"] == NZ
    assert abs(back["R0"] - R0) < 1e-6
    assert abs(back["Btor_at_R0"] - B0) < 1e-6
    assert abs(back["PsiedgeVs"] - psi_edge) < 1e-6
    assert np.allclose(back["fprof"], R0 * B0, atol=1e-6)
    assert np.allclose(back["PsiVs"], eqdata["PsiVs"], atol=1e-6 * psi_edge)

    write_bmod_n(os.path.join(THIS, "bmod_n.dat"), R, Z)
    write_convexwall(os.path.join(THIS, "convexwall.dat"))
    write_profile_poly(os.path.join(THIS, "profile_poly.in"))

    print("wrote circ.eqdsk, bmod_n.dat, convexwall.dat, profile_poly.in")
    print(f"  R0={R0} a={A} B0={B0} psi_edge={psi_edge:.4f} "
          f"q_axis={eqdata['qprof'][0]:.3f} q_edge={eqdata['qprof'][-1]:.3f} "
          f"Ip={eqdata['Ip']:.3e} A")


if __name__ == "__main__":
    main()
