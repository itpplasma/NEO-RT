#!/usr/bin/env python3
"""Generate the public-safe circular tokamak inputs for the POTATO resonance gate.

The resonance regression that motivates this test (a dummy-procedure trampoline
corrupting the threadprivate get_rescond by-products under OpenMP) is
equilibrium-INDEPENDENT: it corrupts root weights regardless of the field. So a
synthetic concentric-circle equilibrium reproduces it without shipping any
experimental AUG #30835 g-file or kinetic profiles, which must not enter public
NEO-RT.

The equilibrium is a finite-aspect geometric circular orbit/interpolation
fixture, not a Grad-Shafranov force-balance equilibrium. Fortsym owns the
normalization, inverse maps, profiles, geometry, and every serialized sample.
Everything is in SI (m, T, Vs); POTATO's libneo reader converts to its
internal cgs grid.

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

from circular_flux_continuation_generated import (
    A,
    B0,
    BMOD_IM,
    BMOD_RAD_CM,
    BMOD_RE,
    BMOD_ZET_CM,
    BOUNDARY_POINTS,
    D_VOLUME_DS_TOR_AXIS,
    D_VOLUME_DS_TOR_PROFILE,
    DPRESS_DPSIPROF,
    FFPRIMEPROF,
    FPROF,
    IP,
    LCFS,
    LIMITER,
    NR,
    NZ,
    PSI_AXIS,
    PSI_EDGE,
    PSI_PROFILE,
    PSI_VS,
    PROFILE_POLY,
    QPROF,
    TOROIDAL_FLUX,
    PTOTPROF,
    Q0,
    QA,
    R0,
    R_AXIS,
    RBOX_LEFT,
    RBOX_LENGTH,
    R_GRID,
    Z_AXIS,
    ZBOX_LENGTH,
    ZBOX_MID,
    Z_GRID,
    WALL,
    WALL_POINTS,
)

# Round-trip helper lives in the libneo python package fetched by the build.
THIS = os.path.dirname(os.path.abspath(__file__))
LIBNEO_PY = os.path.join(
    THIS, "..", "..", "build", "_deps", "libneo-src", "python"
)
if os.path.isdir(LIBNEO_PY):
    sys.path.insert(0, LIBNEO_PY)
from libneo.eqdsk_base import read_eqdsk, write_eqdsk  # noqa: E402

# --- circular equilibrium parameters (public-safe synthetic) -----------------
# The generated q table is sheared so the resonance gate has isolated lines.
def build_eqdsk():
    # The Fortsym generator owns the actual grid and every physical sample.
    # This writer only arranges those generated arrays into the EQDSK record.
    rboxleft = RBOX_LEFT
    rboxlength = RBOX_LENGTH
    zboxmid = ZBOX_MID
    zboxlength = ZBOX_LENGTH

    R = R_GRID.copy()
    Z = Z_GRID.copy()

    psi_edge = PSI_EDGE
    PsiVs = PSI_VS.copy()

    # Generated toroidal-field arrays are copied into the EQDSK record.
    fprof = FPROF.copy()
    fdfdpsiprof = FFPRIMEPROF.copy()

    ptotprof = PTOTPROF.copy()
    dpressdpsiprof = DPRESS_DPSIPROF.copy()

    # q samples are Fortsym-generated values on the uniform physical grid.
    qprof = QPROF.copy()

    Ip = IP
    lcfs = LCFS.copy()
    lim = LIMITER.copy()

    eqdata = {
        "header": "POTATO circular synthetic gate  CIRC 00000        ",
        "nrgr": NR,
        "nzgr": NZ,
        "rboxlength": rboxlength,
        "zboxlength": zboxlength,
        "R0": R_AXIS,
        "rboxleft": rboxleft,
        "zboxmid": zboxmid,
        "Rpsi0": R_AXIS,
        "Zpsi0": Z_AXIS,
        "PsiaxisVs": PSI_AXIS,
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
    if PSI_VS.shape != (NZ, NR) or PSI_PROFILE.size != NR:
        raise ValueError("Fortsym circular fixture dimensions are inconsistent")
    if QPROF.size != NR or not np.allclose(qprof, QPROF, rtol=0.0, atol=0.0):
        raise ValueError("Fortsym circular q profile was altered")
    if fprof.size != NR or fdfdpsiprof.size != NR:
        raise ValueError("Fortsym circular F profiles are inconsistent")
    if ptotprof.size != NR or dpressdpsiprof.size != NR:
        raise ValueError("Fortsym circular pressure profiles are inconsistent")
    if lcfs.shape != (BOUNDARY_POINTS, 2) or lim.shape != (BOUNDARY_POINTS, 2):
        raise ValueError("Fortsym circular boundary arrays are inconsistent")
    if BMOD_RAD_CM.size != NR or BMOD_ZET_CM.size != NZ:
        raise ValueError("Fortsym perturbation coordinates are inconsistent")
    if BMOD_RE.shape != (NR, NZ) or BMOD_IM.shape != (NR, NZ):
        raise ValueError("Fortsym perturbation arrays are inconsistent")
    if WALL.shape != (WALL_POINTS, 2) or PROFILE_POLY.shape != (4, 10):
        raise ValueError("Fortsym wall/profile arrays are inconsistent")
    if not np.all(np.isfinite(BMOD_RE)) or not np.all(np.isfinite(BMOD_IM)):
        raise ValueError("Fortsym perturbation arrays contain non-finite values")
    if D_VOLUME_DS_TOR_PROFILE.size != NR:
        raise ValueError("Fortsym circular volume profile is inconsistent")
    if not np.isfinite(D_VOLUME_DS_TOR_AXIS):
        raise ValueError("invalid Fortsym circular volume Jacobian")
    if not np.all(np.isfinite(D_VOLUME_DS_TOR_PROFILE)) or np.any(
        D_VOLUME_DS_TOR_PROFILE <= 0.0
    ):
        raise ValueError("invalid Fortsym circular volume Jacobian profile")
    if not np.all(np.isfinite(PsiVs)) or not np.all(np.isfinite(qprof)):
        raise ValueError("Fortsym circular samples contain non-finite values")
    return eqdata, R, Z


def write_bmod_n(path):
    """Serialize the Fortsym-generated perturbation POTATO's pertham reads.

    Format (Fortran unformatted sequential, default record markers): generated
    dimensions, coordinates, and real/imaginary arrays on the field grid.
    """
    rad_cm = BMOD_RAD_CM.copy()
    zet_cm = BMOD_ZET_CM.copy()
    bre = BMOD_RE.copy()
    bim = BMOD_IM.copy()
    nrad = rad_cm.size
    nzet = zet_cm.size

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
    """Serialize the Fortsym-generated circular stretch-coords wall (cm).

    stretch_coords reads generated (R,Z) pairs in the field-grid units (cm).
    """
    with open(path, "w") as f:
        for r_value, z_value in WALL:
            f.write(f"{r_value:24.16e}{z_value:24.16e}\n")


def write_profile_poly(path):
    """Serialize Fortsym-generated descending-power profile coefficients.

    POTATO reads ten generated coefficients per array. The reader takes line 1
    as density, line 2 as dummy, line 3 as temperature, and line 4 as
    potential; lines 1-2 of the file are comment/dummy header.
    """
    with open(path, "w") as f:
        f.write("% Public-safe synthetic circular-case profiles "
                "(NOT AUG #30835)\n")
        f.write("% ten coefficients per line, descending powers of s_pol; "
                "order: density, dummy, temperature, potential\n")
        for row in PROFILE_POLY:
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
    assert np.allclose(back["fprof"], TOROIDAL_FLUX, atol=1e-6)
    assert np.allclose(back["PsiVs"], eqdata["PsiVs"], atol=1e-6 * psi_edge)

    write_bmod_n(os.path.join(THIS, "bmod_n.dat"))
    write_convexwall(os.path.join(THIS, "convexwall.dat"))
    write_profile_poly(os.path.join(THIS, "profile_poly.in"))

    print("wrote circ.eqdsk, bmod_n.dat, convexwall.dat, profile_poly.in")
    print(f"  R0={R0} a={A} B0={B0} psi_edge={psi_edge:.4f} "
          f"q_axis={eqdata['qprof'][0]:.3f} q_edge={eqdata['qprof'][-1]:.3f} "
          f"Ip={eqdata['Ip']:.3e} A")


if __name__ == "__main__":
    main()
