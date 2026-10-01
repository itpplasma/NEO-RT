#!/usr/bin/env python3
"""Regenerate the Boozer-side fixtures for an analytic GEQDSK.

    gen_analytic_boozer.py EFIT_TO_BOOZER_X FIXTURE_DIR NAME

FIXTURE_DIR must contain NAME.geqdsk (written by gen_analytic_eqdsk.py).
EFIT_TO_BOOZER_X is libneo's converter at the pinned LIBNEO_REF (target
efit_to_boozer.x in the NEO-RT build tree).  The script runs it in a
temporary directory and writes

  NAME_boozer.bc  axisymmetric Boozer file of the same field (inp_swi=9)
  NAME_pert.bc    synthetic n=-2 Boozer perturbation on the same s grid

The perturbation amplitudes are an arbitrary smooth closed form; only the
fact that both backends read the same file matters for the comparison.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

NSURF = 48
MPOL = 16
NTOR_PERT = -2
M_PERT = np.arange(-6, 13)
B_REF = 2.0  # T, scale of the perturbation amplitudes
EFIT_INP = """3600      nstep
500       nlabel
256       ntheta
1000      nsurfmax
{nsurf}        nsurf
{mpol}        mpol
{psimax:.10e} psimax
"""
FIELD_INP = """0           ipert
1           iequil
1.00        ampl
72          ntor
0.99        cutoff
4           icftype
'equilibrium.geqdsk'  gfile
'unused'    pfile
'convexwall.dat' convexfile
'unused'    fluxdatapath
0           nwindow_r
0           nwindow_z
1           ieqfile
"""


def psi_boundary_cgs(geqdsk):
    lines = Path(geqdsk).read_text().splitlines()
    sibry = float(lines[2][48:64])
    return sibry * 1.0e8


def write_wall(path, r_axis_cm=165.0, radius_cm=80.0):
    theta = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False)
    wall = np.column_stack([r_axis_cm + radius_cm * np.cos(theta),
                            radius_cm * np.sin(theta)])
    np.savetxt(path, wall, fmt="%24.16e")


def run_efit_to_boozer(exe, fixture_dir, name):
    geqdsk = fixture_dir / f"{name}.geqdsk"
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        shutil.copy(geqdsk, work / "equilibrium.geqdsk")
        write_wall(work / "convexwall.dat")
        (work / "field_divB0.inp").write_text(FIELD_INP)
        psimax = psi_boundary_cgs(geqdsk)
        (work / "efit_to_boozer.inp").write_text(
            EFIT_INP.format(nsurf=NSURF, mpol=MPOL, psimax=psimax))
        subprocess.run([exe], cwd=work, check=True, stdout=subprocess.DEVNULL)
        shutil.copy(work / "fromefit.bc", fixture_dir / f"{name}_boozer.bc")


def surface_params(bc_path):
    lines = Path(bc_path).read_text().splitlines()
    header = lines[:6]
    params = [lines[i + 2] for i in range(len(lines))
              if lines[i].strip().startswith("s ")]
    return header, params


def perturbation_amplitudes(s, m):
    envelope = 1.0e-3 * B_REF * s ** (0.5 * abs(m)) * (1.0 - s)
    shape = 1.0 / (1.0 + (m - 4.0) ** 2 / 8.0)
    return envelope * shape, 0.5 * envelope * shape * ((m % 3) - 1)


def write_perturbation(fixture_dir, name):
    header, params = surface_params(fixture_dir / f"{name}_boozer.bc")
    nm = M_PERT.size
    out = header[:4]
    out.append(" m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]")
    out.append(f"{nm - 1:6d}{0:6d}{len(params):6d}{1:6d}  0.000000e+00   0.000000e+00   0.000000e+00")
    for row in params:
        s = float(row.split()[0])
        out.append("        s               iota           Jpol/nper          Itor"
                   "            pprime         sqrt g(0,0)")
        out.append("                                          [A]           [A]"
                   "             [Pa]         (dV/ds)/nper")
        out.append(row)
        out.append("    m    n      rmnc [m]         rmns [m]         zmnc [m]"
                   "         zmns [m]         vmnc [ ]         vmns [ ]"
                   "         bmnc [T]         bmns [T]")
        for m in M_PERT:
            bc, bs = perturbation_amplitudes(s, m)
            out.append(f"{m:5d}{NTOR_PERT:5d}" + "".join(
                f"{v:17.8e}" for v in (0, 0, 0, 0, 0, 0, bc, bs)))
    (fixture_dir / f"{name}_pert.bc").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    fixtures = Path(sys.argv[2])
    run_efit_to_boozer(sys.argv[1], fixtures, sys.argv[3])
    write_perturbation(fixtures, sys.argv[3])
