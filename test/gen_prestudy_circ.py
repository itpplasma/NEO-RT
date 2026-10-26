#!/usr/bin/env python3
"""Boozer file of the rmp-proposal pre-study's circular model field.

    gen_prestudy_circ.py EFIT_TO_BOOZER_X OUTPUT_DIR

The pre-study (prestudies/neort-realspace/rsdrive_fields.py) uses R0 = 1,
F = R0 B0 = 1 and psi = B0/(2 Q2) log(1 + Q2 r^2/Q0) with Q0 = 1.5, Q2 = 8,
so q = (Q0 + Q2 r^2)/sqrt(1 - r^2).  This is the 'circ' family of
gen_analytic_eqdsk.py with other constants (SI units, R0 = 1 m, B0 = 1 T);
the boundary is r = 0.4 R0.  Writes prestudy.geqdsk and prestudy_boozer.bc.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


import gen_analytic_boozer as booz
import gen_analytic_eqdsk as eq

eq.R0, eq.A, eq.B0, eq.Q0 = 1.0, 0.4, 1.0, 1.5
eq.C = 8.0
eq.QA = eq.Q0 + eq.C * eq.A**2


def main(exe, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    geqdsk = outdir / "prestudy.geqdsk"
    eq.main("circ", geqdsk)
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        shutil.copy(geqdsk, work / "equilibrium.geqdsk")
        booz.write_wall(work / "convexwall.dat", r_axis_cm=100.0, radius_cm=55.0)
        (work / "field_divB0.inp").write_text(booz.FIELD_INP)
        (work / "efit_to_boozer.inp").write_text(booz.EFIT_INP.format(
            nsurf=booz.NSURF, mpol=booz.MPOL,
            psimax=booz.psi_boundary_cgs(geqdsk)))
        subprocess.run([exe], cwd=work, check=True, stdout=subprocess.DEVNULL)
        shutil.copy(work / "fromefit.bc", outdir / "prestudy_boozer.bc")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
