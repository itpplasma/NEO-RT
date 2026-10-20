#!/usr/bin/env python3
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


INPUT_NAMES = (
    "circ.eqdsk",
    "convexwall.dat",
    "field_divB0.inp",
    "profile_poly.in",
)
EXPECTED_COLUMNS = (
    "id",
    "R_cm",
    "Z_cm",
    "p",
    "xi",
    "sigma",
    "H0",
    "J_perp",
    "psi_star",
    "psi_axis",
    "psi_edge",
    "phi_elec",
    "v0_cm_s",
    "status",
)


def write_input(path: Path, itest_type: int = 4, scan: bool = False) -> None:
    scan_input = """
  freq_Rmin = 145d0
  freq_Rmax = 175d0
  freq_n = 3
""" if scan else ""
    path.write_text(
        f"""&potato_nml
  itest_type = {itest_type}
  E_alpha = 5d3
  A_alpha = 2d0
  Z_alpha = 1d0
  rho_pol_max = 0.65d0
  Rmax_orbit = 250d0
  ntimstep = 2000
  npoicut = 1000
  orbit_Rstart = 180d0
  orbit_Zstart = 0d0
  orbit_lambda = 0.3d0
  profile_file = 'profile_poly.in'
{scan_input}
/
"""
    )


def zero_electric_potential(path: Path) -> None:
    lines = path.read_text().splitlines()
    lines[-1] = " ".join(["0.0"] * 10)
    path.write_text("\n".join(lines) + "\n")


def continue_circular_flux_outside_lcfs(path: Path) -> None:
    """Give the synthetic circular field a monotone exterior flux continuation.

    The stored fixture clamps psi at its r=.5m LCFS. An exterior field is
    legitimate, while its plasma profile remains supported only on s=0..1.
    Keep every interior grid value and the LCFS/profile unchanged.
    """
    lines = path.read_text().splitlines()
    nr, nz = map(int, lines[0].split()[-2:])
    header = [float(value) for line in lines[1:5] for value in line.split()]
    first = 1 + (20 + 4*nr)//5
    flux = [float(value) for line in lines[first:first+nr*nz//5]
            for value in line.split()]
    for j in range(nz):
        z = header[4] - header[1]/2 + j*header[1]/(nz-1)
        for i in range(nr):
            r = header[3] + i*header[0]/(nr-1)
            radius = math.hypot(r-header[5], z-header[6])
            if radius > .5:
                flux[j*nr+i] = header[8] + .4*(radius-.5)
    lines[first:first+nr*nz//5] = [
        "".join(f"{value:16.8E}" for value in flux[k:k+5])
        for k in range(0, len(flux), 5)
    ]
    path.write_text("\n".join(lines)+"\n")


def read_handoff(path: Path) -> tuple[list[str], list[list[float]]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if lines[0] != "# potato-simple-invariants-v1":
        raise AssertionError(f"unexpected handoff version: {lines[0]}")
    columns = lines[1].removeprefix("# ").split()
    rows = [[float(value) for value in line.split()] for line in lines[2:]]
    return columns, rows


def require_valid_rows(columns: list[str], rows: list[list[float]]) -> None:
    if tuple(columns) != EXPECTED_COLUMNS:
        raise AssertionError(f"unexpected columns: {columns}")
    for values in rows:
        if len(values) != len(EXPECTED_COLUMNS):
            raise AssertionError(f"unexpected row width: {len(values)}")
        row = dict(zip(columns, values))
        if not all(math.isfinite(value) for value in values):
            raise AssertionError("handoff contains non-finite values")
        if row["status"] != 0 or row["sigma"] != 1:
            raise AssertionError(f"unexpected status/sign: {row}")
        if not math.isclose(row["H0"], 1.0, rel_tol=0.0, abs_tol=1e-14):
            raise AssertionError(f"unexpected H0: {row['H0']}")
        if row["J_perp"] <= 0.0 or row["psi_edge"] == row["psi_axis"]:
            raise AssertionError(f"invalid invariants: {row}")
        if abs(row["phi_elec"]) > 1e-12 or row["v0_cm_s"] <= 0.0:
            raise AssertionError(f"invalid normalization metadata: {row}")


def require_frequency_output(path: Path, expected_rows: int) -> None:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    expected_header = "# R_start[cm] rho_pol omega_b[1/s] omega_phi[1/s] taub delphi ierr"
    if lines[0] != expected_header:
        raise AssertionError(f"unexpected frequency header: {lines[0]}")
    if len(lines) != expected_rows + 1:
        raise AssertionError(f"unexpected frequency row count: {len(lines) - 1}")
    if any(len(line.split()) != 7 for line in lines[1:]):
        raise AssertionError("frequency output schema changed")


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_invariant_handoff.py <potato.x> <fixture-dir>")
        return 2

    executable = Path(sys.argv[1]).resolve()
    fixture = Path(sys.argv[2]).resolve()
    with tempfile.TemporaryDirectory(prefix="potato-invariants-") as tmp:
        work = Path(tmp)
        for name in INPUT_NAMES:
            shutil.copy(fixture / name, work / name)
        zero_electric_potential(work / "profile_poly.in")
        write_input(work / "potato.in")
        subprocess.run([executable], cwd=work, check=True, timeout=120)
        columns, rows = read_handoff(work / "potato_invariants.dat")
        require_valid_rows(columns, rows)
        if len(rows) != 1:
            raise AssertionError(f"unexpected single-orbit row count: {len(rows)}")

        write_input(work / "potato.in", itest_type=5, scan=True)
        subprocess.run([executable], cwd=work, check=True, timeout=120)
        columns, rows = read_handoff(work / "potato_invariants.dat")
        require_valid_rows(columns, rows)
        if [row[0] for row in rows] != [1.0, 2.0, 3.0]:
            raise AssertionError("frequency handoff identifiers changed")
        require_frequency_output(work / "freq_scan.dat", 3)

        continue_circular_flux_outside_lcfs(work / "circ.eqdsk")
        # Axis initialization first probes R=1 cm outside this circular grid.
        # It needs only B and must not evaluate the strictly supported profile.
        (work / "profile_poly.in").write_text(
            "POTATO_PROFILE_TABLE_V1\n5\n"
            + "".join(f"{s} 2e13 5000 0\n" for s in (0, .25, .5, .75, 1))
        )
        write_input(work / "potato.in")
        subprocess.run([executable], cwd=work, check=True, timeout=120)
        columns, rows = read_handoff(work / "potato_invariants.dat")
        require_valid_rows(columns, rows)
        if len(rows) != 1:
            raise AssertionError("tabulated-profile orbit handoff is missing")
    return 0


if __name__ == "__main__":
    sys.exit(main())
