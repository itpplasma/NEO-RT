#!/usr/bin/env python3
import math
import re
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
    "rho_pol",
    "R_gc",
    "Z_gc",
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
    "taub",
    "delphi",
    "omega_b",
    "omega_phi",
    "ierr",
)


def zero_electric_potential(path: Path) -> None:
    lines = path.read_text().splitlines()
    lines[-1] = " ".join(["0.0"] * 10)
    path.write_text("\n".join(lines) + "\n")


def write_input(path: Path, extra: str) -> None:
    path.write_text(
        f"""&potato_nml
  itest_type = 4
  E_alpha = 5d3
  A_alpha = 2d0
  Z_alpha = 1d0
  rho_pol_max = 0.65d0
  Rmax_orbit = 250d0
  ntimstep = 2000
  npoicut = 1000
  profile_file = 'profile_poly.in'
{extra}
/
"""
    )


def read_contour(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    columns = lines[0].removeprefix("# ").split()
    rows = [
        dict(zip(columns, (float(value) for value in line.split())))
        for line in lines[1:]
    ]
    return columns, rows


def parse_single_orbit(text: str) -> tuple[float, float, int]:
    match = re.search(
        r"single orbit: taub=\s*(\S+)\s+delphi=\s*(\S+)\s+ierr=(\d+)",
        text,
    )
    if match is None:
        raise AssertionError("single-orbit summary is missing")
    return float(match.group(1)), float(match.group(2)), int(match.group(3))


def assert_close(actual: float, expected: float, tolerance: float, name: str) -> None:
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{name}: expected {expected}, got {actual}")


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "usage: test_fixed_energy_contour.py "
            "<contour.x> <potato.x> <fixture-dir>"
        )
        return 2

    contour_executable = Path(sys.argv[1]).resolve()
    potato_executable = Path(sys.argv[2]).resolve()
    fixture = Path(sys.argv[3]).resolve()
    with tempfile.TemporaryDirectory(prefix="potato-fixed-contour-") as tmp:
        work = Path(tmp)
        for name in INPUT_NAMES:
            shutil.copy(fixture / name, work / name)
        zero_electric_potential(work / "profile_poly.in")
        write_input(
            work / "potato.in",
            """  contour_rho_min = 0.45d0
  contour_rho_max = 0.55d0
  contour_nrho = 2
  contour_jperp_min = 2d-5
  contour_jperp_max = 3d-5
  contour_njperp = 2
  contour_enkin = 1d0
  contour_sigma = 1
""",
        )
        subprocess.run([contour_executable], cwd=work, check=True, timeout=120)
        columns, rows = read_contour(work / "potato_resonance_contour.dat")

        if tuple(columns) != EXPECTED_COLUMNS:
            raise AssertionError(f"unexpected columns: {columns}")
        if len(rows) != 4:
            raise AssertionError(f"expected a rectangular 2x2 grid, got {len(rows)}")
        if {row["rho_pol"] for row in rows} != {0.45, 0.55}:
            raise AssertionError("configured rho_pol endpoints are missing")
        if {row["J_perp"] for row in rows} != {2.0e-5, 3.0e-5}:
            raise AssertionError("configured J_perp endpoints are missing")

        successful = [row for row in rows if row["ierr"] == 0.0]
        if not successful:
            raise AssertionError("fixed-energy contour has no successful orbit")
        row = successful[0]
        assert_close(row["H0"], 1.0, 1.0e-12, "fixed normalized energy")
        assert_close(row["p"], 1.0, 1.0e-12, "fixed normalized momentum")
        if row["psi_axis"] == row["psi_edge"]:
            raise AssertionError("contour flux gauge is degenerate")
        if row["phi_elec"] != 0.0 or row["v0_cm_s"] <= 0.0:
            raise AssertionError("contour normalization metadata is invalid")
        assert_close(
            row["omega_phi"] / row["omega_b"],
            row["delphi"] / (2.0 * math.pi),
            1.0e-8,
            "frequency normalization",
        )

        write_input(
            work / "potato.in",
            f"""  orbit_Rstart = {row['R_gc']:.17e}
  orbit_Zstart = {row['Z_gc']:.17e}
  orbit_lambda = {row['xi']:.17e}
""",
        )
        single = subprocess.run(
            [potato_executable],
            cwd=work,
            check=True,
            timeout=120,
            capture_output=True,
            text=True,
        )
        taub, delphi, ierr = parse_single_orbit(single.stdout)
        if ierr != 0:
            raise AssertionError(f"equivalent single orbit failed with ierr={ierr}")
        assert_close(taub, row["taub"], 2.0e-6, "single-orbit taub")
        assert_close(delphi, row["delphi"], 2.0e-6, "single-orbit delphi")

    return 0


if __name__ == "__main__":
    sys.exit(main())
