#!/usr/bin/env python3
"""The bounce period of an orbit must not depend on where it is started.

A reference orbit is traced from a start point on the outboard midplane of the
circular fixture (itest_type=4). Points are then taken from that trajectory and
used as start points of the frequency scan (itest_type=5). Every one of them
lies on the same orbit, so the scan must return the reference bounce time and
toroidal shift. Start points off POTATO's Poincare cut used to close the orbit
after half a transit when the orbit moved towards the cut.
"""
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
R_REF = 180.0
PITCHES = (0.6, -0.6, 0.15, -0.15)
FRACTIONS = (0.15, 0.35, 0.65, 0.85)
MIN_OFFSET_CM = 1.0
REL_TOL = 1e-5


def write_input(path: Path, itest_type: int, r: float, z: float, xi: float) -> None:
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
  orbit_Rstart = {r:.17e}
  orbit_Zstart = {z:.17e}
  orbit_lambda = {xi:.17e}
  freq_Rmin = {r:.17e}
  freq_Rmax = {r:.17e}
  freq_n = 1
  profile_file = 'profile_poly.in'
/
"""
    )


def zero_electric_potential(path: Path) -> None:
    lines = path.read_text().splitlines()
    lines[-1] = " ".join(["0.0"] * 10)
    path.write_text("\n".join(lines) + "\n")


def run(executable: Path, work: Path, itest_type: int, r: float, z: float,
        xi: float) -> str:
    write_input(work / "potato.in", itest_type, r, z, xi)
    done = subprocess.run([executable], cwd=work, check=True, timeout=120,
                          capture_output=True, text=True)
    return done.stdout


def reference_orbit(executable: Path, work: Path, xi: float):
    stdout = run(executable, work, 4, R_REF, 0.0, xi)
    match = re.search(r"single orbit: taub=\s*(\S+)\s+delphi=\s*(\S+)\s+ierr=(\d+)",
                      stdout)
    if match is None or int(match.group(3)) != 0:
        raise AssertionError(f"reference orbit failed for xi={xi}")
    rows = []
    for line in (work / "fort.100").read_text().splitlines():
        values = line.split()
        if len(values) == 6 and "NaN" not in line:
            rows.append([float(v) for v in values])
    return float(match.group(1)), float(match.group(2)), rows


def scan_point(executable: Path, work: Path, r: float, z: float, xi: float):
    run(executable, work, 5, r, z, xi)
    lines = [line for line in (work / "freq_scan.dat").read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    values = lines[0].split()
    if int(values[6]) != 0:
        raise AssertionError(f"scan failed at R={r} Z={z} xi={xi}")
    return float(values[4]), float(values[5])


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: test_freq_scan_start.py <potato.x> <fixture-dir>")
        return 2

    executable = Path(sys.argv[1]).resolve()
    fixture = Path(sys.argv[2]).resolve()
    failures = []
    with tempfile.TemporaryDirectory(prefix="potato-freq-start-") as tmp:
        work = Path(tmp)
        for name in INPUT_NAMES:
            shutil.copy(fixture / name, work / name)
        zero_electric_potential(work / "profile_poly.in")
        for xi in PITCHES:
            taub_ref, delphi_ref, rows = reference_orbit(executable, work, xi)
            for fraction in FRACTIONS:
                r, _, z, _, xi_point, _ = rows[int(fraction * (len(rows) - 1))]
                if abs(z) < MIN_OFFSET_CM:
                    raise AssertionError(f"start point too close to the cut: Z={z}")
                taub, delphi = scan_point(executable, work, r, z, xi_point)
                ratio = taub / taub_ref
                dphi = abs(delphi - delphi_ref) / abs(delphi_ref)
                status = "ok"
                if abs(ratio - 1.0) > REL_TOL or dphi > REL_TOL:
                    status = "FAIL"
                    failures.append((xi, fraction))
                print(f"xi0={xi:+.2f} R={r:9.4f} Z={z:+8.4f} xi={xi_point:+.4f} "
                      f"taub/taub_ref={ratio:.8f} ddelphi={dphi:.2e} {status}")
    if failures:
        print(f"{len(failures)} start points disagree with the reference orbit")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
