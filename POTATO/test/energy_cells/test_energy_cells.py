#!/usr/bin/env python3
"""Every cell of the midpoint energy grid is integrated by the resonant torque.

The energy integral of resonant_int.f90 is a midpoint rule over nenerg cells of the total-energy range
[toten_min, toten_max]:  toten_i = toten_min + (toten_max - toten_min) (i - 1/2) / nenerg,  i = 1..nenerg.
An earlier loop started at i = 2 and silently dropped the lowest cell (a 17 % error of the torque at nenerg = 12,
because that cell carries the opposite sign of the thermodynamic-force weight).

Oracle: the grid formula above, applied to the "minimum/maximum total energy" values potato.x prints, and the
energy slices of the merged resonance-line file fort.31415 (first column).  Usage:
  test_energy_cells.py POTATO_EXE     (default: POTATO/build/potato.x)
Runs the committed synthetic circular case of golden_record_resonance with nenerg = 3 (about 10 s).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CASE = Path(__file__).resolve().parent.parent / "golden_record_resonance"
FILES = ("potato.in", "field_divB0.inp", "convexwall.dat", "profile_poly.in", "circ.eqdsk", "bmod_n.dat")


def main() -> int:
    exe = Path(sys.argv[1] if len(sys.argv) > 1 else CASE.parent.parent / "build" / "potato.x").resolve()
    if not exe.is_file():
        print(f"potato.x not found: {exe}")
        return 77
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for name in FILES:
            shutil.copy(CASE / name, work / name)
        nml = (work / "potato.in").read_text()
        nenerg = int(re.search(r"nenerg\s*=\s*(\d+)", nml).group(1))
        run = subprocess.run([str(exe)], cwd=work, env=dict(os.environ, OMP_NUM_THREADS="2"),
                             capture_output=True, text=True, timeout=300)
        if run.returncode != 0:
            print(run.stdout[-2000:], run.stderr[-1000:])
            return 1
        out = run.stdout
        tmin = float(re.search(r"minimum total energy =\s*(\S+)", out).group(1))
        tmax = float(re.search(r"maximum total energy =\s*(\S+)", out).group(1))
        expected = [tmin + (tmax - tmin) * (i - 0.5) / nenerg for i in range(1, nenerg + 1)]
        used = sorted(float(x) for x in re.findall(r"^toten =\s*(\S+)", out, re.M))
        slices = sorted({float(line.split()[0]) for line in (work / "fort.31415").read_text().splitlines() if line.strip()})
    ok = len(used) == nenerg and all(abs(u - e) <= 1e-6 * abs(tmax) for u, e in zip(used, expected))
    ok = ok and len(slices) == nenerg and all(abs(u - e) <= 1e-6 * abs(tmax) for u, e in zip(slices, expected))
    print("expected", expected, "integrated", used, "resonance-line slices", slices)
    print("PASS" if ok else "FAIL: an energy cell of the midpoint rule is missing")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
