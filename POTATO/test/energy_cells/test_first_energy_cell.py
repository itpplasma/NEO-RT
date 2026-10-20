#!/usr/bin/env python3
"""A one-cell thermal integral equals its fixed-energy slice times Delta H/T.

At zero potential, nenerg=1 and thermen_max=3 put the midpoint at H/T=1.5.
The fixed-energy path supplies an independent sampling/normalization oracle:
it must find the same resonances, and the thermal torque must be three times
the delta-normalized slice. This also checks that cell 1 reaches output merging.
"""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


CASE = Path(__file__).resolve().parent.parent / "golden_record_resonance"
FILES = ("potato.in", "field_divB0.inp", "convexwall.dat",
         "profile_poly.in", "circ.eqdsk", "bmod_n.dat")


def run(executable: Path, work: Path, energy_ratio: float):
    work.mkdir()
    for name in FILES:
        shutil.copy(CASE / name, work / name)
    profile = work / "profile_poly.in"
    lines = profile.read_text().splitlines()
    lines[-1] = " ".join(["0.0"] * 10)
    profile.write_text("\n".join(lines) + "\n")
    namelist = work / "potato.in"
    text = re.sub(r"nenerg\s*=\s*\d+", "nenerg = 1", namelist.read_text())
    text = re.sub(r"thermen_max\s*=\s*[^\s,]+", "thermen_max = 3d0", text)
    text = text.replace("/", f"  monoenergetic_x = {energy_ratio}\n/")
    namelist.write_text(text)
    result = subprocess.run(
        [str(executable)], cwd=work,
        env=dict(os.environ, OMP_NUM_THREADS="2"),
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-1000:]
    torque_path = work / "integral_torque.dat"
    assert torque_path.is_file(), (
        "first energy cell did not reach torque output:\n" + result.stdout[-1000:]
    )
    torque = float(torque_path.read_text())
    assert math.isfinite(torque) and torque != 0.0, f"invalid torque: {torque}"
    resonances = sorted((work / "fort.31415").read_text().splitlines())
    assert resonances, "first energy cell has no resonances"
    boxes = [float(line.split()[0]) for line in
             (work / "boxcounted_torque.dat").read_text().splitlines()]
    return torque, resonances, boxes


def main() -> None:
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="potato-first-cell-") as directory:
        root = Path(directory)
        thermal, thermal_lines, thermal_boxes = run(executable, root / "thermal", -1.0)
        fixed, fixed_lines, fixed_boxes = run(executable, root / "fixed", 1.5)
    assert thermal_lines == fixed_lines, "thermal/fixed-energy resonances differ"
    assert math.isclose(thermal, 3.0*fixed, rel_tol=1e-10), (thermal, fixed)
    assert len(thermal_boxes) == len(fixed_boxes)
    for actual, reference in zip(thermal_boxes, fixed_boxes):
        assert math.isclose(actual, 3.0*reference, rel_tol=1e-10, abs_tol=1e-12)
    print(f"PASS: cell 1 torque {thermal} = 3 * fixed-energy torque {fixed}")


if __name__ == "__main__":
    main()
