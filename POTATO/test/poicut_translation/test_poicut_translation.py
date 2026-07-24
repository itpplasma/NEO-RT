#!/usr/bin/env python3
"""Behavioral gate for the POTATO Poincare-cut equilibrium seed.

A rigid vertical translation of an axisymmetric equilibrium cannot change any
radial cut coordinate. Every vertical cut coordinate must change by exactly the
same translation. This catches origin-dependent seeds such as the historical
hard-coded Z=0 without encoding the implementation used to find the cut.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


SHIFT_M = 2.0
SHIFT_CM = 100.0 * SHIFT_M


def _parse_fixed_reals(lines: list[str]) -> list[float]:
    values: list[float] = []
    for line in lines:
        for start in range(0, len(line.rstrip("\n")), 16):
            token = line[start : start + 16].strip()
            if token:
                values.append(float(token.replace("D", "E")))
    return values


def _write_fixed_reals(stream, values: list[float]) -> None:
    for start in range(0, len(values), 5):
        stream.write("".join(f"{value:16.9E}" for value in values[start : start + 5]))
        stream.write("\n")


def translate_geqdsk(source: Path, destination: Path, shift_m: float) -> None:
    lines = source.read_text().splitlines(keepends=True)
    header = lines[0]
    fields = header.split()
    nw, nh = int(fields[-2]), int(fields[-1])
    prefix_count = 20 + 5 * nw + nw * nh

    # Locate the conventional integer boundary-count record independently of
    # line wrapping in the floating-point payload.
    float_lines = (prefix_count + 4) // 5
    count_line_index = 1 + float_lines
    count_fields = lines[count_line_index].split()
    if len(count_fields) != 2:
        raise AssertionError("missing GEQDSK boundary-count record")
    n_boundary, n_limiter = map(int, count_fields)

    equilibrium = _parse_fixed_reals(lines[1:count_line_index])
    if len(equilibrium) != prefix_count:
        raise AssertionError("unexpected GEQDSK equilibrium record length")
    geometry = _parse_fixed_reals(lines[count_line_index + 1 :])
    if len(geometry) != 2 * (n_boundary + n_limiter):
        raise AssertionError("unexpected GEQDSK boundary geometry length")

    # GEQDSK duplicates Zmaxis in the third and fifth scalar records.
    equilibrium[4] += shift_m
    equilibrium[6] += shift_m
    equilibrium[15] += shift_m
    geometry[1::2] = [value + shift_m for value in geometry[1::2]]

    with destination.open("w") as stream:
        stream.write(header)
        _write_fixed_reals(stream, equilibrium)
        stream.write(f"{n_boundary:5d}{n_limiter:5d}\n")
        _write_fixed_reals(stream, geometry)


def translate_wall(source: Path, destination: Path, shift_cm: float) -> None:
    wall = np.loadtxt(source)
    wall[:, 1] += shift_cm
    np.savetxt(destination, wall, fmt="%24.16e")


def run_probe(executable: Path, case_dir: Path) -> np.ndarray:
    completed = subprocess.run(
        [str(executable)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        timeout=120,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"Poincare-cut probe failed with {completed.returncode}:\n"
            f"{completed.stdout}"
        )
    records = [
        line for line in completed.stdout.splitlines() if line.startswith("POTATO_POICUT ")
    ]
    if len(records) != 1:
        raise AssertionError(f"missing unique POTATO_POICUT record:\n{completed.stdout}")
    values = np.fromstring(records[0].removeprefix("POTATO_POICUT "), sep=" ")
    if values.size != 10 or not np.all(np.isfinite(values)):
        raise AssertionError(f"invalid POTATO_POICUT record: {records[0]}")
    return values


def stage_case(fixture: Path, destination: Path) -> None:
    destination.mkdir()
    for name in ("potato.in", "field_divB0.inp", "profile_poly.in"):
        shutil.copy2(fixture / name, destination / name)


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_poicut_translation.py PROBE FIXTURE_DIR")
    executable = Path(sys.argv[1]).resolve()
    fixture = Path(sys.argv[2]).resolve()

    with tempfile.TemporaryDirectory(prefix="potato-poicut-") as tmp:
        root = Path(tmp)
        base = root / "base"
        shifted = root / "shifted"
        stage_case(fixture, base)
        stage_case(fixture, shifted)
        shutil.copy2(fixture / "circ.eqdsk", base / "circ.eqdsk")
        shutil.copy2(fixture / "convexwall.dat", base / "convexwall.dat")
        translate_geqdsk(fixture / "circ.eqdsk", shifted / "circ.eqdsk", SHIFT_M)
        translate_wall(
            fixture / "convexwall.dat", shifted / "convexwall.dat", SHIFT_CM
        )

        original = run_probe(executable, base)
        translated = run_probe(executable, shifted)

    radial = np.array([0, 2, 4, 6, 8])
    vertical = np.array([1, 3, 5, 7, 9])
    np.testing.assert_allclose(translated[radial], original[radial], rtol=0.0, atol=2e-7)
    np.testing.assert_allclose(
        translated[vertical] - original[vertical],
        SHIFT_CM,
        rtol=0.0,
        atol=2e-7,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
