#!/usr/bin/env python3
"""
Create minimal input files for golden record tests.

Reads full input files and creates minimal versions containing only
the data points needed for the test s values (0.1, 0.2, ..., 0.9).

Usage:
    1. Place full input files in test/golden_record/input/
    2. Run: python create_minimal_input.py
    3. Minimal files replace the originals in input/
    4. Commit the minimal files to git
"""

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

# Target s values for different file types
# Boozer files need more points for accurate magnetic field integration
S_TARGET_BOOZER = np.linspace(0.0, 1.0, 101)
# Plasma and profile files can use fewer points
S_TARGET_PROFILES = np.linspace(0.0, 1.0, 21)

INPUT_DIR = Path(__file__).parent / "input"


def read_plasma_in(path: Path) -> tuple[dict, np.ndarray]:
    """Read plasma.in file.

    Returns header dict and data array (nrows, 6).
    """
    with open(path) as f:
        lines = f.readlines()

    # Line 2 has: nplasma am1 am2 Z1 Z2
    parts = lines[1].split()
    header = {
        "nplasma": int(parts[0]),
        "am1": float(parts[1]),
        "am2": float(parts[2]),
        "Z1": float(parts[3]),
        "Z2": float(parts[4]),
    }

    # Data starts at line 4 (0-indexed: line 3)
    data = np.loadtxt(path, skiprows=3)
    return header, data


def write_plasma_in(path: Path, header: dict, data: np.ndarray) -> None:
    """Write plasma.in file."""
    with open(path, "w") as f:
        f.write("% N am1 am2 Z1 Z2\n")
        f.write(
            f"{len(data)} {header['am1']:.15g} {header['am2']:.15g} "
            f"{header['Z1']:.1f} {header['Z2']:.1f}\n"
        )
        f.write("% s ni_1[cm^-3] ni_2[cm^-3] Ti_1[eV] Ti_2[eV] Te[eV]\n")
        for row in data:
            f.write(
                f" {row[0]:.10e}  {row[1]:.10e}  {row[2]:.10e}  "
                f"{row[3]:.10e}  {row[4]:.10e}  {row[5]:.10e}\n"
            )


def process_plasma() -> None:
    """Process plasma.in file."""
    path = INPUT_DIR / "plasma.in"
    print(f"Processing {path}...")

    header, data = read_plasma_in(path)
    s_orig = data[:, 0]

    # Clamp target s to original data range
    s_min, s_max = s_orig.min(), s_orig.max()
    s_target = np.clip(S_TARGET_PROFILES, s_min, s_max)

    # Build splines for each column and evaluate
    new_data = np.zeros((len(s_target), 6))
    new_data[:, 0] = s_target
    for i in range(1, 6):
        spl = CubicSpline(s_orig, data[:, i], bc_type="natural")
        new_data[:, i] = spl(s_target)

    write_plasma_in(path, header, new_data)
    print(f"  Reduced from {len(data)} to {len(new_data)} rows")


def read_profile_in(path: Path) -> np.ndarray:
    """Read profile.in file (no header, 3 columns)."""
    return np.loadtxt(path)


def write_profile_in(path: Path, data: np.ndarray) -> None:
    """Write profile.in file."""
    np.savetxt(path, data, fmt=" %.10e")


def process_profile() -> None:
    """Process profile.in file."""
    path = INPUT_DIR / "profile.in"
    print(f"Processing {path}...")

    data = read_profile_in(path)
    s_orig = data[:, 0]

    # Clamp target s to original data range
    s_min, s_max = s_orig.min(), s_orig.max()
    s_target = np.clip(S_TARGET_PROFILES, s_min, s_max)

    # Build splines for each column and evaluate
    new_data = np.zeros((len(s_target), data.shape[1]))
    new_data[:, 0] = s_target
    for i in range(1, data.shape[1]):
        spl = CubicSpline(s_orig, data[:, i], bc_type="natural")
        new_data[:, i] = spl(s_target)

    write_profile_in(path, new_data)
    print(f"  Reduced from {len(data)} to {len(new_data)} rows")


def read_boozer(path: Path) -> tuple[dict, list[str], list[dict]]:
    """Read Boozer coordinate file.

    Returns:
        header: dict with m0b, n0b, nsurf, nper, flux, a, R0
        comment_lines: list of 4 CC comment lines
        surfaces: list of dicts, each with s, iota, curr_pol, curr_tor,
                  pprime, sqrt_g, modes (array of shape (nmode, 10))
    """
    with open(path) as f:
        lines = f.readlines()

    # First 4 lines are comments
    comment_lines = lines[:4]

    # Line 5 is header labels, line 6 has values
    header_parts = lines[5].split()
    header = {
        "m0b": int(header_parts[0]),
        "n0b": int(header_parts[1]),
        "nsurf": int(header_parts[2]),
        "nper": int(header_parts[3]),
        "flux": float(header_parts[4]),
        "a": float(header_parts[5]),
        "R0": float(header_parts[6]),
    }

    nmode = (header["m0b"] + 1) * (2 * header["n0b"] + 1)
    if header["n0b"] == 0:
        nmode = header["m0b"] + 1

    surfaces = []
    line_idx = 6  # Start after header

    for _ in range(header["nsurf"]):
        # Skip 2 header lines (labels + units)
        line_idx += 2
        # Read surface parameters
        surf_parts = lines[line_idx].split()
        surface = {
            "s": float(surf_parts[0]),
            "iota": float(surf_parts[1]),
            "curr_pol": float(surf_parts[2]),
            "curr_tor": float(surf_parts[3]),
            "pprime": float(surf_parts[4]),
            "sqrt_g": float(surf_parts[5]),
        }
        line_idx += 1

        # Skip mode header line
        line_idx += 1

        # Read mode data
        modes = np.zeros((nmode, 10))
        for j in range(nmode):
            mode_parts = lines[line_idx].split()
            modes[j, 0] = int(mode_parts[0])  # m
            modes[j, 1] = int(mode_parts[1])  # n
            for k in range(8):
                modes[j, 2 + k] = float(mode_parts[2 + k])
            line_idx += 1

        surface["modes"] = modes
        surfaces.append(surface)

    return header, comment_lines, surfaces


def write_boozer(
    path: Path, header: dict, comment_lines: list[str], surfaces: list[dict]
) -> None:
    """Write Boozer coordinate file."""
    with open(path, "w") as f:
        # Write comment lines
        for line in comment_lines:
            f.write(line)

        # Write header
        f.write(
            " m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n"
        )
        f.write(
            f"  {header['m0b']:2d}   {header['n0b']:3d}  {len(surfaces):4d}   "
            f"{header['nper']:3d}  {header['flux']:14.8e}   {header['a']:.8f}   "
            f"{header['R0']:.8f}\n"
        )

        for surf in surfaces:
            # Write surface header
            f.write(
                "        s               iota           Jpol/nper          Itor"
                "            pprime         sqrt g(0,0)\n"
            )
            f.write(
                "                                          [A]           [A]"
                "             [Pa]         (dV/ds)/nper\n"
            )
            f.write(
                f"   {surf['s']:.8e}   {surf['iota']:.8e}   {surf['curr_pol']:.8e}"
                f"   {surf['curr_tor']:.8e}  {surf['pprime']:.8e}  "
                f"{surf['sqrt_g']:.8e}\n"
            )

            # Write mode header
            f.write(
                "    m    n      rmnc [m]         rmns [m]         zmnc [m]"
                "         zmns [m]         vmnc [ ]         vmns [ ]"
                "         bmnc [T]         bmns [T]\n"
            )

            # Write mode data
            for mode in surf["modes"]:
                m, n = int(mode[0]), int(mode[1])
                f.write(f"  {m:3d}  {n:3d}")
                for k in range(8):
                    f.write(f"   {mode[2 + k]:.8e}")
                f.write("\n")


def process_boozer(name: str) -> None:
    """Process a Boozer coordinate file."""
    path = INPUT_DIR / name
    print(f"Processing {path}...")

    header, comment_lines, surfaces = read_boozer(path)
    nmode = len(surfaces[0]["modes"])

    # Collect all s values
    s_orig = np.array([s["s"] for s in surfaces])

    # Clamp target s to original data range
    s_min, s_max = s_orig.min(), s_orig.max()
    s_target = np.clip(S_TARGET_BOOZER, s_min, s_max)

    # Build splines for surface quantities
    iota_orig = np.array([s["iota"] for s in surfaces])
    curr_pol_orig = np.array([s["curr_pol"] for s in surfaces])
    curr_tor_orig = np.array([s["curr_tor"] for s in surfaces])
    pprime_orig = np.array([s["pprime"] for s in surfaces])
    sqrt_g_orig = np.array([s["sqrt_g"] for s in surfaces])

    spl_iota = CubicSpline(s_orig, iota_orig, bc_type="natural")
    spl_curr_pol = CubicSpline(s_orig, curr_pol_orig, bc_type="natural")
    spl_curr_tor = CubicSpline(s_orig, curr_tor_orig, bc_type="natural")
    spl_pprime = CubicSpline(s_orig, pprime_orig, bc_type="natural")
    spl_sqrt_g = CubicSpline(s_orig, sqrt_g_orig, bc_type="natural")

    # Build splines for each mode's 8 floating-point values
    mode_splines = []
    for mode_idx in range(nmode):
        mode_data = np.array([s["modes"][mode_idx, 2:] for s in surfaces])
        splines = [
            CubicSpline(s_orig, mode_data[:, k], bc_type="natural") for k in range(8)
        ]
        mode_splines.append(splines)

    # Get mode numbers (m, n) from first surface (they don't change)
    mode_mn = surfaces[0]["modes"][:, :2].astype(int)

    # Create new surfaces at target s values
    new_surfaces = []
    for s in s_target:
        surf = {
            "s": s,
            "iota": float(spl_iota(s)),
            "curr_pol": float(spl_curr_pol(s)),
            "curr_tor": float(spl_curr_tor(s)),
            "pprime": float(spl_pprime(s)),
            "sqrt_g": float(spl_sqrt_g(s)),
        }
        modes = np.zeros((nmode, 10))
        modes[:, :2] = mode_mn
        for mode_idx in range(nmode):
            for k in range(8):
                modes[mode_idx, 2 + k] = mode_splines[mode_idx][k](s)
        surf["modes"] = modes
        new_surfaces.append(surf)

    write_boozer(path, header, comment_lines, new_surfaces)
    print(f"  Reduced from {len(surfaces)} to {len(new_surfaces)} surfaces")


def main() -> None:
    """Process all input files."""
    print("Creating minimal input files for golden record tests\n")

    process_plasma()
    process_profile()
    process_boozer("in_file")
    process_boozer("in_file_pert")

    print("\nMinimal input files created successfully in input/")


if __name__ == "__main__":
    main()
