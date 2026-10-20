#!/usr/bin/env python3
"""Qualify the native Boozer field probe against manufactured type-9 data.

Invocation: qualify_probe.py BUILD/neo_rt_field_probe.x --output-dir DIR
The circular geometry fixes chart orientation only; this interpolation fixture
does not claim to solve an equilibrium. All sampled physical coefficients are
linear in s, so spline interpolation must reproduce their analytic values.
Conventions: NEO-RT src/do_magfie_standalone.f90 and
test/test_boozer_handedness.f90 (SI .bc input, cgs native field output).
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


COLUMNS = ("s_tor", "phi_B_rad", "theta_B_rad", "Bmod_G", "Bphcov_Gcm",
           "Bthcov_Gcm", "iota", "sqrtg_cm3", "psi_pr_Gcm2")
SURFACES = (0.1, 0.3, 0.5, 0.7, 0.9)
TOLERANCE = 1.0e-11
POINTS = [(s, phi, theta) for s in (0.5, 0.41) for phi in (0.0, 0.73)
          for theta in (0.0, 0.37, math.pi / 2, math.pi, -1.23, 2 * math.pi + 0.37)]


def parameters(s, iota_sign):
    return iota_sign * (0.23 + 0.07 * s), 2.5e7 * (1 + 0.2 * s), 3.0e5 * (0.5 + 0.3 * s)


def volume_derivative(s, orientation, flux_sign, iota_sign):
    """Return 4*pi**2*<sqrt(g)> in m^3, the consistent (dV/ds)/nper column.

    Uses <1/(c + r*cos(t - d))**2> = c/(c**2 - r**2)**1.5 for |B| below.
    """
    iota, jpol, itor = parameters(s, iota_sign)
    bph, bth = -0.2 * orientation * jpol, -0.2 * orientation * itor
    psi = flux_sign * 2.5e7
    r = s * math.hypot(0.1, 0.03)
    inv_b2_avg = 1.0e-8 * 2 / (4 - r**2) ** 1.5
    return 4 * math.pi**2 * (-psi * (bph + iota * bth)) * inv_b2_avg * 1.0e-6


def write_fixture(path, orientation, flux_sign, iota_sign):
    flux = flux_sign * math.pi / 2  # Wb: psi_pr = +/-2.5e7 G cm^2.
    lines = ["CC manufactured analytic field", "CC geometry determines orientation",
             "CC m0b n0b nsurf nper flux a R", "CC", "CC",
             f"1 0 {len(SURFACES)} 1 {flux:.17e} 0.5 5.0"]
    for s in SURFACES:
        iota, jpol, itor = parameters(s, iota_sign)
        radius = 0.5 * math.sqrt(s)
        dvds = volume_derivative(s, orientation, flux_sign, iota_sign)
        lines.extend(("CC s iota Jpol Itor pprime sqrtg00", "CC units A A Pa m3",
                      f"{s:.17e} {iota:.17e} {jpol:.17e} {itor:.17e} 0 {dvds:.17e}",
                      "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns",
                      "0 0 5 0 0 0 0 0 2 0",
                      f"1 0 {radius:.17e} 0 0 {orientation * radius:.17e} 0 0 "
                      f"{0.1 * s:.17e} {0.03 * s:.17e}"))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def expected(point, orientation, flux_sign, iota_sign):
    s, phi, theta = point
    iota, jpol, itor = parameters(s, iota_sign)
    # Ampere -> G cm: mu0/(2*pi)*1e6 = 0.2. Geometry gives h=-orientation.
    bph, bth = -0.2 * orientation * jpol, -0.2 * orientation * itor
    bmod = 1.0e4 * (2 + 0.1 * s * math.cos(theta) + 0.03 * s * math.sin(theta))
    psi = flux_sign * 2.5e7
    # Boozer Jacobian identity with NEO-RT sign_theta=-1.
    sqrtg = -psi * (bph + iota * bth) / bmod**2
    return (s, phi, theta, bmod, bph, bth, iota, sqrtg, psi)


def write_points(path, points):
    path.write_text("".join(" ".join(f"{x:.17e}" for x in p) + "\n" for p in points),
                    encoding="utf-8")


def run_probe(probe, bc, points, output, timeout):
    output.unlink(missing_ok=True)
    command = [str(probe), str(bc), "9", str(points), str(output)]
    result = subprocess.run(command, cwd=output.parent, text=True, capture_output=True,
                            timeout=timeout, check=False)
    output.with_suffix(".stdout").write_text(result.stdout, encoding="utf-8")
    output.with_suffix(".stderr").write_text(result.stderr, encoding="utf-8")
    return result


def read_output(path):
    with path.open(encoding="utf-8") as stream:
        reader = csv.reader(line for line in stream if line.strip() and not line.lstrip().startswith("#"))
        header = tuple(item.strip() for item in next(reader, ()))
        if header != COLUMNS:
            raise ValueError(f"{path.name}: incorrect CSV columns {header!r}")
        rows = []
        for index, row in enumerate(reader, 1):
            if len(row) != len(COLUMNS):
                raise ValueError(f"{path.name}: row {index} has {len(row)} columns")
            values = tuple(float(value.replace("D", "E").replace("d", "e")) for value in row)
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path.name}: nonfinite value in row {index}")
            rows.append(values)
    return rows


def check_rows(path, orientation, flux_sign, iota_sign):
    rows = read_output(path)
    if len(rows) != len(POINTS):
        raise ValueError(f"{path.name}: expected {len(POINTS)} rows, got {len(rows)}")
    largest = 0.0
    for index, (point, actual) in enumerate(zip(POINTS, rows), 1):
        reference = expected(point, orientation, flux_sign, iota_sign)
        for column, value, target in zip(COLUMNS, actual, reference):
            # All field references are nonzero; coordinates also include exact zero.
            error = abs(value - target) / (abs(target) if target != 0 else 1.0)
            largest = max(largest, error)
            if error > TOLERANCE:
                raise ValueError(f"{path.name}: row {index} {column}: got {value:.17e}, "
                                 f"expected {target:.17e}, relative error {error:.3e}")
    return largest


def check_wrapper(probe, wrapper, bc, directory, timeout, report):
    points, output = directory / "wrapper_points.csv", directory / "wrapper_output.csv"
    with points.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(COLUMNS[:3])
        writer.writerows(POINTS)
    output.unlink(missing_ok=True)
    metadata = output.with_suffix(".csv.json")
    metadata.unlink(missing_ok=True)
    command = [sys.executable, str(wrapper), "--probe", str(probe), "--bc", str(bc),
               "--points", str(points), "--output", str(output), "--inp-swi", "9"]
    result = subprocess.run(command, cwd=directory, text=True, capture_output=True,
                            timeout=timeout, check=False)
    output.with_suffix(".stdout").write_text(result.stdout, encoding="utf-8")
    output.with_suffix(".stderr").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise ValueError(f"wrapper exited {result.returncode}; see wrapper_output.stderr")
    error = check_rows(output, -1, -1, -1)
    sidecar = json.loads(metadata.read_text(encoding="utf-8"))
    paths = {"bc": bc, "points": points, "probe": probe, "output": output,
             "query": output.with_suffix(".csv.points.txt")}
    for key, path in paths.items():
        if sidecar.get(key + "_sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError(f"wrapper sidecar {key} hash mismatch")
    if (sidecar.get("bc_path") != str(bc) or sidecar.get("points_path") != str(points)
            or sidecar.get("inp_swi") != 9 or sidecar.get("bfac") != 1
            or sidecar.get("coordinates") != list(COLUMNS[:3])):
        raise ValueError("wrapper sidecar carrier/points/convention metadata mismatch")
    report["wrapper"] = {"path": str(wrapper), "max_relative_error": error,
                         "sidecar": metadata.name, "exact_file_hashes": "PASS"}


def qualify(probe, wrapper, directory, timeout, report):
    points = directory / "points.txt"
    write_points(points, POINTS)
    for orientation in (-1, 1):
        for flux_sign in (-1, 1):
            for iota_sign in (-1, 1):
                name = f"orientation{orientation:+d}_flux{flux_sign:+d}_iota{iota_sign:+d}"
                bc, output = directory / (name + ".bc"), directory / (name + ".csv")
                write_fixture(bc, orientation, flux_sign, iota_sign)
                result = run_probe(probe, bc, points, output, timeout)
                case = {"fixture": bc.name, "returncode": result.returncode,
                        "orientation": orientation, "flux_sign": flux_sign, "iota_sign": iota_sign}
                report["cases"].append(case)
                if result.returncode != 0:
                    raise ValueError(f"{name}: native probe exited {result.returncode}; see stderr")
                case["rows"] = len(POINTS)
                case["max_relative_error"] = check_rows(output, orientation, flux_sign, iota_sign)
    # The wrapper must reject radial extrapolation, although native magfie clamps it.
    bc = directory / report["cases"][0]["fixture"]
    for label, invalid_s in (("below", 0.05), ("above", 0.95)):
        bad_points = directory / f"outside_{label}.txt"
        write_points(bad_points, [POINTS[0], (invalid_s, 0.73, 0.37)])
        result = run_probe(probe, bc, bad_points, directory / f"outside_{label}.csv", timeout)
        report["rejections"].append({"side": label, "s_tor": invalid_s, "returncode": result.returncode})
        if result.returncode <= 0:
            raise ValueError(f"outside_{label}: expected a normal nonzero rejection, got {result.returncode}")
    if wrapper:
        check_wrapper(probe, wrapper, bc, directory, timeout, report)
        check_wrapper_preserves_inputs(probe, wrapper, bc, directory, timeout)


def check_wrapper_preserves_inputs(probe, wrapper, bc, directory, timeout):
    carrier = directory / "protected.csv.json"
    carrier.write_bytes(bc.read_bytes())
    before = carrier.read_bytes()
    points = directory / "points.txt"
    command = [sys.executable, str(wrapper), "--probe", str(probe), "--bc", str(carrier),
               "--points", str(points), "--output", str(directory / "protected.csv")]
    result = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    if result.returncode == 0 or carrier.read_bytes() != before:
        raise ValueError("wrapper overwrote its input carrier with the metadata sidecar")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path, help="native neo_rt_field_probe.x executable")
    parser.add_argument("--wrapper", type=Path, help="also qualify run_field_probe.py and exact-file sidecar")
    parser.add_argument("--output-dir", type=Path, help="directory retaining fixtures and evidence")
    parser.add_argument("--timeout", type=float, default=30.0, help="seconds allowed per native invocation")
    args = parser.parse_args()
    probe = args.probe.expanduser().resolve()
    if not probe.is_file():
        parser.error(f"probe does not exist: {probe}")
    wrapper = args.wrapper.expanduser().resolve() if args.wrapper else None
    if wrapper and not wrapper.is_file():
        parser.error(f"wrapper does not exist: {wrapper}")
    if args.timeout <= 0 or not math.isfinite(args.timeout):
        parser.error("timeout must be finite and positive")
    directory = args.output_dir.expanduser().resolve() if args.output_dir else Path(tempfile.mkdtemp(prefix="neort-probe-"))
    directory.mkdir(parents=True, exist_ok=True)
    report = {"status": "FAIL", "probe": str(probe), "input_switch": 9,
              "relative_tolerance": TOLERANCE, "radial_surfaces": SURFACES,
              "cases": [], "rejections": [],
              "oracle": "B[T]=2+0.1*s*cos(theta)+0.03*s*sin(theta); "
                        "Bcov[G cm]=-0.2*orientation*I[A]; psi_pr=flux[Wb]*1e8/(2*pi); "
                        "sqrtg=-psi_pr*(Bphcov+iota*Bthcov)/Bmod^2"}
    try:
        qualify(probe, wrapper, directory, args.timeout, report)
        report["status"] = "PASS"
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        report["error"] = str(error)
    (directory / "qualification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if report["status"] != "PASS":
        print(f"FAIL: {report['error']}\nEvidence: {directory}", file=sys.stderr)
        return 1
    maximum = max(case["max_relative_error"] for case in report["cases"])
    print(f"PASS: {len(report['cases'])} analytic cases, {len(POINTS)} points each; "
          f"max relative error {maximum:.3e}; both radial bounds rejected.\nEvidence: {directory}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
