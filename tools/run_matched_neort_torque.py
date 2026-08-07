#!/usr/bin/env python3
"""Run radial NEO-RT torque scans for the matched circular comparison."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import shutil
import subprocess
import time
from pathlib import Path

import netCDF4
import numpy as np


def radial_map(chartmap: Path) -> tuple[np.ndarray, np.ndarray]:
    with netCDF4.Dataset(chartmap) as data:
        s_tor = np.asarray(data.variables["s"][:], dtype=float)
        a_phi = np.asarray(data.variables["A_phi"][:], dtype=float)
    s_pol = (a_phi - a_phi[0]) / (a_phi[-1] - a_phi[0])
    if np.any(np.diff(s_pol) <= 0.0):
        raise RuntimeError("chartmap poloidal-flux map is not monotone")
    return s_tor, s_pol


def write_profiles(directory: Path, s_tor: np.ndarray, s_pol: np.ndarray) -> None:
    density = 5.0e13 - 2.5e13 * s_pol
    temperature = 2.0e3 - 1.4e3 * s_pol
    mass_kg = 2.0 * 1.66053906660e-27
    speed_cm_s = np.sqrt(2.0 * temperature * 1.602176634e-19 / mass_kg) * 100.0
    with (directory / "plasma.in").open("w") as stream:
        stream.write("% N am1 am2 Z1 Z2\n")
        stream.write(f"{len(s_tor)} 2.0 2.0 1.0 1.0\n")
        stream.write("% s ni_1 ni_2 Ti_1 Ti_2 Te\n")
        for values in zip(s_tor, density, temperature, strict=True):
            surface, number_density, ion_temperature = values
            stream.write(
                f"{surface:.16e} {number_density:.16e} 0.0 "
                f"{ion_temperature:.16e} 1.0 {ion_temperature:.16e}\n"
            )
    with (directory / "profile.in").open("w") as stream:
        for surface, speed in zip(s_tor, speed_cm_s, strict=True):
            stream.write(f"{surface:.16e} 0.0 {speed:.16e}\n")


def deck(surface: float, model: int, vsteps: int, width_scale: float) -> str:
    input_switch = 10 if model == 0 else 11
    return f"""&params
    s = {surface:.16e}
    M_t = 0.0
    qs = 1.0
    ms = 2.0
    epsmn = 1.0e-3
    m0 = 0
    mph = 2
    supban = .false.
    magdrift = .true.
    magdrift_passing = 1
    nopassing = .false.
    noshear = .true.
    pertfile = .false.
    nonlin = .false.
    comptorque = .true.
    bfac = 1.0
    efac = 1.0
    inp_swi = {input_switch}
    frequency_model = {model}
    gc_full_orbit_width_scale = {width_scale:.16e}
    mth_max_abs = 3
    vsteps = {vsteps}
    vmax_over_vth = 4.0
    log_level = 2
/
"""


def read_torque(path: Path) -> tuple[float, ...]:
    values = np.loadtxt(path, comments="#")
    values = np.atleast_2d(values)
    if values.shape != (1, 6) or not np.all(np.isfinite(values)):
        raise RuntimeError(f"invalid torque output {path}: shape={values.shape}")
    return tuple(float(value) for value in values[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--chartmap", type=Path, required=True)
    parser.add_argument("--eqdsk", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument("--width-scale", type=float, default=1.0)
    parser.add_argument("--vsteps", type=int, default=128)
    parser.add_argument("--surface-count", type=int, default=25)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.width_scale < 0.0:
        parser.error("--width-scale must be nonnegative")

    args.output.mkdir(parents=True, exist_ok=False)
    map_s_tor, map_s_pol = radial_map(args.chartmap)
    target_s_pol = np.linspace(0.04, 0.64, args.surface_count)
    target_s_tor = np.interp(target_s_pol, map_s_pol, map_s_tor)
    write_profiles(args.output, map_s_tor, map_s_pol)

    start_all = time.perf_counter()

    def run_surface(item: tuple[int, float, float]) -> dict[str, float]:
        index, surface, poloidal = item
        work = args.output / f"surface-{index:03d}"
        work.mkdir()
        shutil.copy2(args.chartmap if args.model == 0 else args.eqdsk, work / "in_file")
        shutil.copy2(args.output / "plasma.in", work / "plasma.in")
        shutil.copy2(args.output / "profile.in", work / "profile.in")
        (work / "torque.in").write_text(
            deck(float(surface), args.model, args.vsteps, args.width_scale)
        )
        started = time.perf_counter()
        completed = subprocess.run(
            [str(args.executable.resolve()), "torque"],
            cwd=work,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env={
                **os.environ,
                "OMP_NUM_THREADS": "1",
                "OMP_DYNAMIC": "false",
                "OMP_PROC_BIND": "false",
            },
            check=False,
        )
        (work / "run.log").write_text(completed.stdout)
        if completed.returncode != 0:
            raise RuntimeError(f"surface {surface} failed; see {work / 'run.log'}")
        native_s, dvds, mach, co, counter, trapped = read_torque(
            work / "torque_torque.out"
        )
        if not math.isclose(native_s, surface, rel_tol=0.0, abs_tol=1.0e-12):
            raise RuntimeError("NEO-RT returned a different radial surface")
        total = co + counter + trapped
        density_si = total * 1.0e-7 / (dvds * 1.0e-6)
        return {
            "s_pol": poloidal,
            "s_tor": surface,
            "rho_tor": math.sqrt(surface),
            "dVds_cm3": dvds,
            "mach": mach,
            "Tco_dyn_cm": co,
            "Tcounter_dyn_cm": counter,
            "Ttrapped_dyn_cm": trapped,
            "Ttotal_dyn_cm": total,
            "torque_density_Nm_m3": density_si,
            "wall_seconds": time.perf_counter() - started,
        }

    items = [
        (index, float(surface), float(poloidal))
        for index, (surface, poloidal) in enumerate(
            zip(target_s_tor, target_s_pol, strict=True)
        )
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(run_surface, items))

    with (args.output / "torque.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (args.output / "wall_seconds.txt").write_text(
        f"{time.perf_counter() - start_all:.9f}\n"
    )


if __name__ == "__main__":
    main()
