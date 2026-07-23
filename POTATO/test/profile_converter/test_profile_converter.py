#!/usr/bin/env python3
"""Regression for explicit NEO-RT-to-POTATO profile and frequency signs."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


TOOL = Path(__file__).parents[2] / "tools" / "neo_rt_profiles_to_potato.py"
SPEC = importlib.util.spec_from_file_location("neo_rt_profiles_to_potato", TOOL)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        profile = root / "profile.in"
        plasma = root / "plasma.in"
        geometry = root / "components.npz"
        profile.write_text(
            "0.0 0.02\n"
            "0.5 0.03\n"
            "1.0 0.04\n"
        )
        plasma_grid = np.linspace(0.0, 1.0, 239)
        plasma.write_text(
            "% N am1 am2 Z1 Z2\n"
            f"{plasma_grid.size} 2.0 3.0 1.0 1.0\n"
            "% s n1 n2 T1 T2 Te\n"
            + "".join(
                f"{s:.16e} {(1.0-0.4*s)*1.0e14:.16e} 0.0 "
                f"{(1.0-0.4*s)*1.0e4:.16e} 1.0 "
                f"{(1.2-0.4*s)*1.0e4:.16e}\n"
                for s in plasma_grid
            )
        )
        s = np.linspace(0.0, 1.0, 11)
        np.savez(geometry, s_sqrt_poloidal=np.sqrt(s), s_toroidal=s)

        metadata = {}
        for sign in (-1, 0, 1):
            output = root / f"profile_{sign:+d}.in"
            metadata[sign] = MODULE.convert(
                profile, plasma, geometry, output, r0_cm=100.0,
                psi_span_tm2=2.0, relation_sign=sign, degree=3,
            )
            if metadata[sign]["selected_ion_index_one_based"] != 1:
                raise AssertionError("physical ion was not selected by charge and density")
            if metadata[sign]["plasma_grid_count"] != 239:
                raise AssertionError("non-default NEO-RT grid count was not preserved")

        plus = metadata[1]["potential_statv_edge"]
        minus = metadata[-1]["potential_statv_edge"]
        zero = metadata[0]["potential_statv_edge"]
        if not plus > 0.0 or not np.isclose(minus, -plus) or zero != 0.0:
            raise AssertionError("potential-gradient relation signs were not preserved")
        zero_coefficients = np.asarray(
            MODULE.numeric_rows(root / "profile_+0.in")[3]
        )
        if not np.array_equal(zero_coefficients, np.zeros(4)):
            raise AssertionError("zero-frequency control produced a nonzero potential")

        profile_with_vth = root / "profile_with_vth.in"
        temperatures = (1.0-0.4*plasma_grid)*1.0e4
        vth = np.sqrt(
            2.0*temperatures*MODULE.EV_CGS/(2.0*MODULE.AMU_CGS)
        )
        mach = np.interp(plasma_grid, [0.0, 0.5, 1.0], [0.02, 0.03, 0.04])
        profile_with_vth.write_text(
            "".join(
                f"{s:.16e} {m:.16e} {speed:.16e}\n"
                for s, m, speed in zip(plasma_grid, mach, vth, strict=True)
            )
        )
        with_vth = MODULE.convert(
            profile_with_vth, plasma, geometry, root / "profile_with_vth.out",
            r0_cm=100.0, psi_span_tm2=2.0, relation_sign=1, degree=3,
        )
        if with_vth["optional_profile_vth_max_relative_error"] > 1.0e-5:
            raise AssertionError("equivalent optional vth column did not close")
        if not np.isclose(
            with_vth["potential_statv_edge"],
            metadata[1]["potential_statv_edge"],
            rtol=1.0e-12,
        ):
            raise AssertionError("two- and three-column profile schemas differ")

        truncated = root / "plasma_truncated.in"
        truncated.write_text("\n".join(plasma.read_text().splitlines()[:-1]) + "\n")
        try:
            MODULE.convert(
                profile, truncated, geometry, root / "invalid.in",
                r0_cm=100.0, psi_span_tm2=2.0, relation_sign=1, degree=3,
            )
        except ValueError as error:
            if "declares 239 rows, found 238" not in str(error):
                raise
        else:
            raise AssertionError("truncated NEO-RT plasma grid was accepted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
