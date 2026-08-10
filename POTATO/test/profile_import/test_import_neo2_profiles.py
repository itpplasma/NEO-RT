#!/usr/bin/env python3
"""Check the NEO-2 HDF5 profile importer uses POTATO's flux gauge."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np


TOOL = Path(__file__).parents[2] / "python" / "import_neo2_profiles.py"
SPEC = importlib.util.spec_from_file_location("import_neo2_profiles", TOOL)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def main() -> int:
    surfaces = np.linspace(0.0, 1.0, 11)
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        neo2 = root / "neo2_out.h5"
        output = root / "profile_poly.in"
        with h5py.File(neo2, "w") as handle:
            handle["boozer_s"] = surfaces
            handle["aiota"] = np.ones_like(surfaces)
            handle["n_spec"] = np.ones((surfaces.size, 2))
            handle["T_spec"] = np.full((surfaces.size, 2), 2.0)
            handle["m_spec"] = np.full((surfaces.size, 2), 2.0)
            handle["MtOvR"] = np.full((surfaces.size, 2), 3.0)
            handle["Bref"] = np.ones_like(surfaces)
            handle["psi_pr_hat"] = np.full_like(surfaces, 5.0)
            handle["R0"] = np.array([100.0])

        class FakeConverter:
            def __init__(self, qprof, s_pol):
                assert np.array_equal(qprof, [1.0])
                assert np.array_equal(s_pol, [0.0, 1.0])

            def stor2spol(self, values):
                return np.asarray(values)**2

        MODULE.FluxConverter = FakeConverter
        MODULE.eqdsk_base.read_eqdsk = lambda _: {
            "qprof": np.array([1.0]),
            "s_pol": np.array([0.0, 1.0]),
            "PsiedgeVs": 2.0,
            "PsiaxisVs": 0.0,
        }

        mapped = MODULE.stor_to_spol(np.array([0.25, 0.5]), root / "fake.eqdsk")
        if not np.allclose(mapped, [0.0625, 0.25]):
            raise AssertionError("s_tor to s_pol conversion was not used")

        original_argv = sys.argv
        try:
            sys.argv = [
                str(TOOL), "--neo2", str(neo2), "--eqdsk", str(root / "fake.eqdsk"),
                "--output", str(output),
            ]
            MODULE.main()
        finally:
            sys.argv = original_argv

        rows = [
            [float(value) for value in line.split()]
            for line in output.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        phi_coefficients = np.asarray(rows[3])
        # The test has constant Er, so phi_e is linear.  Its magnitude must use
        # the EQDSK poloidal span (2e8 Mx), not psi_pr_hat*Bref (5 Mx here).
        er = -2.0e8 * np.sqrt(2.0) * 3.0 / MODULE.C_CGS
        expected_phi = -er * surfaces**2
        if not np.allclose(np.polyval(phi_coefficients, surfaces**2), expected_phi,
                           rtol=1.0e-10, atol=1.0e-12):
            raise AssertionError("profile potential did not use the EQDSK flux gauge")

    return 0


if __name__ == "__main__":
    sys.exit(main())
