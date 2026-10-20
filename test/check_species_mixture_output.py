"""Read the native mixture file and compare physically additive torque outputs."""
from pathlib import Path
import sys

import h5py
import numpy as np


prefix = Path(sys.argv[1])
text = np.loadtxt(str(prefix) + "_species.out", dtype=str)
with h5py.File(str(prefix) + ".h5", "r") as output:
    assert int(np.asarray(output["format_version"]).squeeze()) == 2
    assert int(np.asarray(output["n_species"]).squeeze()) == 3
    assert set(output["species"]) == {"D", "T", "e"}
    summed = np.zeros(3)
    for row in text:
        name, *values = row
        values = np.asarray(values, dtype=float)
        species = output[f"species/{name}"]
        np.testing.assert_allclose(species["mass_kg"][()], values[0], rtol=1e-14)
        np.testing.assert_allclose(species["charge_e"][()], values[1], rtol=1e-14)
        np.testing.assert_allclose(species["magfie/Om_tE"][()], -700.0, rtol=1e-14)
        np.testing.assert_allclose(species["magfie/T"][()], 2000.0, rtol=1e-14)
        np.testing.assert_allclose(species["config/qs"][()], values[1], rtol=1e-14)
        torque = np.asarray([species[f"torque/{key}"][()].item()
                             for key in ("Tco", "Tctr", "Tt")])
        np.testing.assert_allclose(torque, values[9:12], rtol=1e-14)
        summed += torque
    np.testing.assert_allclose(output["mixture/torque_by_population"][()], summed, rtol=1e-14)
    np.testing.assert_allclose(output["mixture/torque_total"][()], summed.sum(), rtol=1e-14)
    total_text = np.loadtxt(str(prefix) + "_mixture_torque.out")
    np.testing.assert_allclose(total_text[1:4], summed, rtol=1e-14)
    np.testing.assert_allclose(total_text[4], summed.sum(), rtol=1e-14)
print("Native HDF5 and text D/T/e species metadata and torque sums agree")
