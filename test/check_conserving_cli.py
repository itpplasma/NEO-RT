"""Actual native common-T Landau CLI: physical EM work and source receipts."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np

executable, field, scratch = map(Path, sys.argv[1:4])
scratch.mkdir(parents=True, exist_ok=True)
shutil.copyfile(field, scratch / "in_file")
species_input = """&species_count n_species=3 /
&species reference_s=.5,omega_e=-700,domega_e_ds=0,coulomb_log=20,
species_name='D','T','e',charge_e=1,1,-1,
mass_kg=3.3435837724e-27,5.0073567446e-27,9.1093837139e-31,
density_m3=5e18,5e18,1e19,temperature_ev=2000,2000,2000,
density_gradient_m3=-2e18,-2e18,-4e18,
temperature_gradient_ev=400,400,400 /
"""
deck = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,noshear=T,inp_swi=8,log_level=-1,
output_format='both',conserving_fp_collisions=T,conserving_speed_steps=2,
conserving_pitch_steps=2,conserving_theta_order=8,conserving_velocity_order=4,
conserving_cm_order=2,conserving_rel_order=2,conserving_polar_order=3,
conserving_azimuth_order=6,conserving_rate_scale=.6 /
"""
(scratch / "species.in").write_text(species_input)
(scratch / "conserving.in").write_text(deck)
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"}


def run(name, error=None):
    result = subprocess.run([str(executable.resolve()), name], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=45)
    output = result.stdout + result.stderr
    if error is None:
        assert result.returncode == 0, output[-4000:]
    else:
        assert result.returncode != 0 and error in output, output[-4000:]


def scalar(group, path):
    return group[path][()].item()


run("conserving")
with h5py.File(scratch / "conserving.h5", "r") as output:
    assert "conserving_matrix_receipt" not in output["mixture"]
(scratch / "receipt.in").write_text(deck.replace(
    "conserving_rate_scale=.6", "conserving_rate_scale=.6,conserving_matrix_receipt=T"))
run("receipt")
from check_conserving_receipt import check_receipt
check_receipt(scratch / "conserving.h5", scratch / "receipt.h5")
total = np.zeros(3)
elementary_charge = 1.602176634e-19
with h5py.File(scratch / "conserving.h5", "r") as output:
    assert set(output["species"]) == {"D", "T", "e"}
    for name in ("D", "T", "e"):
        species = output[f"species/{name}"]
        config, magfie, transport = (species[key] for key in
                                     ("config", "magfie", "transport"))
        assert scalar(config, "conserving_fp_collisions")
        np.testing.assert_equal(scalar(config, "conserving_rate_scale"), .6)
        assert scalar(transport, "conserving_active")
        assert not scalar(transport, "susceptibility_valid")
        assert scalar(transport, "source_feature").decode().strip() == (
            "all-pair-linearized-Landau-commonT-v1")
        assert "infinity" in scalar(transport, "conserving_velocity_domain").decode()
        assert scalar(transport, "conserving_rank") == 51
        assert scalar(transport, "conserving_nullity") == 0
        assert scalar(transport, "conserving_residual") < 1e-10
        assert scalar(transport, "conserving_entropy") >= 0
        energy = 2000 * elementary_charge
        np.testing.assert_allclose(scalar(magfie, "fp_thermal_energy_J"), energy,
                                   rtol=2e-14)
        np.testing.assert_allclose(scalar(magfie, "fp_vth_m_s"),
                                   np.sqrt(2 * energy / scalar(species, "mass_kg")),
                                   rtol=2e-14)
        work = (transport["conserving_work_real"][()].ravel() +
                1j * transport["conserving_work_imag"][()].ravel())
        cross = (transport["conserving_cross_work_real"][()] +
                 1j * transport["conserving_cross_work_imag"][()])
        np.testing.assert_array_equal(
            transport["conserving_cross_work_qualified"][()], 1)
        np.testing.assert_allclose(cross.sum(axis=0), work, rtol=2e-10)
        expected = -3 * np.pi * work.imag
        actual = np.array([scalar(species, f"torque/{key}")
                           for key in ("Tt", "Tco", "Tctr")])
        assert np.all(np.isfinite(actual)) and np.linalg.norm(actual) > 0
        np.testing.assert_allclose(actual, expected, rtol=3e-13)
        total += actual[[1, 2, 0]]
        for key in ("D11", "D12"):
            np.testing.assert_array_equal(transport[key][()], 0)
        metadata = (scratch / f"conserving_{name}_collision_model.out").read_text()
        assert "all-pair-linearized-Landau-commonT-v1" in metadata
        assert "infinity" in metadata
    np.testing.assert_allclose(output["mixture/torque_by_population"][()].ravel(),
                               total, rtol=3e-13)

# A physical source along a conserved density moment cannot have a stationary
# zero-detuning solution. Reject it instead of adding a collision/rate floor.
(scratch / "species.in").write_text(species_input.replace("omega_e=-700", "omega_e=0"))
(scratch / "incompatible.in").write_text(deck)
run("incompatible", "incompatible or unresolved conserving FP response")
print("Native coupled common-T Landau CLI work and source-compatibility checks passed")
