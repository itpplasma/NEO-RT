"""Native D/T/e dispatch consumes the certified noncosine well in every mode."""
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
(scratch / "species.in").write_text("""&species_count n_species=3 /
&species reference_s=.5,omega_e=-700,domega_e_ds=0,coulomb_log=20,
species_name='D','T','e',charge_e=1,1,-1,
mass_kg=3.3435837724e-27,5.0073567446e-27,9.1093837139e-31,
density_m3=5e15,5e15,1e16,temperature_ev=2000,2000,2000,
density_gradient_m3=0,0,0,temperature_gradient_ev=0,0,0 /
""")
base = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,noshear=T,inp_swi=8,vsteps=1,vmax_over_vth=2,
log_level=-1,output_format='both', {controls} /
"""
models = {
    "lorentz": "lorentz_pitch_collisions=T,lorentz_rate_source=1,lorentz_pitch_steps=2",
    "tp": "tp_fp_collisions=T,fp_speed_steps=2,fp_pitch_steps=2,fp_umax=4",
    "landau": """conserving_fp_collisions=T,conserving_speed_steps=2,
conserving_pitch_steps=2,conserving_theta_order=4,conserving_velocity_order=2,
conserving_cm_order=2,conserving_rel_order=2,conserving_polar_order=2,
conserving_azimuth_order=4""",
}
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"}
original = (scratch / "in_file").read_text()
epsilon = .5 * np.sqrt(.5) / 5
expected_min = 2e4 * (1 - epsilon + .1 * epsilon)
expected_max = 2e4 * (1 + epsilon + .1 * epsilon)


def run(name, controls, expected_error=None):
    (scratch / f"{name}.in").write_text(base.format(controls=controls))
    result = subprocess.run([str(executable.resolve()), name], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=45)
    output = result.stdout + result.stderr
    if expected_error is not None:
        assert result.returncode != 0 and expected_error in output, output[-3000:]
        return
    assert result.returncode == 0, output[-3000:]
    with h5py.File(scratch / f"{name}.h5") as data:
        assert set(data["species"]) == {"D", "T", "e"}
        for label in ("D", "T", "e"):
            species = data[f"species/{label}"]
            np.testing.assert_allclose(1 / species["magfie/etadt"][()].item(),
                                       expected_min, rtol=2e-14)
            np.testing.assert_allclose(1 / species["magfie/etatp"][()].item(),
                                       expected_max, rtol=2e-14)
            torques = [species[f"torque/{key}"][()].item()
                       for key in ("Tco", "Tctr", "Tt")]
            assert np.isfinite(torques).all() and np.linalg.norm(torques) > 0


for name, controls in models.items():
    run(name, controls)

# This negative second harmonic produces extra maxima while leaving the
# minimum at theta=0. It must fail the strict analytic curvature certificate.
lines = []
for line in original.splitlines():
    numbers = line.split()
    if len(numbers) == 6 and numbers[:2] == ["2", "0"]:
        numbers[-1] = "-.04"
        line = " ".join(numbers)
    lines.append(line)
(scratch / "in_file").write_text("\n".join(lines) + "\n")
for name, controls in models.items():
    run(f"extra_wells_{name}", controls, "status=-14")
(scratch / "in_file").write_text(original)
print("Native two-harmonic Lorentz/TP/Landau species readback and multiwell rejection pass")
