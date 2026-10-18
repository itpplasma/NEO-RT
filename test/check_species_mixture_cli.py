"""One native invocation consumes an explicit D/T/e input without legacy profiles."""
from pathlib import Path
import os
import subprocess
import sys

import h5py
import numpy as np


executable, field, scratch = map(Path, sys.argv[1:])
scratch.mkdir(parents=True, exist_ok=True)
field_link = scratch / "in_file"
if not field_link.exists():
    field_link.symlink_to(field.resolve())
(scratch / "mixture.in").write_text("""&params
 s=.09, species_file='species.in', epsmn=.001, m0=0, mph=3,
 comptorque=.true., supban=.true., magdrift=.false., nopassing=.true.,
 nonlin=.false., inp_swi=8, vsteps=2, mth_max_abs=0, vmax_over_vth=1,
 log_level=-1, krook_collisions=.true., krook_nu=1000,
 krook_eta_steps=4, krook_rtol=1e-4, output_format='both'
/
""")
(scratch / "species.in").write_text("""&species_count n_species=3 /
&species
 reference_s=.09, omega_e=-700, domega_e_ds=0, coulomb_log=20,
 species_name='D','T','e', charge_e=1,1,-1,
 mass_kg=3.3435837724e-27,5.0073567446e-27,9.1093837139e-31,
 density_m3=5e15,5e15,1e16, temperature_ev=2000,2000,2000,
 density_gradient_m3=0,0,0, temperature_gradient_ev=0,0,0
/
""")
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OPENBLAS_NUM_THREADS": "1"}
completed = subprocess.run([str(executable.resolve()), "mixture"], cwd=scratch,
                           env=environment, capture_output=True, text=True, timeout=45)
if completed.returncode:
    raise AssertionError(completed.stdout[-4000:] + completed.stderr[-4000:])
with h5py.File(scratch / "mixture.h5", "r") as output:
    assert set(output["species"]) == {"D", "T", "e"}
    for name in ("D", "T", "e"):
        species = output[f"species/{name}"]
        np.testing.assert_allclose(species["magfie/Om_tE"][()], -700.0, rtol=1e-14)
        mach = species["magfie/M_t"][()].item()
        speed = species["magfie/vth"][()].item()
        radius = species["magfie/R0"][()].item()
        np.testing.assert_allclose(mach * speed / radius, -700.0, rtol=1e-14)
        assert np.isfinite(species["torque/Tt"][()]).all()
        assert (scratch / f"mixture_{name}_torque.out").is_file()
    np.testing.assert_allclose(output["mixture/coulomb_log"][()], 20.0, rtol=1e-14)
print("One native invocation returns D/T/e with shared physical rotation")
