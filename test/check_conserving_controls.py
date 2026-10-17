"""Native conserving input rejections and disabled-mode behavioral parity."""
from pathlib import Path
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np

executable, config_executable, field, scratch = map(Path, sys.argv[1:5])
scratch.mkdir(parents=True, exist_ok=True)
shutil.copyfile(field, scratch / "in_file")
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"}
for mode, message in (
    ("rate", "conserving_rate_scale must be finite"),
    ("speed", "steps must be at least 2"), ("pitch", "steps must be at least 2"),
    ("umax", "conserving_umax must be finite"),
    ("quadrature", "quadrature orders must be at least 2"),
    ("azimuth", "pair azimuth order must be at least 4"),
    ("rank", "conserving_rank_rtol must be finite"),
    ("nan-rank", "conserving_rank_rtol must be finite"),
    ("compatibility", "conserving_compatibility_rtol must be finite"),
    ("tp", "select only one of TP FP"),
    ("lorentz", "select only one of TP FP"),
    ("krook", "select only one of TP FP"),
    ("drift-phase", "linear common-phase thin orbits"),
    ("fow", "linear common-phase thin orbits"),
    ("nopassing", "linear common-phase thin orbits"),
):
    result = subprocess.run([str(config_executable.resolve()), mode], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0 and message in result.stdout + result.stderr, (
        mode, result.stdout, result.stderr)

species = """&species_count n_species=3 /
&species reference_s=.5,omega_e=-700,domega_e_ds=0,coulomb_log=20,
species_name='D','T','e',charge_e=1,1,-1,
mass_kg=3.345228e-27,5.017842e-27,9.1091e-31,
density_m3=5e15,5e15,1e16,temperature_ev=2000,2000,2000,
density_gradient_m3=0,0,0,temperature_gradient_ev=0,0,0 /
"""
(scratch / "species.in").write_text(species)
legacy = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,supban=T,nopassing=T,noshear=T,inp_swi=8,
vsteps=1,vmax_over_vth=1,mth_max_abs=0,log_level=-1,output_format='both',
krook_collisions=T,krook_nu=1000,krook_eta_steps=3,krook_rtol=1e-4 {extra} /
"""
missing = "&params conserving_fp_collisions=T,qs=1,ms=2 /\n"
(scratch / "missing_mixture.in").write_text(missing)
result = subprocess.run([str(executable.resolve()), "missing_mixture"], cwd=scratch,
                        env=environment, capture_output=True, text=True, timeout=10)
assert result.returncode != 0 and "requires an explicit species mixture" in (
    result.stdout + result.stderr), (result.stdout, result.stderr)

for name, extra in (("omitted", ""), ("disabled", ",conserving_fp_collisions=F")):
    (scratch / f"{name}.in").write_text(legacy.format(extra=extra))
    result = subprocess.run([str(executable.resolve()), name], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, (name, result.stdout[-3000:], result.stderr[-3000:])
    assert not (scratch / f"{name}_D_collision_model.out").exists()
for path in scratch.glob("omitted*.out"):
    assert path.read_bytes() == path.with_name(path.name.replace("omitted", "disabled", 1)).read_bytes()
with h5py.File(scratch / "omitted.h5", "r") as old, \
        h5py.File(scratch / "disabled.h5", "r") as new:
    def compare(name, obj):
        if name == "created":
            return
        assert name in new, name
        if isinstance(obj, h5py.Dataset):
            np.testing.assert_array_equal(new[name][()], obj[()], err_msg=name)
        for key, value in obj.attrs.items():
            np.testing.assert_array_equal(new[name].attrs[key], value)
    old.visititems(compare)
    assert set(old.keys()) == set(new.keys())
    for name in ("D", "T", "e"):
        assert "conserving_fp_collisions" not in new[f"species/{name}/config"]
        assert "source_feature" not in new[f"species/{name}/transport"]
print("Native conserving controls reject unsupported inputs; disabled output parity passes")
