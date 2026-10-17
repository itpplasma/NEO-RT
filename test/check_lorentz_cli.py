"""Exercise native Lorentz input, common-phase dispatch and serialized outputs."""
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
mass_kg=3.345228e-27,5.017842e-27,9.1091e-31,
density_m3=5e15,5e15,1e16,temperature_ev=1800,2600,2100,
density_gradient_m3=0,0,0,temperature_gradient_ev=0,0,0 /
""")
base = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,noshear=T,inp_swi=8,vsteps=1,vmax_over_vth=2,
mth_max_abs=0,log_level=-1,output_format='both',
lorentz_pitch_collisions=T,lorentz_rate_source=1,lorentz_rate_scale=.7,
lorentz_pitch_steps=3 {extra} /
"""
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"}


def run(name, deck, expected_error=None, binary=executable):
    (scratch / f"{name}.in").write_text(deck)
    result = subprocess.run([str(binary.resolve()), name], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=45)
    output = result.stdout + result.stderr
    if expected_error is None:
        assert result.returncode == 0, output[-4000:]
    else:
        assert result.returncode != 0, f"{name}: invalid input was accepted"
        assert expected_error in output, f"{name}: unexpected rejection {output[-2000:]}"


run("lorentz", base.format(extra=""))
population_total = np.zeros(3)
with h5py.File(scratch / "lorentz.h5", "r") as output:
    assert set(output["species"]) == {"D", "T", "e"}
    for name in ("D", "T", "e"):
        species = output[f"species/{name}"]
        config = species["config"]
        assert config["lorentz_pitch_collisions"][()].item()
        assert config["lorentz_rate_source"][()].item() == 1
        np.testing.assert_equal(config["lorentz_rate_scale"][()].item(), .7)
        assert config["lorentz_pitch_steps"][()].item() == 3
        assert b"nonadiabatic_g" in config["lorentz_relaxed_variable"][()].tobytes()
        assert b"common_phase_bounce_averaged" in config["lorentz_response_model"][()].tobytes()
        np.testing.assert_allclose(species["magfie/Om_tE"][()].item(), -700, rtol=2e-14)
        labels = np.asarray(species["harmonics/mth"]).ravel()
        np.testing.assert_array_equal(labels, [0, -6])
        torque = []
        for component in ("Tco", "Tctr", "Tt"):
            values = np.asarray(species[f"harmonics/{component}"]).ravel()
            total = species[f"torque/{component}"][()].item()
            assert np.isfinite(total) and total != 0
            np.testing.assert_allclose(values.sum(), total, rtol=2e-14)
            assert values[0 if component != "Tt" else 1] == 0
            torque.append(total)
        population_total += torque
    assert output["species/e/charge_e"][()].item() == -1
    np.testing.assert_allclose(np.asarray(output["mixture/torque_by_population"]).ravel(),
                               population_total, rtol=2e-14)
    assert b"Lorentz" in np.asarray(
        output["mixture/torque_observable"].attrs["comment"]).astype("S").tobytes()
text = np.loadtxt(scratch / "lorentz_mixture_torque.out")
np.testing.assert_allclose(text[1:4], population_total, rtol=2e-14)

for name, replacement, message in (
    ("negative_rate", ("lorentz_rate_scale=.7", "lorentz_rate_scale=-1"),
     "lorentz_rate_scale must be finite"),
    ("nan_rate", ("lorentz_rate_scale=.7", "lorentz_rate_scale=NaN"),
     "lorentz_rate_scale must be finite"),
    ("bad_source", ("lorentz_rate_source=1", "lorentz_rate_source=2"),
     "lorentz_rate_source must be"),
    ("bad_pitch", ("lorentz_pitch_steps=3", "lorentz_pitch_steps=1"),
     "lorentz_pitch_steps must be"),
    ("zero_rate_nonzero_drive", ("lorentz_rate_scale=.7", "lorentz_rate_scale=0"),
     "status=-17"),
    ("missing_mixture", ("species_file='species.in'", "species_file=''"),
     "all-background Lorentz rates require"),
):
    run(name, base.format(extra="").replace(*replacement), message)
for name, extra, message in (
    ("simultaneous", ",krook_collisions=T,krook_nu=1", "either Lorentz diffusion or Krook"),
    ("drift_phase", ",drift_phase=T", "linear thin orbits and all branches"),
    ("fow", ",fow_radial=T", "linear thin orbits and all branches"),
    ("trapped_only", ",nopassing=T", "linear thin orbits and all branches"),
    ("shaing", ",supban=T", "linear thin orbits and all branches"),
    ("nonlinear", ",nonlin=T", "general mixtures require linear transport"),
):
    run(name, base.format(extra=extra), message)

# The source fixture is analytic: alter its q profile without changing B.
original_field = (scratch / "in_file").read_text()
for mode, status in (("noninteger", -13), ("shear", -12)):
    lines = original_field.splitlines()
    for index, line in enumerate(lines):
        numbers = line.split()
        if len(numbers) == 6 and numbers[0] not in ("CC", "0", "1"):
            radius = float(numbers[0])
            if radius in (.25, .5, .75):
                q = 2.1 if mode == "noninteger" else 2 + .1 * (radius - .5)
                numbers[1] = repr(1 / q)
                lines[index] = " ".join(numbers)
    (scratch / "in_file").write_text("\n".join(lines) + "\n")
    run(mode, base.format(extra=""), f"status={status}")
(scratch / "in_file").write_text(original_field)
# A weak second B harmonic satisfies the analytic single-well certificate.
shaped_lines = []
for line in original_field.splitlines():
    numbers = line.split()
    if len(numbers) == 7 and numbers[0] == "1":
        numbers[0] = "2"
        line = " ".join(numbers)
    shaped_lines.append(line)
    if len(numbers) == 6 and numbers[:2] == ["1", "0"]:
        shaped_lines.append("2 0 0 0 0 .01")
(scratch / "in_file").write_text("\n".join(shaped_lines) + "\n")
run("shaped_well", base.format(extra=""))
(scratch / "in_file").write_text(original_field)
# A small m=33 term aliases m=1 at every node of the old 32-point guard.
# Keep R/Z circular and enough low-harmonic B for the native orbit initializer.
alias_lines = []
for line in original_field.splitlines():
    numbers = line.split()
    if len(numbers) == 7 and numbers[0] == "1":
        numbers[0] = "33"
        line = " ".join(numbers)
    alias_lines.append(line)
    if len(numbers) == 6 and numbers[:2] == ["1", "0"]:
        alias_lines.extend(f"{m} 0 0 0 0 {-1e-5 if m == 33 else 0}"
                           for m in range(2, 34))
(scratch / "in_file").write_text("\n".join(alias_lines) + "\n")
run("aliased_shape", base.format(extra=""), expected_error="status=-14")
(scratch / "in_file").write_text(original_field)
(scratch / "zero_rotation.in").write_text((scratch / "species.in").read_text()
    .replace("omega_e=-700", "omega_e=0"))
for name, rate in (("zero_drive", ".7"), ("zero_rate_drive", "0")):
    run(name, base.format(extra="").replace("species.in", "zero_rotation.in")
        .replace("lorentz_rate_scale=.7", f"lorentz_rate_scale={rate}"))
    with h5py.File(scratch / f"{name}.h5", "r") as output:
        for label in ("D", "T", "e"):
            species = output[f"species/{label}"]
            assert not species["transport/susceptibility_valid"][()].item()
            for component in ("Tco", "Tctr", "Tt"):
                assert species[f"torque/{component}"][()].item() == 0
            assert "susceptibility_valid=F" in (
                scratch / f"{name}_{label}_collision_model.out").read_text()
run("no_torque_singular", base.format(extra="").replace("species.in", "zero_rotation.in")
    .replace("comptorque=T", "comptorque=F"), expected_error="status=2")
(scratch / "incompatible.in").write_text((scratch / "zero_rotation.in").read_text()
    .replace("density_gradient_m3=0,0,0", "density_gradient_m3=-1e15,-1e15,-2e15"))
run("null_mode", base.format(extra="").replace("species.in", "incompatible.in"),
    expected_error="status=2")

# Omitted/off controls preserve a legacy trapped-only model that Lorentz rejects.
legacy = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,supban=T,nopassing=T,noshear=T,inp_swi=8,
vsteps=1,vmax_over_vth=1,mth_max_abs=0,log_level=-1,output_format='both',
krook_collisions=T,krook_nu=1000,krook_eta_steps=3,krook_rtol=1e-4 {extra} /
"""
run("legacy_omitted", legacy.format(extra=""))
run("legacy_off", legacy.format(extra=",lorentz_pitch_collisions=F,lorentz_nu=87"))
suffixes = ["_species.out", "_mixture_torque.out"]
for name in ("D", "T", "e"):
    suffixes += [f"_{name}{suffix}" for suffix in
                 ("_magfie_param.out", "_magfie.out", ".out", "_torque.out",
                  "_integral.out", "_torque_integral.out")]
for suffix in suffixes:
    left, right = scratch / f"legacy_omitted{suffix}", scratch / f"legacy_off{suffix}"
    assert left.read_bytes() == right.read_bytes(), suffix
parent_path = sys.argv[4] if len(sys.argv) == 5 else os.environ.get("NEORT_PARENT_EXECUTABLE")
if parent_path:
    parent = Path(parent_path)
    run("legacy_parent", legacy.format(extra=""), binary=parent)
    for suffix in suffixes:
        left, right = scratch / f"legacy_omitted{suffix}", scratch / f"legacy_parent{suffix}"
        assert left.read_bytes() == right.read_bytes(), f"parent {suffix}"
    with h5py.File(scratch / "legacy_parent.h5", "r") as parent_output, \
            h5py.File(scratch / "legacy_omitted.h5", "r") as current_output:
        def compare_parent(name, obj):
            if name == "created":
                return
            assert name in current_output, f"parent dataset/group missing: {name}"
            current = current_output[name]
            if isinstance(obj, h5py.Dataset):
                np.testing.assert_array_equal(current[()], obj[()], err_msg=name)
            for attribute, value in obj.attrs.items():
                np.testing.assert_array_equal(current.attrs[attribute], value,
                                              err_msg=f"{name}/{attribute}")
        parent_output.visititems(compare_parent)
print("Native Lorentz D/T/e CLI, metadata, guards and disabled legacy parity pass")
