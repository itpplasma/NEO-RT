"""Exercise native test-particle FP inputs, physical work and serialized state."""
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
mass_kg=3.345228e-27,5.017842e-27,9.1091e-31,
density_m3=5e15,5e15,1e16,temperature_ev=2000,2000,2000,
density_gradient_m3=0,0,0,temperature_gradient_ev=0,0,0 /
"""
(scratch / "species.in").write_text(species_input)
base = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,noshear=T,inp_swi=8,vsteps=1,vmax_over_vth=2,
mth_max_abs=0,log_level=-1,output_format='both',tp_fp_collisions=T,
fp_rate_scale=.7,fp_speed_steps=6,fp_pitch_steps=6,fp_umax=4 {extra} /
"""
environment = os.environ | {"OMP_NUM_THREADS": "3", "OMP_THREAD_LIMIT": "4",
                            "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1"}


def run(name, deck, expected_error=None, binary=executable):
    (scratch / f"{name}.in").write_text(deck)
    result = subprocess.run([str(binary.resolve()), name], cwd=scratch,
                            env=environment, capture_output=True, text=True, timeout=45)
    output = result.stdout + result.stderr
    if expected_error is None:
        assert result.returncode == 0, f"{name}: {output[-4000:]}"
    else:
        assert result.returncode != 0, f"{name}: invalid input was accepted"
        assert expected_error in output, f"{name}: unexpected rejection {output[-2000:]}"


def scalar(group, path):
    return group[path][()].item()


def assert_zero_d(group):
    for location in ("transport", "harmonics"):
        for name in ("D11co", "D11ctr", "D11t", "D11", "D12co", "D12ctr", "D12t", "D12"):
            np.testing.assert_array_equal(group[f"{location}/{name}"][()], 0)


run("fp", base.format(extra=""))
population_total = np.zeros(3)
elementary_charge = 1.602176634e-19
with h5py.File(scratch / "fp.h5", "r") as output:
    assert set(output["species"]) == {"D", "T", "e"}
    for name in ("D", "T", "e"):
        species = output[f"species/{name}"]
        config, magfie = species["config"], species["magfie"]
        assert scalar(config, "tp_fp_collisions")
        assert not scalar(config, "lorentz_pitch_collisions")
        assert not scalar(config, "krook_collisions")
        for key, value in (("fp_rate_scale", .7), ("fp_speed_steps", 6),
                           ("fp_pitch_steps", 6), ("fp_umax", 4)):
            np.testing.assert_equal(scalar(config, key), value)
        energy = 2000 * elementary_charge
        thermal_speed = np.sqrt(2 * energy / scalar(species, "mass_kg"))
        np.testing.assert_allclose(scalar(magfie, "fp_thermal_energy_J"), energy, rtol=2e-14)
        np.testing.assert_allclose(scalar(magfie, "fp_vth_m_s"), thermal_speed, rtol=2e-14)
        np.testing.assert_allclose(scalar(magfie, "vth"), thermal_speed * 100, rtol=2e-14)
        np.testing.assert_allclose(scalar(magfie, "T"), 2000, rtol=2e-14)
        np.testing.assert_allclose(scalar(magfie, "fp_charge_over_c"),
                                   scalar(species, "charge_e") * elementary_charge / 10,
                                   rtol=2e-14)
        assert not scalar(species, "transport/susceptibility_valid")
        assert scalar(species, "transport/fp_auxiliary_valid")
        assert not scalar(species, "transport/fp_physical_response")
        assert_zero_d(species)
        work = (species["transport/fp_source_work_real"][()] +
                1j * species["transport/fp_source_work_imag"][()])
        assert work.shape == (3, 2)  # HDF5 reverses the Fortran (2,3) dimensions.
        assert np.all(np.isfinite(work)) and np.all(np.abs(work) > 0)
        # With zero profile gradients, beta*A1=-n*OmegaE/Terg. This independently
        # checks the physical torque normalization, signed charge and branch order.
        density_cgs = scalar(species, "density_m3") / 1e6
        flux_scale = abs(scalar(magfie, "psi_pr") / scalar(magfie, "q"))
        prefactor = -3 * np.sqrt(np.pi) * density_cgs * thermal_speed * 100 * flux_scale
        physical_work = 1j * 3 * scalar(magfie, "Om_tE") / (energy * 1e7) * work[:, 0]
        expected = prefactor * physical_work.imag
        expected = expected[[1, 2, 0]]  # Native work order: trapped, co, counter.
        np.testing.assert_array_equal(species["harmonics/mth"][()].ravel(), [0, -6])
        torque = []
        for index, component in enumerate(("Tco", "Tctr", "Tt")):
            values = species[f"harmonics/{component}"][()].ravel()
            total = scalar(species, f"torque/{component}")
            assert np.isfinite(total) and total != 0
            np.testing.assert_allclose(total, expected[index], rtol=3e-13)
            np.testing.assert_allclose(values.sum(), total, rtol=2e-14)
            assert values[1 if component == "Tt" else 0] == 0
            torque.append(total)
        population_total += torque
    assert scalar(output, "species/e/charge_e") == -1
    np.testing.assert_allclose(output["mixture/torque_by_population"][()].ravel(),
                               population_total, rtol=2e-14)
    np.testing.assert_allclose(scalar(output, "mixture/torque_total"),
                               population_total.sum(), rtol=2e-14)
text = np.loadtxt(scratch / "fp_mixture_torque.out")
np.testing.assert_allclose(text[1:4], population_total, rtol=2e-14)
for name in ("D", "T", "e"):
    metadata = (scratch / f"fp_{name}_collision_model.out").read_text()
    assert "susceptibility_valid=F" in metadata
    assert "placeholder" in metadata.lower()

# Equal reference temperatures can become unequal at the evaluated surface.
for name, contents in (
    ("unequal_temperature", species_input.replace("2000,2000,2000", "2000,2100,2000")),
    ("unequal_local_temperature", species_input.replace("reference_s=.5", "reference_s=.25")
     .replace("temperature_gradient_ev=0,0,0", "temperature_gradient_ev=0,100,0")),
):
    (scratch / f"{name}.species").write_text(contents)
    run(name, base.format(extra="").replace("species.in", f"{name}.species"),
        "requires exactly equal temperatures")
for name, replacement, message in (
    ("zero_rate", ("fp_rate_scale=.7", "fp_rate_scale=0"), "fp_rate_scale must be finite"),
    ("negative_rate", ("fp_rate_scale=.7", "fp_rate_scale=-1"), "fp_rate_scale must be finite"),
    ("nan_rate", ("fp_rate_scale=.7", "fp_rate_scale=NaN"), "fp_rate_scale must be finite"),
    ("bad_umax", ("fp_umax=4", "fp_umax=0"), "fp_umax must be finite"),
    ("bad_speed", ("fp_speed_steps=6", "fp_speed_steps=1"), "steps must be at least 2"),
    ("bad_pitch", ("fp_pitch_steps=6", "fp_pitch_steps=1"), "steps must be at least 2"),
    ("missing_mixture", ("species_file='species.in'", "species_file=''"),
     "test-particle FP requires an explicit species mixture"),
):
    run(name, base.format(extra="").replace(*replacement), message)
for name, extra, message in (
    ("simultaneous_krook", ",krook_collisions=T,krook_nu=1", "select only one of TP FP"),
    ("simultaneous_lorentz", ",lorentz_pitch_collisions=T", "select only one of TP FP"),
    ("drift_phase", ",drift_phase=T", "linear common-phase thin orbits and all branches"),
    ("fow", ",fow_radial=T", "linear common-phase thin orbits and all branches"),
    ("trapped_only", ",nopassing=T", "linear common-phase thin orbits and all branches"),
    ("shaing", ",supban=T", "linear common-phase thin orbits and all branches"),
    ("nonlinear", ",nonlin=T", "general mixtures require linear transport"),
):
    run(name, base.format(extra=extra), message)

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

# R=0 needs the physical combined source and one global density gauge. Exact
# zero forcing is compatible even though the two auxiliary sources are not.
(scratch / "zero_rotation.species").write_text(species_input.replace("omega_e=-700", "omega_e=0"))
run("zero_drive", base.format(extra="").replace("species.in", "zero_rotation.species"))
with h5py.File(scratch / "zero_drive.h5", "r") as output:
    for name in ("D", "T", "e"):
        species = output[f"species/{name}"]
        assert not scalar(species, "transport/susceptibility_valid")
        assert not scalar(species, "transport/fp_auxiliary_valid")
        assert scalar(species, "transport/fp_physical_response")
        assert_zero_d(species)
        for component in ("Tco", "Tctr", "Tt"):
            assert scalar(species, f"torque/{component}") == 0
        for component in ("real", "imag"):
            np.testing.assert_array_equal(species[f"transport/fp_source_work_{component}"][()], 0)

incompatible = species_input.replace("omega_e=-700", "omega_e=0").replace(
    "density_gradient_m3=0,0,0", "density_gradient_m3=-5e15,-5e15,-1e16")
(scratch / "incompatible_source.species").write_text(incompatible)
run("incompatible_source", base.format(extra="").replace(
    "species.in", "incompatible_source.species"), "status=2")

# Omitted and disabled controls preserve an existing finite-Krook model.
legacy = """&params s=.5,species_file='species.in',epsmn=.001,m0=-6,mph=3,
comptorque=T,magdrift=F,supban=T,nopassing=T,noshear=T,inp_swi=8,
vsteps=1,vmax_over_vth=1,mth_max_abs=0,log_level=-1,output_format='both',
krook_collisions=T,krook_nu=1000,krook_eta_steps=3,krook_rtol=1e-4 {extra} /
"""
run("legacy_omitted", legacy.format(extra=""))
run("legacy_off", legacy.format(extra=",tp_fp_collisions=F"))
suffixes = ["_species.out", "_mixture_torque.out"]
for name in ("D", "T", "e"):
    suffixes += [f"_{name}{suffix}" for suffix in
                 ("_magfie_param.out", "_magfie.out", ".out", "_torque.out",
                  "_integral.out", "_torque_integral.out")]


def compare_runs(reference, current):
    for suffix in suffixes:
        assert (scratch / f"{reference}{suffix}").read_bytes() == (
            scratch / f"{current}{suffix}").read_bytes(), suffix
    with h5py.File(scratch / f"{reference}.h5", "r") as old, \
            h5py.File(scratch / f"{current}.h5", "r") as new:
        def compare(name, obj):
            if name == "created":
                return
            assert name in new, f"missing old HDF5 entry: {name}"
            if isinstance(obj, h5py.Dataset):
                np.testing.assert_array_equal(new[name][()], obj[()], err_msg=name)
            for attribute, value in obj.attrs.items():
                np.testing.assert_array_equal(new[name].attrs[attribute], value,
                                              err_msg=f"{name}/{attribute}")
        old.visititems(compare)


compare_runs("legacy_omitted", "legacy_off")
parent_path = sys.argv[4] if len(sys.argv) == 5 else os.environ.get("NEORT_PARENT_EXECUTABLE")
if parent_path:
    parent = Path(parent_path)
    run("legacy_parent", legacy.format(extra=""), binary=parent)
    compare_runs("legacy_parent", "legacy_omitted")
    lorentz = base.format(extra="").replace("tp_fp_collisions=T", "lorentz_pitch_collisions=T")
    lorentz = lorentz.replace("fp_rate_scale=.7,fp_speed_steps=6,fp_pitch_steps=6,fp_umax=4",
                              "lorentz_rate_source=1,lorentz_rate_scale=.7,lorentz_pitch_steps=3")
    run("lorentz_current", lorentz)
    run("lorentz_parent", lorentz, binary=parent)
    compare_runs("lorentz_parent", "lorentz_current")
print("Native TP FP D/T/e CLI, exact SI state, physical work, guards and disabled parity pass")
