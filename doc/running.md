---
title: Running NEO-RT
---

# Running NEO-RT

This guide summarises the required inputs, run workflow, and produced outputs for `neo_rt.x`. The solver evaluates neoclassical transport for a single flux surface per execution. Batch scans can be scripted to cover ranges of flux surfaces.

## Build outputs

After `make` the build directory contains two executables:

- `neo_rt.x` – main solver.
- `neo_rt_diag.x` – diagnostics driver that links against the solver library and the `fortplot` plotting backend.

Ensure that the build directory is on your `PATH` or provide the full path to these executables when running from another folder.

## Preparing input data

A simulation expects the following files in the working directory:

### Parameter file `<runname>.in`

The parameter file is a Fortran namelist `&params` with the fields listed below. An example is provided in `examples/driftorbit.in`.

| Name | Meaning | Notes |
| --- | --- | --- |
| `s` | Normalised toroidal flux of the target surface. | Boozer `s` coordinate. |
| `M_t` | Toroidal Mach number. | Scaled internally by `efac/bfac`. |
| `qs` | Ion charge in units of the elementary charge. | Used to compute `qi`. |
| `ms` | Ion mass in atomic mass units. | Used to compute `mi`. |
| `vth` | Thermal velocity `[cm/s]`. | Used when no `plasma.in` is provided. |
| `epsmn` | Perturbation amplitude `B₁/B₀`. | Used when `pertfile=.false.`. |
| `m0` | Poloidal perturbation harmonic. | Integer. |
| `mph` | Toroidal perturbation harmonic. | Positive integer when using analytic perturbations. |
| `magdrift` | Include magnetic drift (`.true.`/`.false.`). | Controls bounce-averaged drifts. |
| `nopassing` | Skip passing resonances. | When `.true.` only trapped contributions are computed. |
| `noshear` | Neglect magnetic shear. | Propagated to orbit module. |
| `pertfile` | Read perturbation from `in_file_pert`. | Otherwise analytic perturbation is used. |
| `nonlin` | Enable nonlinear corrections. | Requires `plasma.in` and `profile.in`. |
| `comptorque` | Write torque diagnostics. | Produces `_torque.out` and `_torque_integral.out`. |
| `bfac` | Magnetic-field scaling factor. | Multiplies `B`. |
| `efac` | Electric-field scaling factor. | Multiplies `E`. |
| `inp_swi` | Boozer input switch. | Passed to `do_magfie`. |
| `vsteps` | Number of velocity grid points. | Set to `0` for adaptive quadrature. |
| `log_level` | Verbosity level for the logger. | Defined in `src/logging.f90`. |

### Magnetic field data

- `in_file` – Axisymmetric Boozer magnetic-field data used by `do_magfie_standalone`.
- `in_file_pert` – Optional perturbation file; required when `pertfile=.true.`.
- `thetafun_inp.dat` – Optional data for attenuation-factor diagnostics (`attenuation_factor.f90`).

Example Boozer files are provided in `examples/base`.

### Thermodynamic profiles (`plasma.in`)

Required when `nonlin=.true.` or `comptorque=.true.`. The file format matches the expectations of `neort_profiles:init_plasma_input`:

1. Comment line (ignored).
2. `nflux  am1  am2  Z1  Z2` – number of flux surfaces, ion masses (amu), and charges.
3. Comment line (ignored).
4. `nflux` data rows with six columns:
   - $s$ – Normalised toroidal flux.
   - $n_1$ – Density of species 1 (1/cm³).
   - $n_2$ – Density of species 2 (1/cm³).
   - $T_1$ – Temperature of species 1 (eV).
   - $T_2$ – Temperature of species 2 (eV).
   - $T_e$ – Electron temperature (eV).

The solver creates spline coefficients from these values to evaluate densities, temperatures, and their radial derivatives at the requested flux surface.

### Rotation profile (`profile.in`)

Required for nonlinear runs and useful when scanning multiple flux surfaces. Each line contains three columns:

1. `s` – Normalised toroidal flux.
2. `M_t` – Toroidal Mach number profile.
3. `vth` – Thermal velocity profile (retained for compatibility; the solver evaluates `vth` from `plasma.in` when available).

`python/run_driftorbit.py` reads this file to populate template inputs for scans.

## Running single-surface calculations

From the directory containing the input files run:

```bash
./build/neo_rt.x <runname>
```

The program reads `<runname>.in`, `in_file`, and optional files as described above, then writes all outputs with the same `<runname>` prefix into the working directory. Progress information and warnings are printed to standard output.

## Batch processing helper

The helper script automates scans over consecutive flux surfaces:

```bash
python3 python/run_driftorbit.py --exe ./build/neo_rt.x <lower_index> [<upper_index>]
```

- `--exe` selects the solver executable. If omitted, the script uses the value of the `NEORT_EXECUTABLE` environment variable or defaults to `neo_rt.x`.
- `<lower_index>` selects the first flux surface index (zero-based).
- `<upper_index>` is optional; when provided the script processes indices `[lower, upper)`. Without it only `lower_index` is run.

The script expects `driftorbit.in.template` and `profile.in` in the working directory. Log and error output for each run are saved as `<runname>.log` and `<runname>.err`.

## Generated outputs

The solver writes the files described below:

- `<runname>.out` – Transport coefficients split into trapped and passing contributions.
- `<runname>_integral.out` – Velocity-space integrals and resonance bounds for each poloidal harmonic.
- `<runname>_magfie_param.out` – Magnetic-field parameters evaluated on the flux surface.
- `<runname>_magfie.out` – Field-line samples, basis vectors, and perturbation amplitudes along the poloidal angle.
- `<runname>_torque.out` – Torque density components (only when `comptorque=.true.`).
- `<runname>_torque_integral.out` – Harmonic-resolved torque contributions (only when `comptorque=.true.`).

Detailed column descriptions are provided in [`doc/file_formats.md`](file_formats.md). Utility scripts in `examples`, such as `plot_torque.py`, offer quick post-processing templates.

## Diagnostics program

`neo_rt_diag.x` links the solver library with visualisation helpers. The drivers in `src/diag` expose additional CLI programs for debugging orbit trajectories, inspecting harmonic contributions, and plotting attenuation maps. Refer to the comments in `src/diag/*.f90` for usage examples.
