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

## Python helper environment

If you plan to use `python/run_driftorbit.py`, plotting scripts from
`examples/`, or `make test`, create the repository-local virtual environment
first:

```bash
./setup-venv.sh
source .venv/bin/activate
```

This keeps the Python dependencies out of the system interpreter and avoids PEP
668 issues on managed machines. Later, re-enter the same environment with the
same `source .venv/bin/activate` command.

For the sibling benchmark layout (`../benchmark-simple-potato`, `../SIMPLE`,
`../NEO-RT`), use `source ../benchmark-simple-potato/.venv/bin/activate` so the
benchmark, SIMPLE, and NEO-RT share one environment.

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
| `epsmn` | Perturbation amplitude `B₁/B₀`. | Used when `pertfile=.false.`; default `1`. |
| `pertfile_scale` | Scale factor for the perturbation read from `in_file_pert`. | Used when `pertfile=.true.`; default `1`. |
| `m0` | Poloidal perturbation harmonic. | Integer. |
| `mph` | Toroidal perturbation harmonic. | Positive integer when using analytic perturbations. |
| `magdrift` | Include magnetic drift (`.true.`/`.false.`). | Controls bounce-averaged drifts. |
| `magdrift_passing` | Include magnetic drift for passing particles (`1`/`0`; `-1` follows `magdrift`). | Independent passing-particle control. With `0`, the passing transit frequency and its `d_Om_ds` derivative omit the magnetic-drift term. |
| `drift_phase` | Keep the oscillating toroidal drift phase `mph*int(Om_tB(theta) - <Om_tB>)dt` in the orbit integral of the perturbation (`.false.` default). | The default is the thin-orbit ordering `Om_d << Om_th`, where this phase is `O(Om_d/Om_th)`. A passing resonance with `mth != 0` or a trapped bounce harmonic `mth != 0` needs `mph*<Om_d> ~ Om_th`, so the modulation index is `O(1)`: enable it there. The sidebands are then Bessel functions of the index. |
| `nopassing` | Skip passing resonances. | When `.true.` only trapped contributions are computed. |
| `noshear` | Neglect magnetic shear. | Propagated to orbit module. |
| `pertfile` | Read perturbation from `in_file_pert`. | Otherwise analytic perturbation is used. |
| `fow_radial` | Sample the perturbation amplitude (`B_n` and `B`) at the first-order guiding-centre radius `s + Delta s` (`.false.` default). | `Delta s = m c q v_par h_phi / (e sign_theta psi_pr)` from conservation of `p_phi`. `s` is the `p_phi` label of the orbit, `(e/c) A_phi(s) = p_phi`, which is also the label of NEO-RT's frequencies (the shear term of `Om_tB` is `q' <Delta s dtheta/dt>`), so the orbit mean of `Delta s` is kept: it vanishes for trapped orbits and is `O(+-rho_pol)` for co/counter-passing orbits. Frequencies and the resonance stay those of the label `s`. For passing contributions combine it with `drift_phase` (and `noshear = .false.`): sampling at `s + Delta s` also shifts the toroidal phase by `mph q' Delta s theta`, whose oscillating part `mph q' int(Delta s dtheta/dt - <.>)dt` is of the same first order as the amplitude correction and enters only through the shear term of the `drift_phase` integrand. Still omitted at first order: the poloidal timing shift (the change of `dtheta/dt` at `s + Delta s` and the poloidal drift), so the orbit phase uses the thin-orbit `theta(t)`. Sign, size and labelling are checked against the full guiding-centre orbit in `test_fow_gc_orbit`. Needs a radially resolved perturbation (`pertfile`) to matter. |
| `nonlin` | Enable nonlinear corrections. | Requires `plasma.in` and `profile.in`. |
| `profile_species` | Select profiles from `plasma.in`. | `1` first ion (default), `2` second ion, `0` electrons. |
| `species_file` | General neutral mixture with shared physical rotation. | Empty by default; an explicit SI species file computes every member in one invocation. See [shared species mixtures](species.md). |
| `krook_collisions` | Enable linear finite-width collision broadening. | Default `.false.`; requires `nonlin=.false.` for positive width. |
| `krook_nu` | Supplied Krook damping rate `[s^-1]`. | Default zero; no `2*pi` conversion. |
| `krook_model` | Rate model. | `0` constant, `1` energy and harmonic, `2` energy only; see [Krook collisions](krook.md) for harmonics, geometry factors and quadrature controls. |
| `comptorque` | Write torque diagnostics. | Produces `_torque.out` and `_torque_integral.out`. |
| `bfac` | Magnetic-field scaling factor. | Multiplies `B`. |
| `efac` | Electric-field scaling factor. | Multiplies `E`. |
| `inp_swi` | Boozer input switch. | Passed to `do_magfie`. |
| `vsteps` | Number of velocity grid points. | Set to `0` for adaptive quadrature. |
| `mth_max_abs` | Maximum absolute orbit-resonance harmonic. | Default `-1` keeps the historical `±ceil(2\|mph q\|)` range; `0` selects only `mth=0`. |
| `vmax_over_vth` | Upper velocity-space cutoff in units of `vth`. | Default `4.0` captures the far-tail resonance (the old `3.0` bound truncated it); set `3.0` to reproduce pre-2026-07-20 results; must be positive. |
| `monoenergetic_x` | Select one normalized kinetic energy, `x=E/T`. | Default `-1` performs the Maxwellian velocity integral; a positive value returns the delta-normalized monoenergetic kernel at `u=v/vth=sqrt(x)`. |
| `inp_swi_pert` | Input format for `in_file_pert` with the standalone field reader. | Default `-1` inherits `inp_swi`; set `9` to combine an axisymmetric chartmap (`inp_swi=10`) with a Strumberger perturbation `.bc`. |
| `log_level` | Verbosity level for the logger. | Defined in `src/logging.f90`. |

`magdrift_passing` is deliberately separate from `magdrift`: setting
`magdrift=.true., magdrift_passing=0` retains trapped magnetic drift while
matching a model whose passing operator has no magnetic-drift contribution.
The frequency and derivative paths use the same switch. In particular,
`d_Om_ds` does not add `bounceavg(3)*v**2` on the passing branch when the
switch is zero; this consistency guard was fixed in commit `0b7c4e8` and is
covered by `test/test_passing_drift_derivative.f90`.

### Magnetic field data

- `in_file` – Axisymmetric Boozer magnetic-field data used by `do_magfie_standalone`.
- `in_file_pert` – Optional perturbation file; required when `pertfile=.true.`.
- `thetafun_inp.dat` – Optional data for attenuation-factor diagnostics (`attenuation_factor.f90`).

Example Boozer files are provided in `examples/base`.

### Thermodynamic profiles (`plasma.in`)

For isolated-species runs, required when `nonlin=.true.` or `comptorque=.true.`.
General mixtures use [the explicit SI species file](species.md). The legacy file
format matches the expectations of `neort_profiles:read_plasma_input`:

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

Set `NEORT_RESONANCE_LEDGER=/path/to/roots.txt` with `comptorque=.true.` to
append one row per collisionless resonance root. The header names every column:
`species` is the one-based mixture index (1 for a single species), `branch` is
`cop`, `ctr`, or `trap`, and `dres_deta` is the derivative of
`mth*Omth + mph*Omph` at the root. `sample_weight` is the actual energy weight
(midpoint Maxwellian weight or 1 for a monoenergetic sample). `dTds` includes
that weight and nonlinear attenuation: sum it by surface, species, harmonic,
and branch to reconstruct native `_torque_integral.out` values. No extra
velocity step or attenuation factor is needed. Frequencies, bounce time,
perturbation weight, attenuation, and integrator status are also recorded.

OpenMP workers buffer rows and serialize each complete append transaction.
Use a separate pathname for each process; cross-process appends are unsupported.
Finite-width Krook transport has no discrete root ledger and rejects an enabled
ledger request. With the variable unset, no ledger is written.

`neo_rt_diag.x bounce_debug <runname> [mth] [all|cop|ctr|trap]` filters the bounce
probe by harmonic and branch and reports both frequencies, their resonance
residual, and its eta derivative. The speed scan now uses the configured
`vmax_over_vth` (default 4), replacing the former fixed upper limit of 3.
