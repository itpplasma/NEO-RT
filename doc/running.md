# Running NEO-RT

## Required Files

### Parameter File
Namelist file according to `examples/neort.in`.

### Boozer Files for Background Field and Perturbation
Boozer type files, i.e., text files with magnetic field data in a specific layout. Names are hardcoded as `in_file` and `in_file_pert`.

### Rotation profile (`profile.in`)
Text file with data in three columns:
1. Radial coordinate (Boozer `s`).
2. Toroidal Mach number.
3. Thermal velocity in ?.

There is an Octave function `neo2_to_NEO-RT` as part of NEO-2, to create a NEO-RT profile file from a NEO-2 flux surface scan output.

### Thermodynamic profiles (`plasma.in`)
Text file with a special format (and fixed name `plasma.in`?):
- First and third lines are comments (starting with `%`), indicating quantities in the following line(s).
- Second line contains five values:
  1. Number of flux surfaces.
  2. Mass of the first ion species (in atomic mass units).
  3. Mass of the second ion species.
  4. Charge of the first ion species (in elementary charge units).
  5. Charge of the second ion species.
- Starting from the fourth line, values for quantities that vary across flux surfaces are listed. The number of these lines should match the value given in the second line.
- Each of these lines has six columns:
  1. Flux surface label (Boozer `s`).
  2. Density of the first ion species (in `1/m^3`).
  3. Density of the second ion species.
  4. Temperature of the first ion species (in `eV`).
  5. Temperature of the second ion species.
  6. Temperature of electrons.

## Running
- For a single run, simply execute NEO-RT in the folder with the input files.
- For scans, you can use the (executable) Python3 script `run_driftorbit.py`, which requires the first (inclusive) and last (exclusive) flux surface that should be used.

## Output
**Note:** For simplicity, no flux surface numbering is done here. The number would appear after the runname, e.g., `runname42.out` or `runname42_torque.out`.

### `runname.out`
Toroidal Mach number and transport coefficients. There are nine columns:
1. Toroidal Mach number.
2. Part of D11 due to co-passing particles.
3. Part of D11 due to counter-passing particles.
4. Part of D11 due to trapped particles.
5. Sum of the three parts of D11.
6-9. Analogous for D12.

### `runname_torque.out`
Boozer coordinate and...

## Collecting Results
Can be done with the (executable) Python3 script `collect_data_from_individual_runs.py`.

# Running Potato

## Input
- Requires an EFIT file, which can be created from a Boozer file via the program `boozer_to_efit` (also part of NEO-RT).
- Profile input is needed, which can be created from NEO-2 output via the script `convert_profiles.m`. It creates coefficients of a polynomial fit to the profiles for various degrees. One should be selected and provided as input with the name `profile_polyIn`. We suggest linking the selected file.

## Output
File `eqmagprofs.dat`:
1. Poloidal radius.
2. Toroidal radius.
3. Poloidal flux.
4. Toroidal momentum.
5. Safety factor.
6. Averaged nabla psi.
7. (Flux?) surface area.
