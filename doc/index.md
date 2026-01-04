---
title: NEO-RT documentation
---

# NEO-RT documentation

Welcome to the NEO-RT documentation hub. This site accompanies the FORD-generated API reference and collects the practical information needed to configure and run the solver.

## Quick start

1. Install the build dependencies (`gfortran`, CMake, Ninja, BLAS/LAPACK, SuiteSparse, NetCDF) and optional documentation tools (`ford`, Graphviz).
2. Configure and build the project:
   ```bash
   make
   ```
3. Prepare a namelist `<runname>.in`, Boozer files `in_file` and `in_file_pert` (optional), and any required profile inputs.
4. Run a single flux-surface calculation:
   ```bash
   ./build/neo_rt.x <runname>
   ```
5. Inspect the generated outputs (`<runname>.out`, `<runname>_integral.out`, etc.) or invoke the diagnostics executable `neo_rt_diag.x` for additional analyses.

## Input reference

The solver consumes the files summarised below. See [Running NEO-RT](running.md) for detailed descriptions and formatting requirements.

- `<runname>.in` – Fortran namelist defining magnetic and kinetic parameters.
- `in_file` / `in_file_pert` – Boozer-formatted magnetic-field data.
- `plasma.in` – Thermodynamic profiles for density and temperature (mandatory for torque and nonlinear runs).
- `profile.in` – Rotation profile used for scans and nonlinear corrections.
- `driftorbit.in.template` – Template filled by the batch helper script.

The helper script `python/run_driftorbit.py` accepts `--exe`, `--template`, `--profile`, and `--prefix` options and defaults to the executable named in `NEORT_EXECUTABLE` or `neo_rt.x`.

## Output overview

The solver writes a consistent set of text files for each run:

- `<runname>.out` – Transport coefficients with trapped/co-/counter-passing breakdown.
- `<runname>_integral.out` – Harmonic-resolved velocity integrals and resonance bounds.
- `<runname>_magfie_param.out` – Flux-surface magnetic parameters for verification.
- `<runname>_magfie.out` – Magnetic-field samples along the Boozer poloidal angle.
- `<runname>_torque.out` / `<runname>_torque_integral.out` – Torque diagnostics when enabled.

Column definitions and units are tabulated in [Output file formats](file_formats.md).

## Further reading

- [Running NEO-RT](running.md) – Complete workflow for preparing inputs, executing runs, and collecting results.
- [Library Interface](library.md) – Using NEO-RT as a library for embedding transport calculations in other codes.
- `src/diag/*.f90` – Diagnostics drivers linked by `neo_rt_diag.x` for harmonic inspection and attenuation studies.
- `examples/` – Sample input decks and Python utilities for plotting transport coefficients and torque profiles.

The FORD-generated API pages include module, type, and procedure documentation extracted from the source code for deeper technical reference.
