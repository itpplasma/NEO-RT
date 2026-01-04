# NEO-RT

NEO-RT computes neoclassical toroidal viscosity (NTV) torque in resonant transport regimes using a Hamiltonian approach. The code reads Boozer-formatted magnetic field data, solves guiding-centre trajectories including resonant interactions, and evaluates transport coefficients and torque densities for individual flux surfaces.

## Key capabilities

- Modular Fortran implementation with a `neo_rt.x` solver and `neo_rt_diag.x` diagnostics executable.
- Standalone magnetic-field interface (`do_magfie_standalone`) with optional coupling to NEO-2 through CMake options.
- Support for trapped and passing particle contributions, optional torque integration, and nonlinear corrections that depend on thermodynamic profiles.
- Batch-processing helpers and Python utilities for scanning multiple flux surfaces and plotting output torque profiles.

## Dependencies

The project is developed with GNU Fortran and CMake. The following dependencies are required to build the solver and diagnostics:

- GNU Fortran (gfortran) and a C/C++ toolchain (`build-essential`).
- CMake 3.22+ and Ninja.
- BLAS and LAPACK implementations.
- SuiteSparse (used by the magnetic-field routines).
- NetCDF (Fortran interface, provides `nf-config`).
- Python 3 with NumPy (for the regression tests and helper scripts).
- Optional: [FORD](https://github.com/Fortran-FOSS-Programmers/ford) and Graphviz for documentation generation.

On Debian/Ubuntu systems you can install the toolchain with:

```bash
make deps
```

This installs the runtime dependencies required for the default standalone build. NetCDF must provide the `nf-config` helper that CMake queries while configuring the project.

## Building

The default build uses CMake and Ninja via the provided `Makefile` wrapper. From the repository root run:

```bash
make
```

This creates `build/neo_rt.x` and `build/neo_rt_diag.x` together with intermediate static libraries. You can control the build configuration through `CONFIG=Release|Debug|Fast` when invoking `make`. To rebuild after changing configuration flags run `make clean` followed by `make CONFIG=Debug`.

CMake also exposes an option `USE_STANDALONE` (enabled by default) that uses the standalone magnetic-field reader. Set `-DUSE_STANDALONE=OFF` when configuring to link against an external NEO-2 checkout if required.

## Running simulations

The solver operates on a single flux surface at a time. After building, execute the main program as

```bash
./build/neo_rt.x <runname>
```

The solver expects a namelist file `<runname>.in` that specifies the simulation parameters (see [`doc/running.md`](doc/running.md) for the complete description). Magnetic-field data are read from the Boozer files `in_file` (axisymmetric background) and, if `pertfile=.true.`, `in_file_pert` (non-axisymmetric perturbation) located in the working directory. Example input decks are provided in `examples/base`.

### Batch scans

For scans over several flux surfaces, use the helper script:

```bash
python3 python/run_driftorbit.py --exe ./build/neo_rt.x 0 5
```

The script fills in template values from `driftorbit.in.template` using `profile.in` data and runs the solver for each requested flux surface. The executable path can also be supplied via the environment variable `NEORT_EXECUTABLE`.

## Output overview

Each run produces a set of plain-text files summarised below:

| File | Description |
| --- | --- |
| `<runname>.out` | Transport coefficients (`D11`, `D12`) separated into trapped, co-passing, and counter-passing contributions. |
| `<runname>_integral.out` | Velocity-space integrals and resonance bounds for each evaluated harmonic. |
| `<runname>_magfie_param.out` | Magnetic-field and geometric parameters sampled on the target surface. |
| `<runname>_magfie.out` | Magnetic-field line samples including perturbation amplitudes. |
| `<runname>_torque.out` | Torque density components (written when `comptorque=.true.`). |
| `<runname>_torque_integral.out` | Harmonic-resolved torque contributions (written when `comptorque=.true.`). |

Detailed descriptions of all inputs and outputs are collected in [`doc/running.md`](doc/running.md) and [`doc/file_formats.md`](doc/file_formats.md).

## Diagnostics utilities

The build also produces `neo_rt_diag.x`, which links against the main library and the [fortplot](https://github.com/lazy-fortran/fortplot) plotting backend. Diagnostic drivers in `src/diag` provide additional post-processing routines for debugging, harmonic inspection, and attenuation maps.

## Library interface

NEO-RT can be embedded into other codes via the `neort_lib` module. The API supports thread-safe parallel computation across flux surfaces:

```fortran
use neort_lib

call neort_init("config", "in_file", "in_file_pert")
call neort_prepare_splines_from_files("plasma.in", "profile.in")

!$omp parallel do
do i = 1, n
    call neort_compute_at_s(s_array(i), results(i))
end do
```

For time-evolving profiles, call `neort_prepare_splines` with updated arrays before each batch of computations. See [`doc/library.md`](doc/library.md) for the complete API reference and usage examples.

## Testing

Regression and unit tests are provided for both the Fortran and Python components. After building the solver, run

```bash
make test
```

which executes the CTest suite together with the Python regression test in `test/ripple_plateau`.

## Documentation

FORD-based API documentation can be generated locally with

```bash
make doc
```

This command runs `ford doc.md` and writes the HTML output to `build/doc`. The same command is executed on GitHub Actions to publish documentation to GitHub Pages. Ensure that the `ford` Python package and Graphviz binaries are available before running the target locally.

## License

NEO-RT is distributed under the terms of the [MIT License](LICENSE).
