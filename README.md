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

On GNU/Linux with GNU Fortran, the build now prefers OpenBLAS when it is
available so standalone NEO-RT builds match the shared benchmark venv instead
of picking up a stale system MKL installation.

CMake exposes `USE_STANDALONE` (enabled by default) for the supported standalone
magnetic-field reader. That default does not fetch NEO-2 or require an external
sparse factorizer. The optional `USE_STANDALONE=OFF` dependency is pinned to
NEO-2 `e9e2877d852a6a29b228cdf7bc8f805ec3a47ee3`, whose default sparse backend is
native CPU FortSparse. An explicit `NEO-2_REF` still overrides that pin.

The optional adapter delegates field interpolation and perturbation Fourier
coefficients to NEO-2's explicit provider lifecycle. It supports the existing
axisymmetric Boozer formats 8/9, independent perturbation-file format selection,
CGS field scaling and file-header minor radius. Provider calls are serialized
because NEO-2 retains one global field instance; different flux surfaces may be
evaluated by NEO-RT worker threads, but independent NEO-2 calculations cannot
share that state concurrently. This does not add distributed field ownership.
Chartmap input (`inp_swi=10`) remains a standalone capability.

The existing zero-radial-covariant-field approximation remains explicit, and the
toroidal component of curl is marked unavailable. The pinned provider fixes the
laboratory-10 Jacobian orientation and computes the complete CGS covariant curl
using derivatives of the same cached field-value polynomials. NEO-RT retains its
declared approximation and does not enable the raw provider's general 3-D field
capability. The provider's cubic angular interpolation uses 1025
poloidal and 9 toroidal points; the independent fixtures qualify the stated
checks, not a rigorous error enclosure for arbitrary input spectra.

Configure `-DUSE_STANDALONE=OFF -DNEORT_BUILD_TESTING=ON`, build `neo_rt.x`, and
run the focused CTest cases `test_boozer_handedness`, `test_hcurl_axisymmetric`,
`test_perturbation_fourier_kernel`, `test_neo_provider_parallel` and `test_bounce`.
The optional test configuration registers checks against the shared field
interface; standalone chartmap and storage-layout fixtures stay on the default
path.

To build against a specific libneo branch, tag, or commit, pass `-DLIBNEO_REF=<ref>` at configure time or `LIBNEO_REF=<ref>` on the make command line:

```bash
make LIBNEO_REF=my-branch
# or directly
cmake -S . -B build -DLIBNEO_REF=my-branch
```

Unset (the default) falls back to the automatic ref resolution in `cmake/Util.cmake`.

### Python virtual environment for helpers and tests

The solver itself is Fortran, but the helper scripts, examples, and regression
tests use Python packages such as NumPy, SciPy, Matplotlib, `f90nml`, and
`h5py`. On systems with an externally managed Python, create a repository-local
virtual environment:

```bash
./setup-venv.sh
```

Later, reactivate it with:

```bash
source .venv/bin/activate
```

With `.venv` active, `make test` uses the expected Python interpreter for the
pytest-based checks.

If you are running the benchmark checkout alongside the sibling repositories
`../benchmark-simple-potato`, `../SIMPLE`, and `../NEO-RT`, prefer the shared
environment in `../benchmark-simple-potato/.venv` so both codes use the same
Python installation.

## Running simulations

The solver operates on a single flux surface at a time. After building, execute the main program as

```bash
./build/neo_rt.x <runname>
```

The solver expects a namelist file `<runname>.in` that specifies the simulation parameters (see [`doc/running.md`](doc/running.md) for the complete description). Magnetic-field data are read from the Boozer files `in_file` (axisymmetric background) and, if `pertfile=.true.`, `in_file_pert` (non-axisymmetric perturbation) located in the working directory. Example input decks are provided in `examples/base`.

### Batch scans

For scans over several flux surfaces, use the helper script:

```bash
python python/run_driftorbit.py --exe ./build/neo_rt.x 0 5
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

NEO-RT can be embedded into other codes via the `neort_lib` module. The API exposes both the spline-based mode and a config-only mode while remaining thread-safe for parallel loops:

```fortran
use neort_lib

type(config_t) :: cfg

cfg%s = 0.5d0
cfg%comptorque = .true.

call neort_init(cfg, "in_file", "in_file_pert")
call neort_prepare_splines("plasma.in", "profile.in")

!$omp parallel do
do i = 1, n
    call neort_compute_at_s(s_array(i), results(i))
end do
```

For time-evolving profiles, call `neort_prepare_splines` with updated arrays before each batch of computations. See [`doc/library.md`](doc/library.md) for the complete API reference and usage examples.

For single-point runs such as ripple tests, populate a `config_t` instance, call `neort_init(cfg, "in_file")`, and use `neort_compute_no_splines` instead of preparing splines.

When spline mode is active the `M_t`, `vth`, `qi`, and `mi` values come from the supplied plasma/profile data and overwrite any prior config settings, so reset them if you later switch back to the config-only workflow.

## Examples

The `examples/` directory contains ready-to-run input decks. To run the base example:

```bash
make examples-base
```

To list all available examples:

```bash
make examples-list
```

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
