# NEO-RT

NEO-RT computes neoclassical toroidal viscosity (NTV) torque in resonant transport regimes using a Hamiltonian approach. The code reads Boozer-formatted magnetic field data, solves guiding-centre trajectories including resonant interactions, and evaluates transport coefficients and torque densities for individual flux surfaces.

## Key capabilities

- Modular Fortran implementation with a `neo_rt.x` solver and `neo_rt_diag.x` diagnostics executable.
- Standalone magnetic-field interface (`do_magfie_standalone`).
- Support for trapped and passing particle contributions, optional torque integration, and nonlinear corrections that depend on thermodynamic profiles.
- Batch-processing helpers and Python utilities for scanning multiple flux surfaces and plotting output torque profiles.

`neo_rt_field_probe.x BCF INPUT_SWITCH POINTS OUTPUT` reads the native Boozer
field at whitespace triples `(s_tor, phi_B, theta_B)` without initializing
profiles or transport. Input switches 8 and 9 are supported; angles are radians.
Its CSV reports field magnitude, covariant components, iota, Jacobian and flux
derivative in the native cgs units stated in the header. Out-of-support radial
queries fail. `tools/run_field_probe.py` also accepts a CSV point list and writes
input, executable and output hashes beside the result.

`neo_rt_profile_probe.x BCF INPUT_SWITCH PLASMA MACH EFAC BFAC POINTS OUTPUT`
executes the native plasma, Mach and collision initializers without orbits or
transport. `PLASMA` has the solver's species header and six columns; `MACH` is
its numeric two-column table. `POINTS` contains one `s_tor` per line. Queries
must lie in the common carrier/plasma/Mach support. The CSV reports native
profile values and derivatives, Mach/thermal speed/rotation, actual carrier
R0, orbit masses/charge and collision arrays. Electron density, Zeff and
thermal pressure are diagnostic combinations of the loaded ion background.
Temperatures remain eV, densities cm^-3, velocities cm/s, masses g, pressure
erg/cm^3, and radial derivatives are with respect to `s_tor`.

The header exposes the distinct orbit and legacy collision mass/energy
constants; no constants are unified. Collision arrays retain their native
normalization (`efcolf` inverse free path, `velrat` velocity ratios, `enrat`
energy ratios). Positive background ion species are required; a zero second
ion density is supported. This readback sets `nonlin=false` and grants no
collision-model or torque admission. `tools/run_profile_probe.py` accepts the
same CSV points as the field wrapper, extracts `s_tor` in the supplied order,
and records exact carrier, plasma, Mach, query, executable and output hashes.
Use distinct outputs to retain historical and revised backgrounds.

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

On GNU/Linux with GNU Fortran, the build now prefers OpenBLAS when it is
available so standalone NEO-RT builds match the shared benchmark venv instead
of picking up a stale system MKL installation.

CMake uses the standalone magnetic-field reader (`USE_STANDALONE=ON`, the
default). `-DUSE_STANDALONE=OFF` is currently unsupported and fails immediately
at configure time. Restoring native coupling requires a compatible
[NEO-2](https://github.com/itpplasma/NEO-2) field/perturbation API, an updated
NEO-RT adapter, and CMake integration that provides
[libneo](https://github.com/itpplasma/libneo)'s `LIBNEO::hdf5_tools` target.
Fetching libneo alone does not restore the adapter.

Verified against NEO-2 main at 967ca2ef541955589439ad620f48bd0afe9210e1:
`src/do_magfie_neo.f90:85` imports `neo_magfie_pert_amp`, which NEO-2 does
not export, and its copied kernel accesses NEO-2's private spline arrays
(`NEO-2-QL/neo_magfie_perturbation.f90:35-67`). NEO-RT currently resolves
NEO-2 by branch (normally `main`), rather than pinning a NEO-2 commit.
The native integration changes in
[NEO-2 #212](https://github.com/itpplasma/NEO-2/pull/212),
[#213](https://github.com/itpplasma/NEO-2/pull/213),
[#232](https://github.com/itpplasma/NEO-2/pull/232) and
[NEO-RT #198](https://github.com/itpplasma/NEO-RT/pull/198) remain deferred;
this build guard does not include or depend on them.

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

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

## Documentation

FORD-based API documentation can be generated locally with

```bash
make doc
```

Install the pinned documentation tooling with `pip install -r requirements-docs.txt`
and ensure that Graphviz is available. This command runs `ford doc.md` and writes
the HTML output to `build/doc`, without configuring or compiling the solver.
FORD documents `src` recursively using the default standalone field backend;
the alternative NEO-2 implementation is excluded because it defines the same
module names. GitHub Actions builds the documentation on pull requests and
merge groups, and publishes it to GitHub Pages only on pushes to `main`.

## License

NEO-RT is distributed under the terms of the [MIT License](LICENSE).
