# NEO-RT

NEO-RT computes neoclassical toroidal viscosity (NTV) torque and transport
coefficients in resonant transport regimes of tokamaks with small
non-axisymmetric magnetic perturbations. It uses a Hamiltonian formulation in
action-angle variables and treats superbanana-plateau and drift-orbit resonances
in a unified way. It is developed at the Institute of Theoretical and
Computational Physics, TU Graz.

Status: research code under active development; interfaces may change.

## Features

- Single flux-surface solver `neo_rt.x` reading Boozer-coordinate field data
  (axisymmetric background plus perturbation).
- Trapped, co-passing and counter-passing contributions to the transport
  coefficients `D11`, `D12`, and torque density when `comptorque=.true.`.
- Optional nonlinear attenuation depending on thermodynamic and rotation
  profiles.
- Fortran library interface (`neort_lib`) for embedding in other codes, usable
  in OpenMP-parallel loops over flux surfaces.
- Python helpers for flux-surface scans and plotting, and a diagnostics
  executable `neo_rt_diag.x`.

## Building

Requirements: gfortran with a C/C++ toolchain, CMake 3.22 or newer, Ninja,
BLAS/LAPACK, and Python 3 with NumPy for tests and helper scripts. On
Debian/Ubuntu, `make deps` installs the system packages. CMake fetches libneo,
fortnum, fortplot and the spline library at configure time.

```bash
make                 # Release build in build/
make CONFIG=Debug    # after make clean
```

This produces `build/neo_rt.x` and `build/neo_rt_diag.x`. Build options, the
Python environment and the status of the optional NEO-2 field backend are
described in [doc/building.md](doc/building.md).

## Usage

The solver works on one flux surface per run:

```bash
./build/neo_rt.x <runname>
```

It reads the namelist `<runname>.in` and the Boozer files `in_file` and, if
`pertfile=.true.`, `in_file_pert` from the working directory. Results are
written to `<runname>.out`, `<runname>_integral.out`, `<runname>_torque.out`
and related text files. A bundled example runs with

```bash
make examples-base
```

Scans over several flux surfaces use the helper script:

```bash
python python/run_driftorbit.py --exe ./build/neo_rt.x 0 5
```

## Documentation

- [Running NEO-RT](doc/running.md): inputs, namelist parameters, workflow.
- [Output file formats](doc/file_formats.md): columns and units.
- [Library interface](doc/library.md): `neort_lib` API.
- [Building](doc/building.md) and [diagnostic probes](doc/diagnostics.md).
- Physics notes in `doc/`, e.g. [Krook model](doc/krook.md),
  [species](doc/species.md), [collisional layer](doc/collisional_layer.md),
  [offset rotation](doc/offset_rotation.md).

FORD API documentation is generated with `make doc` (requires
`pip install -r requirements-docs.txt` and Graphviz) and published on GitHub
Pages from `main`.

## Testing

```bash
make test
```

runs the CTest suite, the golden-record check and the pytest-based tests. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the review and merge workflow.

## Related projects

- [libneo](https://github.com/itpplasma/libneo): Boozer field reading and
  shared utilities.
- [NEO-2](https://github.com/itpplasma/NEO-2): drift-kinetic neoclassical
  transport code; source of Boozer and profile inputs.

## Citing

C. G. Albert, M. F. Heyn, S. V. Kasilov, W. Kernbichler, A. F. Martitsch and
A. M. Runov, "Evaluation of toroidal torque by non-resonant magnetic
perturbations in tokamaks for resonant transport regimes using a Hamiltonian
approach", arXiv:1607.04665 (2016).

C. G. Albert, "Hamiltonian Theory of Resonant Transport Regimes in Tokamaks
with Perturbed Axisymmetry", Verlag der TU Graz (2020),
https://doi.org/10.3217/978-3-85125-746-5.

## License

MIT, see [LICENSE](LICENSE).
