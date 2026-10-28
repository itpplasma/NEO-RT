---
title: Building NEO-RT
---

# Building NEO-RT

NEO-RT is built with CMake and Ninja through a thin `Makefile` wrapper.

## Dependencies

- gfortran and a C/C++ toolchain.
- CMake 3.22 or newer and Ninja.
- BLAS and LAPACK. On GNU/Linux with gfortran, OpenBLAS is preferred when
  pkg-config finds it, so that a stale system MKL is not picked up by accident.
  Set `BLA_VENDOR` to override.
- Python 3 with NumPy for the regression tests and helper scripts.
- Optional: [FORD](https://github.com/Fortran-FOSS-Programmers/ford) and
  Graphviz for the API documentation.

On Debian and Ubuntu, `make deps` installs the system packages listed in the
`deps-debian` target of the `Makefile`.

CMake fetches the remaining Fortran dependencies at configure time:
[libneo](https://github.com/itpplasma/libneo) (Boozer field reader and the
Fortio-backed HDF5 interface), fortnum, fortplot (used by `neo_rt_diag.x`) and
the spline library. Network access is therefore needed for the first
configuration.

## Build configurations

```bash
make                 # Release build in build/
make CONFIG=Debug
make CONFIG=Fast     # aggressive optimisation
make clean           # needed before switching configuration
make reconfigure     # rerun CMake in the existing build directory
```

The build produces `build/neo_rt.x`, `build/neo_rt_diag.x`, the probe
executables described in [diagnostics.md](diagnostics.md), and static
libraries for embedding.

## Selecting the libneo revision

By default CMake uses the libneo commit pinned in `CMakeLists.txt`. To build
against another branch, tag or commit:

```bash
make LIBNEO_REF=my-branch
# or
cmake -S . -B build -DLIBNEO_REF=my-branch
```

Only a value given on the `make` command line is forwarded; a `LIBNEO_REF`
variable in the shell environment is ignored.

## Magnetic-field backend

NEO-RT uses its standalone Boozer field reader (`USE_STANDALONE=ON`, the
default). The alternative backend that delegates field evaluation to
[NEO-2](https://github.com/itpplasma/NEO-2) (`-DUSE_STANDALONE=OFF`) is
currently not supported: configuration stops with an error, because the
adapter in `src/do_magfie_neo.f90` needs a NEO-2 field and perturbation API
that is not available on NEO-2 `main`.

## Python environment

The helper scripts, examples and Python tests use NumPy, SciPy, Matplotlib,
`f90nml` and `h5py`. On systems with an externally managed Python, create a
repository-local virtual environment:

```bash
./setup-venv.sh
source .venv/bin/activate
```

With `.venv` active, `make test` uses this interpreter for the pytest-based
checks.

## API documentation

```bash
pip install -r requirements-docs.txt
make doc
```

This runs `ford doc.md` and writes HTML to `build/doc` without compiling the
solver. FORD documents `src/` with the standalone field backend;
`src/do_magfie_neo.f90` is excluded because it defines the same module names.
