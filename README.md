# NEO-RT
Calculates NTV torque in resonant transport regimes using a Hamiltonian approach

## Dependencies
NEO-RT requires the GNU Fortran compiler, make, cmake, ninja, Python with NumPy, 
BLAS/LAPACK, SuiteSparse, NetCDF. On Debian/Ubuntu they will be installed as `sudo`
via
```bash
make deps
```

## Building
Inside the NEO-RT folder run
```bash
make
```
This will create the main binary `neo_rt.x` amongst other files in the `build`
directory.

## Running

Single run with input file: `neo_rt.x <runname>` will run with input file `<runname>.in` .
An example input file `driftorbit.in` can be found in `examples/base`.
Note that the name of the input file is given without file ending, that
is because this is also the name (prefix) to use for output files.

Script to be called for batch runs on multiple flux surfaces: `run_driftorbit.py`

## Input
- `<runname>.in`: Given as first command line argument, contains run parameters 
- `in_file`: Boozer coordinate file of axisymmetric part of the magnetic field
- `in_file_pert`: Boozer coordinate file of non-axisymmetric perturbation
- `plasma.in`: Plasma parameters, must be **equidistant in radius** `s` ! Required for torque and/or nonlinear.
- `profile.in`: Profile parameters, required for derivatives used inside nonlinear runs.

Used by `run_driftorbit.py`:

- `driftorbit.in.template`: Contains placeholders to be filled by profile data
- `profile.in`: Radius `s`, toroidal mach number `Mt` and thermal velocity `vth` (the latter is not used anymore)

## Output

* `driftorbit_magfie_param.out`: magnetic field parameters
* `driftorbit_torque.out`: torque density data

## Postprocessing

* `examples/plot_torque.py` is the most recent plotting utility for toroidal torque
