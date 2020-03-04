# NEO-RT
Calculates NTV torque in resonant transport regimes using a Hamiltonian approach

## Building
* Clone latest NEO-RT version from github
```bash
git clone https://github.com/itpplasma/NEO-RT
```
* Prerequisites: In the directory where the NEO-RT folder is located, run
```bash
git clone https://github.com/itpplasma/spline
mkdir contrib
cd contrib
git clone https://github.com/itpplasma/quadpack
git clone https://github.com/itpplasma/vode
```
This will download a custom spline library as well as copies of quadpack and vode
in Fortran90 and with CMake support.

* Inside the NEO-RT folder run
```bash
mkdir build
cmake ..
make
```
This will create the main binary `driftorbit_test` amongst other files.
From there on, when changing code, one only needs to run `make` inside `build`.

## Running

Single run with input file: `driftorbit_test driftorbit.in`.
An example can be found in `test/base`.

Script to be called for batch runs on multiple flux surfaces: `run_driftorbit.py`

## Input
- `driftorbit.in`: Given as first command line argument, contains run parameters 
- `in_file`: Boozer coordinate file of axisymmetric part of the magnetic field
- `in_file_pert`: Boozer coordinate file of non-axisymmetric perturbation
- `plasma.in`: Plasma parameters, must be **equidistant in radius** `s` !

Used by `run_driftorbit.py`:

- `driftorbit.in.template`: Contains placeholders to be filled by profile data
- `profile.in`: Radius `s`, toroidal mach number `Mt` and thermal velocity `vth` (the latter is not used anymore)

## Output

* `driftorbit_magfie_param.out`: magnetic field parameters
* `driftorbit_torque.out`: torque density data

## Postprocessing

* `test/plot_torque.py` is the most recent plotting utility for toroidal torque
