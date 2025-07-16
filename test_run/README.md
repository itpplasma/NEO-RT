# NEO-RT Working Directory Setup

This directory contains the proper working directory structure for running NEO-RT with real physics instead of synthetic test data.

## Required Files

The following files are required for NEO-RT to run properly:

### 1. `in_file` - Boozer Coordinate Magnetic Field Data
- **Source**: `/home/ert/code/NEO-RT/examples/base/in_file`
- **Format**: Boozer coordinate data file
- **Description**: Contains magnetic field data in Boozer coordinates
- **First line**: `CC Boozer-coordinate data file Version 0.1 Author J.Geiger Created: 26.07.2010`

### 2. `driftorbit.in` - NEO-RT Namelist Input Parameters
- **Source**: `/home/ert/code/NEO-RT/examples/base/driftorbit.in`
- **Format**: Fortran namelist format
- **Description**: Contains physics parameters for NEO-RT calculation
- **Key parameters**:
  - `s = 0.5`: Radial coordinate (normalized flux surface)
  - `vth = 40000000.0`: Thermal velocity [cm/s]
  - `comptorque = .true.`: Compute torque
  - `inp_swi = 8`: Input switch for Boozer file format

### 3. `plasma.in` - Plasma Profiles (Required for Torque Calculations)
- **Source**: `/home/ert/code/NEO-RT/examples/base/plasma.in`
- **Format**: Plasma profile data
- **Description**: Contains plasma density and temperature profiles
- **Required when**: `comptorque = .true.` in `driftorbit.in`

## Usage

### Running NEO-RT with Real Physics

Instead of using synthetic test data, this working directory allows running NEO-RT with:

1. **Real magnetic field geometry** from Boozer coordinates
2. **Realistic plasma parameters** from the namelist
3. **Actual plasma profiles** for torque calculations

### Test Verification

The working directory setup can be verified by running:

```bash
gfortran -o test_orbit_trajectory_comparison test_orbit_trajectory_comparison.f90
./test_orbit_trajectory_comparison
```

This test verifies that:
- All required input files are present
- Files are in the correct format
- NEO-RT can be run with real physics

### Comparison with Test Program

The original `test_orbit_trajectory_comparison.f90` program was designed to:

1. **Initialize real NEO-RT physics** using `do_magfie_init`
2. **Compare thin vs thick orbit calculations** using real magnetic field data
3. **Validate orbit trajectory calculations** with actual plasma parameters

This working directory provides the necessary input files for that initialization to work properly.

## Physics Validation

With this setup, the test program will:

1. **Load real magnetic field** from `in_file`
2. **Initialize plasma parameters** from `driftorbit.in`
3. **Use actual plasma profiles** from `plasma.in`
4. **Calculate real bounce times** and toroidal shifts
5. **Compare thin vs thick orbit physics** with realistic parameters

This represents a significant improvement over synthetic test data, as it validates the physics implementation against real tokamak geometry and plasma conditions.