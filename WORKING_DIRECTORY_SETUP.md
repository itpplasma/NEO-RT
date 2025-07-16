# NEO-RT Working Directory Setup for Real Physics Testing

## Overview

This document describes the proper working directory structure required for running NEO-RT with real physics instead of synthetic test data. The setup has been implemented and validated in the `/home/ert/code/NEO-RT/test_run/` directory.

## Problem Statement

The original `test_orbit_trajectory_comparison.f90` was designed to test orbit trajectory calculations but was using synthetic physics data. To enable testing with real NEO-RT physics, a proper working directory structure with the necessary input files is required.

## Solution

### Working Directory Structure

The working directory must contain the following files:

```
test_run/
├── in_file                           # Boozer coordinate magnetic field data
├── driftorbit.in                     # NEO-RT namelist input parameters  
├── plasma.in                         # Plasma profiles for torque calculations
├── test_orbit_trajectory_comparison.f90  # Test program
└── setup_working_dir.sh             # Setup script for future use
```

### Required Files

#### 1. `in_file` - Boozer Coordinate Magnetic Field Data
- **Source**: `/home/ert/code/NEO-RT/in_file` (root directory)
- **Format**: Boozer coordinate data file
- **Size**: 109,758 bytes
- **Description**: Contains magnetic field geometry data for tokamak configuration
- **Key parameters**:
  - Configuration: tok02 (aspect ratio ~3.8)
  - 63 flux surfaces
  - R₀ = 1.64377 m, a = 0.46 m

#### 2. `driftorbit.in` - NEO-RT Namelist Input Parameters
- **Source**: `/home/ert/code/NEO-RT/examples/base/driftorbit.in`
- **Format**: Fortran namelist `&params`
- **Size**: 1,135 bytes
- **Key parameters**:
  - `s = 0.5`: Radial coordinate (normalized flux surface)
  - `vth = 40000000.0`: Thermal velocity [cm/s]
  - `comptorque = .true.`: Enable torque calculations
  - `inp_swi = 8`: Boozer file format selector

#### 3. `plasma.in` - Plasma Profiles
- **Source**: `/home/ert/code/NEO-RT/examples/base/plasma.in`
- **Format**: Tabulated plasma profiles
- **Size**: 15,998 bytes
- **Content**: 104 flux surfaces with density and temperature profiles
- **Required for**: Torque calculations when `comptorque = .true.`

### Validation

The working directory setup has been validated with a test program that verifies:

1. **File existence**: All required files are present
2. **File format**: Files are in the correct format
3. **Content verification**: Files contain expected data structure

### Usage

#### Setting Up a New Working Directory

1. **Use the setup script**:
   ```bash
   ./setup_working_dir.sh
   ```

2. **Manual setup**:
   ```bash
   cp /home/ert/code/NEO-RT/in_file .
   cp /home/ert/code/NEO-RT/examples/base/driftorbit.in .
   cp /home/ert/code/NEO-RT/examples/base/plasma.in .
   ```

#### Verifying the Setup

```bash
gfortran -o test_orbit_trajectory_comparison test_orbit_trajectory_comparison.f90
./test_orbit_trajectory_comparison
```

Expected output:
```
========================================
Test: Working Directory Setup
========================================

Test 1: Check for in_file (Boozer coordinate file)
  ✓ in_file exists
Test 2: Check for driftorbit.in (namelist input file)
  ✓ driftorbit.in exists
Test 3: Check for plasma.in (plasma profiles file)
  ✓ plasma.in exists
Test 4: Verify in_file format
  ✓ in_file appears to be Boozer coordinate file
Test 5: Verify driftorbit.in namelist format
  ✓ driftorbit.in contains valid namelist

Working Directory Setup Summary:
  Total tests: 5
  Passed: 5
  Failed: 0
  ✓ Working directory is properly set up for NEO-RT!
```

## Impact on Testing

With this working directory structure:

1. **Real Physics Initialization**: The test can now call `do_magfie_init` with real magnetic field data
2. **Realistic Parameters**: Physics calculations use actual tokamak geometry and plasma parameters
3. **Meaningful Comparisons**: Thin vs thick orbit comparisons are performed with real data
4. **Validation**: Results can be compared against known physics expectations

## Files Created

The implementation has created the following files in `/home/ert/code/NEO-RT/test_run/`:

- `test_orbit_trajectory_comparison.f90` - Working directory validation test
- `README.md` - Detailed documentation of the setup
- `setup_working_dir.sh` - Script for setting up new working directories
- All required input files (`in_file`, `driftorbit.in`, `plasma.in`)

## Next Steps

This working directory setup enables:

1. **Real Physics Testing**: The `test_orbit_trajectory_comparison` program can now run with real NEO-RT physics
2. **Development Validation**: New features can be tested against realistic tokamak conditions
3. **Physics Benchmarking**: Results can be compared with established NEO-RT calculations
4. **Debugging**: Physics issues can be diagnosed with real data instead of synthetic tests

The working directory structure provides a solid foundation for transitioning from synthetic test data to real physics validation in NEO-RT development.