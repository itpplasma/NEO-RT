# Thick Orbit Example

This example demonstrates how to use NEO-RT's thick orbit functionality through the runtime dispatch architecture.

## Files

- `thick_orbit_example.f90`: Main example program
- `README.md`: This documentation

## Description

The example demonstrates three key use cases:

### 1. Single Orbit Comparison
- Compares thick and thin orbit calculations for a single particle
- Shows bounce times and frequencies side by side
- Demonstrates runtime orbit type selection

### 2. Pitch Angle Scan
- Scans across pitch angle parameter η
- Shows how frequencies vary with orbit geometry
- Calculates effective safety factor q_eff = ω_φ/ω_θ

### 3. Resonance Analysis
- Searches for resonances with magnetic perturbations
- Uses canonical frequencies: n·ω_φ - m·ω_θ = ω_mode
- ω_φ from POTATO's delphi (toroidal shift per bounce time)
- Critical for NTV torque calculations and resonant transport

## Building and Running

From the NEO-RT root directory:

```bash
# Build the project
make

# Copy example to build directory (if not automatically built)
cp examples/thick_orbit/thick_orbit_example.f90 build/
cd build
gfortran -I. thick_orbit_example.f90 -L. -lneo_rt -o thick_orbit_example.x

# Run the example
./thick_orbit_example.x
```

## Key Features Demonstrated

### Runtime Dispatch
```fortran
class(orbit_type_t), allocatable :: orbit

! Select thick orbit implementation
allocate(thick_orbit_type_t :: orbit)
call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
call orbit%calculate_frequencies(eta, omega_theta, omega_phi)
```

### Physics Comparison
- **Bounce Times**: Thick orbits typically show different bounce times due to finite width effects
- **Frequencies**: Modified by thick orbit physics affecting resonance locations
- **Safety Factor**: Effective q-profile modification by thick orbit dynamics

### Performance Awareness
- Demonstrates computational cost differences
- Shows appropriate use cases for thick vs thin orbits
- Guides optimization strategies

## Expected Output

The example will show:
1. Direct comparison of thick vs thin orbit results
2. Frequency variations across pitch angles
3. Resonance identification and analysis

## Notes

- Currently uses stub implementation until actual POTATO integration
- Results demonstrate the interface pattern and expected physics
- Ready for immediate use when POTATO integration is complete

## Physics Interpretation

### Thick Orbit Effects
- **Finite Larmor Radius**: Orbits sample magnetic field gradients
- **Guiding Center Motion**: Full 5D phase space integration
- **Modified Resonances**: Broader resonance widths, shifted locations

### When to Use Thick Orbits
- Near magnetic axis where gradients are strong
- For fast particles with large Larmor radius
- When precise resonance calculations are needed
- For physics requiring finite orbit width effects