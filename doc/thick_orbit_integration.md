# Thick Orbit Integration Architecture

## Overview

NEO-RT now supports both thin and thick orbit calculations through a runtime dispatch architecture. This allows seamless switching between the traditional thin orbit approximation and the more accurate thick orbit calculations from POTATO.

## Architecture Design

### Abstract Interface (`orbit_types.f90`)

```fortran
type, abstract :: orbit_type_t
contains
    procedure(calculate_bounce_time_interface), deferred :: calculate_bounce_time
    procedure(calculate_frequencies_interface), deferred :: calculate_frequencies
end type orbit_type_t
```

### Runtime Dispatch

Two concrete implementations extend the abstract interface:

1. **`thin_orbit_type_t`**: Wraps existing NEO-RT thin orbit calculations
2. **`thick_orbit_type_t`**: Interfaces with POTATO for thick orbit calculations

### Key Benefits

- **No preprocessor complexity**: Both implementations always available
- **Runtime selection**: Choose orbit type based on physics requirements
- **Easy comparison**: Switch between methods for validation
- **Clean interface**: Same API for both thick and thin orbits

## Usage Example

```fortran
class(orbit_type_t), allocatable :: orbit
real(8) :: taub, bounceavg(nvar)
real(8) :: omega_theta, omega_phi

! Select orbit type at runtime
if (use_thick_orbits) then
    allocate(thick_orbit_type_t :: orbit)
else
    allocate(thin_orbit_type_t :: orbit)
endif

! Use same interface for both
call orbit%calculate_bounce_time(v, eta, taub, bounceavg)
call orbit%calculate_frequencies(eta, omega_theta, omega_phi)
```

## Implementation Status

### Completed

- ✅ Abstract orbit interface (`orbit_type_t`)
- ✅ Thin orbit wrapper implementation
- ✅ Thick orbit interface with POTATO stub
- ✅ Runtime dispatch architecture
- ✅ Comprehensive test framework
- ✅ Comparison test utilities

### In Progress

- 🔄 Replace POTATO stub with actual POTATO integration
- 🔄 Full physics validation tests

## Testing Framework

The test suite includes:

1. **`test_potato_simple.x`**: Basic interface and runtime dispatch testing
2. **`test_orbit_comparison.x`**: Framework for comparing thick vs thin orbits
3. **`test_potato_bounce.x`**: Bounce time calculation tests
4. **`test_potato_frequencies.x`**: Frequency calculation tests

## Important Limitations

### Spline Scaling Invalid for Thick Orbits

The NEO-RT frequency optimization using splines assumes thin orbit scaling with velocity `v`, which breaks down for finite orbit width effects in thick orbits:

- Thick orbits sample different magnetic field regions
- Complex velocity-dependent drift physics invalidate simple scaling
- Must use direct bounce integration for each (v,η) point
- No frequency spline optimization possible for thick orbits

### Performance Considerations

- Thick orbit calculations significantly more expensive than thin orbits
- Each frequency evaluation requires full bounce integration
- Consider hybrid approach: thin orbits for bulk, thick for resonances

## Future Work

1. **POTATO Integration**: Replace stub with actual POTATO `find_bounce()` 
2. **Performance Optimization**: Investigate caching strategies for thick orbits
3. **Hybrid Mode**: Automatic selection based on orbit width criteria
4. **Validation Suite**: Comprehensive physics tests against POTATO standalone