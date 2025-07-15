# POTATO Integration Status Report

## Current Issue: Coordinate System Mismatch

The main problem is that I incorrectly assumed direct physical coordinates when NEO-RT uses Boozer flux coordinates.

### The Problem

In `src/potato_field_bridge.f90:365`, the coordinate conversion is wrong:
```fortran
call convert_neort_to_potato(v, eta, 1.7d0, 0.0d0, 0.0d0, z_eqm, success)
```

This uses fixed R=1.7m as a starting position, but NEO-RT works with normalized Boozer flux coordinates where:
- **s ∈ [0,1]**: Flux surface coordinate (s=0 is magnetic axis, s=1 is separatrix)
- **θ**: Poloidal angle in Boozer coordinates
- **φ**: Toroidal angle

### Coordinate System Differences

| System | Coordinates | Range | Description |
|--------|-------------|-------|-------------|
| NEO-RT | (s, θ, φ) | s ∈ [0,1] | Boozer flux coordinates |
| POTATO | (R, Z, φ) | R ∈ [R_min, R_max] | Physical cylindrical coordinates |

### Required Fix

Need to implement proper coordinate conversion:
1. Take flux surface s from NEO-RT (e.g., s=0.5 for mid-radius)
2. Convert (s, θ) → (R, Z) using equilibrium magnetic field data
3. Ensure (R, Z) is inside plasma separatrix for POTATO integration

## Files Created/Modified

### Core POTATO Integration
- **`src/potato_field_bridge.f90`**: Main interface between NEO-RT and POTATO
  - Status: Needs coordinate conversion fix
  - Current issue: Uses fixed R=1.7m instead of flux coordinates

### Real POTATO Integration  
- **`POTATO/SRC/velo_safe.f90`**: Safe velocity routine with FPE protection
  - Status: Working - integrated into build system
  - Purpose: Prevents floating point exceptions in POTATO orbit integration

### Physics Implementation
- **`src/thick_orbit_drift.f90`**: Complete thick orbit transport calculations
  - Status: Physics complete, needs proper coordinate inputs
  - Features: Drift velocities, perturbed Hamiltonian, transport coefficients

### Test and Example Files
- **`src/test_potato_real_integration.f90`**: Test for real POTATO functionality
  - Status: Needs coordinate conversion updates
  
- **`examples/thick_orbit_example.f90`**: Working example with plots
  - Status: Working with test coordinates

### Configuration Files
- **`field_divB0.inp`**: POTATO configuration for EFIT integration
  - Status: Working - successfully reads ASDEX Upgrade EFIT data
  - Path: Uses `/temp/ert/data/AUG/EQDSK/g30835.3200_ed6`

- **`convexwall.dat`**: ASDEX Upgrade vessel geometry
  - Status: Working - downloaded from libneo repository
  - Size: 252 R,Z coordinate pairs defining vessel boundary

### Build System
- **`Makefile`**: Updated to support thick orbit compilation
  - Status: Working - `make CONFIG=Debug USE_THICK_ORBITS=ON`

- **`CMakeLists.txt`**: Conditional compilation for thick orbits
  - Status: Working - properly links POTATO library

- **`POTATO/SRC/CMakeLists.txt`**: Added velo_safe.f90 to POTATO build
  - Status: Working

## Current Work: Coordinate Conversion Implementation ✅ COMPLETED

### What Was Done
1. **✅ Fixed convert_neort_to_potato()** in `src/potato_field_bridge.f90:276`:
   - Now accepts Boozer coordinates (s_flux, theta_boozer, phi_boozer) instead of fixed R,Z
   - Uses `booz_to_cyl()` from `do_magfie_mod` to convert (s, θ, φ) → (R, Z, φ)
   - Properly handles unit conversion (cm → m) between magfie and POTATO
   - Validates flux surface coordinate s ∈ [0,1]

2. **✅ Updated coordinate conversion calls**:
   - Now uses s=0.5 (mid-radius) as test flux surface
   - Proper Boozer coordinate system integration
   - Removed hardcoded R=1.7m physical coordinates

3. **✅ Enhanced do_magfie_standalone.f90**:
   - Fixed phi handling in `booz_to_cyl()` to use proper lambda function transformation
   - Ready for realistic coordinate conversion with ASDEX Upgrade equilibrium

### Next Steps
1. **Test with realistic coordinates**: Verify POTATO integration with s=0.5 flux surface
2. **Validate physics**: Compare thick vs thin orbit results for finite Larmor radius effects
3. **Production integration**: Use actual NEO-RT flux surface coordinates in production runs

## Integration Status Summary

| Component | Status | Description |
|-----------|--------|-------------|
| ✅ POTATO Build | Complete | Library compiles and links successfully |
| ✅ Field Bridge | Complete | Interface layer for magnetic field data |
| ✅ Safety Layer | Complete | velo_safe prevents floating point exceptions |
| ✅ Physics Calculations | Complete | Drift velocities and transport coefficients |
| ✅ EFIT Integration | Complete | Reads ASDEX Upgrade equilibrium data |
| ✅ Coordinate Conversion | **Complete** | **Boozer to cylindrical conversion implemented** |
| 🔄 Real Integration Test | **In Progress** | Testing with s=0.5 flux surface |
| ⏳ Production Validation | Pending | Full physics validation |

## Next Steps
1. Test POTATO integration with realistic Boozer coordinates (s=0.5)
2. Validate thick orbit physics results vs thin orbit approximation  
3. Integrate with NEO-RT production workflow for realistic transport calculations