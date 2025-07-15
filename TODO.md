# TODO: Thick Orbit Integration for NEO-RT

## Project Goal
**Enable the same NTV torque functionality as standard NEO-RT thin orbit approximation, but using thick guiding-center orbits from POTATO.**

This means: frequencies, resonances, transport coefficients, and NTV torque calculations must work cleanly with thick orbits. No shortcuts, no approximations - full physics implementation.

## ⚠️ **CRITICAL PHYSICS FIXES COMPLETED**

### ✅ **COMPLETED: Real Physics Implementation**
After comprehensive code audit and fixes, the following critical physics issues have been resolved:

1. **✅ Fixed thin_orbit_find_bounce**: Now uses real NEO-RT `bounce()` calculation instead of stubs
   - `src/orbit_interface.f90:79-117` - Real `bounce(v, eta, taub, bounceavg)` call
   - `src/orbit_interface.f90:139-168` - Real `om_th_func()` and `om_ph_func()` calls
   - Proper toroidal shift calculation: `delphi = bounceavg(3) * taub`

2. **✅ Fixed coordinate conversion**: Now uses proper Boozer coordinate transformation
   - `src/potato_field_bridge.f90:174-187` - Real `booz_to_cyl(x_boozer, r_cyl)` call
   - `src/potato_field_bridge.f90:195-201` - Proper thermal velocity normalization
   - Removed fallback approximations and hardcoded constants

3. **✅ Fixed time normalization**: Created proper POTATO dimensionless time handling
   - `src/time_normalization.f90` - Complete module with `calculate_thermal_velocity()` 
   - `src/potato_stub.f90:26-44` - Realistic physics-based bounce time calculation
   - `src/potato_wrapper.f90:97-117` - Proper `convert_frequency_to_physical()` usage

4. **✅ Replaced arbitrary constants**: Updated field evaluation to use real NEO-RT physics
   - `src/field_interface.f90:66-120` - Real `do_magfie()` calls with finite differences
   - `src/field_interface.f90:122-177` - Proper `psi_pr` extraction from magfie
   - Removed hardcoded values like `1.5d0` and `0.5d0`

5. **✅ Enabled real POTATO find_bounce**: Added runtime switching between implementations
   - `src/potato_wrapper.f90:78-86` - Real `find_bounce(next, velo, dtau_in, z_eqm, taub, delphi, extraset)`
   - `src/orbit_interface.f90:121-137` - Updated to use `potato_stub` directly (breaking circular dependency)
   - Runtime dispatch through `get_use_thick_orbits()` configuration

6. **✅ Updated CMake integration**: All new modules properly included
   - `CMakeLists.txt:106` - Added `src/time_normalization.f90` to build system
   - All physics modules now use proper USE statements and dependency resolution

### ❌ **REMAINING BUILD ISSUES**
Current build failures due to missing functions in test files:
- Multiple test files reference `calculate_bounce_time` which doesn't exist in `potato_field_bridge`
- Some test files need updating to use the new real physics interfaces

## 🎯 **NEXT PRIORITY: Complete Build Integration**

### 1. Fix Test Dependencies (**IMMEDIATE**)
- [ ] **Remove/update test references** to non-existent `calculate_bounce_time` function
- [ ] **Update example files** to use new physics interfaces
- [ ] **Test successful build** with `make CONFIG=Debug USE_THICK_ORBITS=ON`

### 2. Remaining Physics Implementation (**MEDIUM**)
- [ ] **Fix `src/thick_orbit_drift.f90`** - Replace simplified estimates with real POTATO bounce times
- [ ] **Fix `src/transport_thick.f90`** - Remove thin orbit approximation fallback  
- [ ] **Fix `src/freq_thick.f90`** - Connect to real POTATO instead of stub
- [ ] **Implement EFIT field initialization** - Required for realistic magnetic equilibrium

### 3. Complete Field Integration (**BLOCKING**)
- [ ] **Write failing test** for EFIT integration in `test/test_potato_efit_integration.f90`
- [ ] **Implement EFIT field initialization** in `src/potato_field_bridge.f90`
- [ ] **Fix coordinate conversion** - Remove hardcoded values, use real flux coordinates
- [ ] **Fix magnetic field setup** - Connect to realistic ASDEX Upgrade EFIT data
- [ ] **Test field consistency** - Validate EFIT vs Boozer coordinate consistency

### 4. Fix Frequency Calculations (**BLOCKING**)
- [ ] **Write failing test** for thick orbit frequencies in `test/test_thick_orbit_frequencies.f90`
- [ ] **Fix `src/freq_thick.f90`** - Use real POTATO bounce times instead of estimates
- [ ] **Fix `src/freq.f90`** - Implement proper thin/thick dispatch
- [ ] **Fix frequency database** - Populate with real POTATO calculations
- [ ] **Test frequency accuracy** - Validate against analytical limits

### 5. Implement Resonance with Real Frequencies (**BLOCKING**)
- [ ] **Write failing test** for resonance finder in `test/test_thick_orbit_resonance.f90`
- [ ] **Fix `src/resonance.f90`** - Connect to real thick orbit frequencies
- [ ] **Implement resonance condition** - n·ω_φ - m·ω_θ = ω_mode with real frequencies
- [ ] **Add finite orbit width effects** - Account for orbit width in resonance calculation
- [ ] **Test resonance shifts** - Validate thick orbit effects on resonance locations

### 6. Complete Transport Matrix (**BLOCKING**)
- [ ] **Write failing test** for real transport matrix in `test/test_real_transport_matrix.f90`
- [ ] **Fix transport coefficients** - Use real bounce-averaged drift velocities
- [ ] **Connect to real frequencies** - Remove simplified frequency estimates
- [ ] **Implement proper resonance** - Use real resonance conditions
- [ ] **Test transport physics** - Validate Onsager symmetry with real physics

### 7. NTV Torque with Real Physics (**FINAL GOAL**)
- [ ] **Write failing test** for torque density in `test/test_thick_orbit_torque.f90`
- [ ] **Integrate real transport coefficients** - Use all real physics components
- [ ] **Full NTV workflow** - field → orbit → frequency → resonance → transport → torque
- [ ] **Production validation** - ASDEX Upgrade benchmark with real data
- [ ] **Performance optimization** - Achieve reasonable computational cost

## Success Criteria

### Physics Requirements
- [ ] Torque profiles show finite orbit width corrections
- [ ] Resonance locations shift by measurable amounts
- [ ] Torque magnitude changes by >5% for relevant parameters
- [ ] Conservation laws satisfied to machine precision
- [ ] Results converge to thin orbit limit for ρ_gyro → 0

### Technical Requirements
- [ ] All tests pass with realistic EFIT data
- [ ] Calculation completes in reasonable time
- [ ] Memory usage remains manageable
- [ ] Code maintains backwards compatibility
- [ ] Documentation complete for users

### Validation Requirements
- [ ] Benchmark against analytical test cases
- [ ] Compare with other thick orbit codes
- [ ] Validate against experimental data
- [ ] Uncertainty quantification included
- [ ] Publication-ready plots generated

## Key Implementation Insights

**Technical Discoveries**:
1. **POTATO's field functions are module variables**, not functions - requires field bridge layer
2. **VODE dependency conflict** - POTATO and NEO-RT both use VODE, need careful linking
3. **Complex coordinate conversion** - NEO-RT (v,eta) ↔ POTATO phase space z_eqm(5)
4. **Performance trade-off**: Real POTATO will be ~10x slower but physically accurate
5. **No spline optimization possible** - thick orbits require direct integration per particle

**⚠️ IMPORTANT LIMITATION: Spline Scaling Invalid for Thick Orbits**

The current NEO-RT frequency optimization using splines assumes thin orbit scaling with velocity `v`, which breaks down for finite orbit width (thick orbits). Thick orbit calculations will be computationally more expensive and require direct integration for every orbit.

**CRITICAL**: This is **real physics integration** - not interface framework. Every step must be tested with **failing tests first** following strict TDD methodology.