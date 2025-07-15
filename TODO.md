# TODO: Thick Orbit Integration for NEO-RT

## Project Goal
**Enable the same NTV torque functionality as standard NEO-RT thin orbit approximation, but using thick guiding-center orbits from POTATO.**

This means: frequencies, resonances, transport coefficients, and NTV torque calculations must work cleanly with thick orbits. No shortcuts, no approximations - full physics implementation.

## ⚠️ **CRITICAL REALITY CHECK: Current Implementation Status**

### ✅ **COMPLETED: Infrastructure Only**
- Build system integration (CMake, linking, dependencies)
- POTATO field bridge framework (structure only)
- Test framework with TDD methodology
- Basic coordinate conversion utilities
- Stub interface for development

### ❌ **INCOMPLETE: Core Physics Implementation**
**After comprehensive code audit, ALL core thick orbit physics is still using stubs/approximations:**

1. **`src/potato_stub.f90`** - Entire stub module returns fake values (**BLOCKING**)
2. **`src/potato_wrapper.f90`** - Still calls stub instead of real POTATO (**BLOCKING**)
3. **`src/potato_field_bridge.f90`** - Real `find_bounce` call commented out (**BLOCKING**)
4. **`src/thick_orbit_drift.f90`** - All calculations use simplified estimates (**BLOCKING**)
5. **`src/transport_thick.f90`** - Falls back to thin orbit approximation (**BLOCKING**)
6. **EFIT field integration** - Not implemented (**BLOCKING**)

### 🚨 **DOCUMENTATION WAS INCORRECT** 
Previous documentation claimed "Real POTATO integration complete" but code audit reveals:
- NO real POTATO find_bounce calls
- NO real bounce-averaged physics
- NO realistic field initialization
- ALL physics calculations use stubs or thin-orbit approximations

## 🎯 **ACTUAL CURRENT PRIORITY: Implement Real POTATO Integration**

### 1. Replace Stub with Real POTATO Integration (**CRITICAL**)
- [ ] **Write failing test** for real POTATO find_bounce in `test/test_potato_find_bounce_real.f90`
- [ ] **Remove `src/potato_stub.f90`** - Replace with real POTATO calls
- [ ] **Fix `src/potato_wrapper.f90`** - Uncomment real POTATO calls (line 74)
- [ ] **Fix `src/potato_field_bridge.f90`** - Uncomment real find_bounce call (line 361)
- [ ] **Implement EFIT field initialization** - Required for realistic magnetic equilibrium
- [ ] **Test real POTATO integration** - Validate against stub behavior for development

### 2. Implement Real Thick Orbit Physics (**CRITICAL**)
- [ ] **Write failing test** for real bounce averaging in `test/test_real_bounce_averaging.f90`
- [ ] **Fix `src/thick_orbit_drift.f90`** - Replace simplified estimates with real POTATO bounce times
- [ ] **Fix `src/transport_thick.f90`** - Remove thin orbit approximation fallback
- [ ] **Fix `src/freq_thick.f90`** - Connect to real POTATO instead of stub
- [ ] **Test real physics** - Validate thick vs thin orbit differences are physical

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