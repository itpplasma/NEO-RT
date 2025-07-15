# POTATO Integration Strategy

## Project Goal
**Enable the same NTV torque functionality as standard NEO-RT thin orbit approximation, but using thick guiding-center orbits from POTATO.**

This means: frequencies, resonances, transport coefficients, and NTV torque calculations must work cleanly with thick orbits. No shortcuts, no approximations - full physics implementation.

## ⚠️ **CRITICAL STATUS UPDATE: Infrastructure Only**

### ✅ **COMPLETED: Framework and Build System**
- **CMake/Makefile**: Conditional compilation with `USE_THICK_ORBITS=ON`
- **Dependency Resolution**: POTATO library builds and links with NEO-RT
- **VODE Conflict**: Resolved linking conflicts between POTATO and NEO-RT VODE
- **Test Framework**: Comprehensive TDD test structure in place

### ❌ **INCOMPLETE: Core Physics Implementation**
**After comprehensive code audit, the following components are NOT implemented:**

### Coordinate System Bridge - **PARTIALLY COMPLETE**
- **`src/potato_field_bridge.f90`**: Framework exists but real `find_bounce` call commented out (line 361)
- **Boozer ↔ Cylindrical**: Basic conversion exists but uses hardcoded values
- **Field Interface**: Structure exists but falls back to stub implementations
- **Unit Conversion**: Framework present but not fully integrated

### Real POTATO Integration - **STUB ONLY**
- **`POTATO/SRC/velo_safe.f90`**: Safe velocity routine exists
- **Orbit Integration**: **STUB ONLY** - `src/potato_wrapper.f90` still calls stub (line 74)
- **Error Handling**: Framework exists but not connected to real POTATO
- **find_bounce calls**: **COMMENTED OUT** - All real POTATO calls are disabled

### Physics Implementation - **APPROXIMATIONS ONLY**
- **`src/thick_orbit_drift.f90`**: Uses simplified estimates, not real POTATO physics
- **Real Bounce Averaging**: **STUB** - No actual bounce averaging implementation
- **Transport Coefficients**: **THIN ORBIT FALLBACK** - Falls back to thin orbit approximation
- **Perturbed Hamiltonian**: Uses simplified formulas, not real orbit integration

### Configuration and Testing - **STUB CONFIGURATION**
- **`field_divB0.inp`**: ASDEX Upgrade EFIT configuration exists
- **`convexwall.dat`**: Vessel geometry integration exists
- **Test Suite**: Tests exist but expect stub behavior
- **Examples**: **STUB MODE** - All examples use stub implementations

## 🎯 **CURRENT PRIORITY: Complete NTV Workflow**

### Immediate Next Steps - **REALITY-BASED PRIORITIES**

1. **UNCOMMENT REAL POTATO CALLS**: Fix `src/potato_wrapper.f90` line 74 and `src/potato_field_bridge.f90` line 361
2. **REMOVE STUB DEPENDENCIES**: Replace `src/potato_stub.f90` with real POTATO integration
3. **IMPLEMENT EFIT FIELD SETUP**: Complete realistic magnetic field initialization
4. **FIX PHYSICS CALCULATIONS**: Replace simplified estimates with real POTATO physics in all modules
5. **VALIDATE REAL INTEGRATION**: Test that thick orbit results differ meaningfully from thin orbit

### Key Technical Insights

**Coordinate Systems**:
- NEO-RT: (s, θ, φ) Boozer flux coordinates, s ∈ [0,1]
- POTATO: (R, Z, φ) physical cylindrical coordinates
- Conversion: `booz_to_cyl()` from `do_magfie_mod`

**Performance Limitations**:
- Thick orbits ~10x slower than thin orbits (no spline optimization)
- Direct integration required for each (v,η) point
- Frequency calculations cannot use existing spline infrastructure

**Physics Validation**:
- Realistic gyroradius: ρ_gyro ~ 1-5 mm for thermal ions
- Finite orbit width effects: Δω/ω ~ (δr/L_B)²
- Resonance shifts due to altered bounce frequencies

## Current Status Summary - **CORRECTED**

| Component | Status | Next Action |
|-----------|--------|-------------|
| ✅ Build System | Complete | Maintain |
| ❌ Field Bridge | **STUB ONLY** | **Uncomment real find_bounce calls** |
| ❌ POTATO Integration | **STUB ONLY** | **Remove potato_stub.f90, fix wrapper** |
| ❌ Physics Calculations | **APPROXIMATIONS** | **Use real POTATO physics** |
| ❌ EFIT Integration | **INCOMPLETE** | **Implement field initialization** |
| ❌ Frequency Module | **STUB BASED** | **Connect to real POTATO** |
| ❌ Resonance Integration | **MISSING** | **Implement with real frequencies** |
| ❌ NTV Torque | **MISSING** | **Full workflow with real physics** |
| ❌ Production Validation | **MISSING** | **Real physics validation** |

## Strategic Approach - **UPDATED**

**Reality Check Policy**: Current implementation is mostly stubs - must be completely rewritten
**No Shortcuts Policy**: Every component must use real thick orbit physics, not approximations
**TDD Methodology**: All implementation must start with failing tests
**Production Ready**: Code must work with realistic EFIT data and production workflows
**Performance Awareness**: Document computational cost but prioritize physics correctness

**CRITICAL**: The foundation framework exists but ALL physics implementation is incomplete. The priority is implementing real POTATO integration, not building additional features on top of stub implementations.