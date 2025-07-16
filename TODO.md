# TODO: Thick Orbit Integration for NEO-RT

## Project Goal
**Enable the same NTV torque functionality as standard NEO-RT thin orbit approximation, but using thick guiding-center orbits from POTATO.**

This means: frequencies, resonances, transport coefficients, and NTV torque calculations must work cleanly with thick orbits. No shortcuts, no approximations - full physics implementation.

## ⚠️ **CRITICAL STATUS: FUNDAMENTAL IMPLEMENTATION GAPS**

### **🔴 BLOCKING ISSUES - NOTHING IS ACTUALLY IMPLEMENTED**

After honest assessment, the thick orbit integration is **NOT FUNCTIONAL**:

1. **🔴 POTATO Integration**: `src/potato_stub.f90` returns fake values
2. **🔴 Field Interface**: Real field evaluation fails with convergence errors
3. **🔴 Physics Calculations**: All "real" API calls fall back to approximations when they fail
4. **🔴 Thick Orbit Functions**: All thick orbit modules use stubs or thin orbit fallbacks
5. **🔴 Visualization**: Programs exist but cannot run due to physics initialization failures

**REALITY**: The thick orbit integration framework exists but contains NO WORKING PHYSICS.

## **MANDATORY IMPLEMENTATION ROADMAP**

### **Phase 1: Fix POTATO Integration Foundation** 
**STATUS: 🔴 BLOCKED - NOTHING WORKS**

#### **1.1 Fix POTATO Stub Implementation**
- [ ] **Remove `src/potato_stub.f90` entirely** - This file returns fake values
- [ ] **Implement real `find_bounce()` in `src/potato_wrapper.f90`** - Currently commented out (line 74)
- [ ] **Connect actual POTATO orbit integration** - No shortcuts, no approximations
- [ ] **Fix undefined symbol errors** - `phielec_of_psi_`, `denstemp_of_psi_`, `field_eq_`
- [ ] **Test POTATO library builds** - Verify all symbols resolve correctly

#### **1.2 Fix Field Interface** 
- [ ] **Fix field evaluation convergence** - Current magnetic field calculations fail
- [ ] **Implement proper EFIT reader** - Real equilibrium data, not synthetic
- [ ] **Fix `field_divB0.inp` initialization** - Proper field data loading
- [ ] **Validate field consistency** - Check ∇·B = 0 and flux surface alignment
- [ ] **Test field evaluation** - Verify `psif`, `dpsidr`, `dpsidz` return physical values

#### **1.3 Fix Physics Initialization**
- [ ] **Fix VODE solver convergence** - Eliminate numerical instabilities
- [ ] **Fix `init_done` initialization** - Physics modules fail to initialize
- [ ] **Fix bounce integral convergence** - Current bounce calculations fail
- [ ] **Test basic physics functions** - Verify `bounce()`, `Om_th()`, `Om_ph()` work
- [ ] **Validate against known results** - Compare with working thin orbit calculations

### **Phase 2: Implement Real Thick Orbit Physics**
**STATUS: 🔴 NOT STARTED - DEPENDS ON PHASE 1**

#### **2.1 Real Frequency Calculations**
- [ ] **Implement `src/freq_thick.f90`** - Real thick orbit frequencies, not stubs
- [ ] **Connect to real POTATO bounce times** - No approximations
- [ ] **Implement `compute_canonical_frequencies_thick()`** - Real physics, not fallbacks
- [ ] **Test frequency calculations** - Verify ω_θ and ω_φ differ from thin orbit
- [ ] **Validate frequency limits** - Check trapped/passing particle limits

#### **2.2 Real Resonance Analysis**
- [ ] **Implement `src/resonance.f90` thick orbit functions** - Real resonance finding
- [ ] **Fix `driftorbit_coarse_thick()`** - Real resonance broadening, not approximations
- [ ] **Implement orbit width parameter** - Real `calculate_orbit_width_parameter()`
- [ ] **Test resonance conditions** - Verify n*ω_φ - m*ω_θ = 0 finding works
- [ ] **Validate resonance broadening** - Check finite orbit width effects

#### **2.3 Real Transport Coefficients**
- [ ] **Implement `src/transport_thick.f90`** - Real transport matrix, not stubs
- [ ] **Connect to real thick orbit drift velocities** - No approximations
- [ ] **Implement `calculate_transport_coefficients_thick()`** - Real physics only
- [ ] **Test Onsager symmetry** - Verify D₁₂ = D₂₁ for real calculations
- [ ] **Validate transport matrix** - Check against analytical limits

#### **2.4 Real NTV Torque Calculation**
- [ ] **Implement `src/torque_thick.f90`** - Real torque density, not stubs
- [ ] **Connect to real velocity space integration** - No shortcuts
- [ ] **Implement `calculate_ntv_torque_density()`** - Real physics only
- [ ] **Test torque conservation** - Verify total torque is conserved
- [ ] **Validate against experiments** - Compare with measured NTV torque

### **Phase 3: Create Working Visualization**
**STATUS: 🔴 BLOCKED - DEPENDS ON PHASE 2**

#### **3.1 Fix Existing Fortran Programs**
- [ ] **Make `plot_bounce_time_real.f90` actually work** - Currently fails due to physics issues
- [ ] **Make `plot_canonical_frequencies_real.f90` actually work** - Currently uses fallbacks
- [ ] **Make `plot_orbit_trajectories_real.f90` actually work** - Currently fails initialization
- [ ] **Make `plot_resonance_analysis_real.f90` actually work** - Currently uses stubs
- [ ] **Make `plot_transport_coefficients_real.f90` actually work** - Currently uses stubs
- [ ] **Make `plot_torque_profiles_real.f90` actually work** - Currently uses stubs

#### **3.2 Generate Real Physics Plots**
- [ ] **Real bounce time comparison** - Show actual thick orbit bounce times
- [ ] **Real frequency comparison** - Show actual ω_θ and ω_φ differences
- [ ] **Real orbit trajectories** - Show actual thick orbit R-Z trajectories
- [ ] **Real resonance analysis** - Show actual resonance condition changes
- [ ] **Real transport matrix** - Show actual D_ij matrix elements
- [ ] **Real torque profiles** - Show actual NTV torque density profiles

### **Phase 4: Validation and Testing**
**STATUS: 🔴 NOT STARTED - DEPENDS ON PHASE 3**

#### **4.1 Physics Validation**
- [ ] **Test thin orbit limit** - Verify thick orbit → thin orbit as ρ* → 0
- [ ] **Test trapped particle limits** - Verify deeply trapped particle physics
- [ ] **Test passing particle limits** - Verify minimal corrections for passing particles
- [ ] **Test conservation laws** - Verify energy and magnetic moment conservation
- [ ] **Test symmetries** - Verify Onsager relations and other symmetries

#### **4.2 Benchmark Against Experiments**
- [ ] **ASDEX Upgrade benchmark** - Compare with experimental NTV torque
- [ ] **Compare with other codes** - Validate against other finite orbit codes
- [ ] **Performance analysis** - Document runtime and memory usage
- [ ] **Convergence studies** - Document numerical convergence properties

## **BLOCKING DEPENDENCIES**

### **Critical Path:**
1. **POTATO Integration** → **Field Interface** → **Physics Initialization**
2. **Physics Initialization** → **Real Thick Orbit Physics**
3. **Real Thick Orbit Physics** → **Working Visualization**
4. **Working Visualization** → **Validation and Testing**

### **Cannot Proceed Without:**
- [ ] **Functional POTATO library** - No stubs, no approximations
- [ ] **Working field evaluation** - No convergence failures
- [ ] **Stable physics initialization** - No VODE errors
- [ ] **Real thick orbit calculations** - No fallbacks to thin orbit

## **SUCCESS CRITERIA**

### **Minimum Viable Product**
- [ ] **All physics functions work** - No stubs, no approximations, no fallbacks
- [ ] **All visualization programs run** - No initialization failures
- [ ] **All plots show real physics** - Actual thick orbit vs thin orbit differences
- [ ] **Quantitative validation** - Results match expected physics

### **Performance Requirements**
- [ ] **Reasonable runtime** - <10x thin orbit calculation time
- [ ] **Memory usage** - Scales linearly with orbit count
- [ ] **Numerical stability** - Convergent for all relevant parameters
- [ ] **Backwards compatibility** - Does not break existing thin orbit functionality

## **CURRENT STATUS SUMMARY**

### **🔴 PHASE 1: BLOCKED**
- **POTATO Integration**: NOT IMPLEMENTED (stubs only)
- **Field Interface**: BROKEN (convergence failures)
- **Physics Initialization**: BROKEN (VODE errors)

### **🔴 PHASE 2: NOT STARTED**
- **Thick Orbit Physics**: NOT IMPLEMENTED (all modules use stubs/fallbacks)

### **🔴 PHASE 3: BLOCKED**
- **Visualization**: EXISTS BUT CANNOT RUN (physics initialization failures)

### **🔴 PHASE 4: NOT STARTED**
- **Validation**: IMPOSSIBLE (no working physics to validate)

## **IMMEDIATE NEXT STEPS**

1. **Fix POTATO stub implementation** - Remove all fake return values
2. **Fix field evaluation convergence** - Resolve numerical instabilities
3. **Fix physics initialization** - Eliminate VODE solver errors
4. **Test basic physics functions** - Verify bounce(), Om_th(), Om_ph() work
5. **Only then** can we proceed to implement real thick orbit physics

**BOTTOM LINE**: Currently, NOTHING works except synthetic examples. The entire thick orbit integration must be implemented from scratch with real physics, no shortcuts, no approximations.