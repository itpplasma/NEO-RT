# TODO: Thick Orbit Integration for NEO-RT

## Project Goal
Extend NEO-RT main code to support full guiding-center (thick) orbit calculations like in POTATO, with automated testing comparing bounce and precession frequencies between thin and thick orbit implementations.

## Phase 1: Interface Compatibility Layer (TDD Foundation)

### 1.1 Create Thick Orbit Abstract Interface
- [x] Write failing test for thick orbit interface in `test/test_thick_orbit_interface.f90`
- [x] Create abstract orbit type `orbit_type_t` in new `src/orbit_types.f90` module
- [x] Define common interface for bounce calculations: `calculate_bounce_time()`, `calculate_frequencies()`
- [x] Implement thin orbit wrapper using existing `bounce()` routine
- [x] Write basic integration test comparing thin orbit wrapper to direct calls

### 1.2 Magnetic Field Interface Extension  
- [x] Write failing test for thick orbit magnetic field interface in `test/test_magfie_thick.f90`
- [x] Create `magfie_thick_mod` module in `src/magfie_thick.f90`
- [x] ~~Implement finite gyroradius magnetic field evaluation~~ (REVISED: Not needed - POTATO uses standard magnetic field)
- [x] ~~Add gyroradius calculation routines from POTATO's approach~~ (REVISED: Thick orbit effects come from guiding-center equations)
- [x] Test magnetic field interface availability and basic functionality

**NOTE**: Phase 1.2 was initially misguided. POTATO doesn't modify the magnetic field interface - it uses the standard magnetic field but integrates full guiding-center equations of motion to naturally include finite Larmor radius effects.

### 1.3 Build System Integration
- [x] Write failing test for CMake thick orbit option in `test/test_build_thick_orbits.cmake`
- [x] Add `USE_THICK_ORBITS` CMake option to `CMakeLists.txt`
- [x] Create conditional compilation flags for thick orbit modules
- [ ] Integrate POTATO source files into build when thick orbits enabled
- [x] Test build works with both thick orbit modes (ON/OFF)

## Phase 2: POTATO Integration (Direct Interface)

**REVISED APPROACH**: Instead of porting POTATO code, create direct interface layer around existing POTATO functionality and add POTATO sources to CMake build system.

### 2.1 POTATO Source Integration  
- [ ] Write failing test for POTATO source availability in `test/test_potato_sources.f90`
- [ ] Add POTATO sources to CMake when `USE_THICK_ORBITS=ON`
- [ ] Create conditional build for `POTATO/SRC/*.f90` files
- [ ] Add POTATO include directories and dependencies
- [ ] Test POTATO library builds successfully with NEO-RT

### 2.2 POTATO Interface Layer
- [ ] Write failing test for POTATO interface in `test/test_potato_interface.f90` 
- [ ] Create `src/potato_interface.f90` wrapper module
- [ ] Implement `thick_orbit_type_t` that calls POTATO's `find_bounce()` directly
- [ ] Add interface to POTATO's `velo()` for guiding-center velocity field
- [ ] Map POTATO's 5D phase space (R, φ, Z, p, λ) to NEO-RT variables
- [ ] Test interface correctly calls POTATO routines

### 2.3 Orbit Integration via POTATO
- [ ] Write failing test for POTATO bounce in `test/test_potato_bounce.f90`
- [ ] Implement `bounce_thick()` that calls POTATO's `find_bounce()` 
- [ ] Use POTATO's existing VODE integration with guiding-center equations
- [ ] Handle POTATO's bounce time calculation and orbit closure detection
- [ ] Convert POTATO output to NEO-RT bounce averaging format
- [ ] Test thick orbit bounce times vs reference values

### 2.4 Frequency Calculation via POTATO  
- [ ] Write failing test for POTATO frequencies in `test/test_potato_frequencies.f90`
- [ ] Interface to POTATO's canonical frequency calculations
- [ ] Use POTATO's `taub` and `delphi` to compute `omega_b` and `omega_phi`
- [ ] Implement frequency caching/spline interface for performance
- [ ] Test frequency accuracy vs POTATO standalone runs

## Phase 3: Automated Testing Framework (Validation)

### 3.1 Frequency Comparison Tests
- [ ] Write comprehensive test in `test/test_frequency_comparison.f90`
- [ ] Compare bounce frequencies: `omega_b_thin` vs `omega_b_thick`
- [ ] Compare precession frequencies: `omega_prec_thin` vs `omega_prec_thick`
- [ ] Test across full pitch angle range: trapped and passing particles
- [ ] Validate convergence: thick → thin as gyroradius → 0
- [ ] Create automated tolerance checking with relative error thresholds

### 3.2 Physics Validation Suite
- [ ] Write test for orbit closure in `test/test_orbit_closure.f90`
- [ ] Validate orbit periodicity: `x(t + T_bounce) = x(t)`
- [ ] Test energy conservation along thick orbits
- [ ] Verify canonical momentum conservation
- [ ] Check first adiabatic invariant conservation
- [ ] Test Poincaré section analysis for orbit classification

### 3.3 Performance Benchmarking
- [ ] Write performance test in `test/test_performance_comparison.f90`
- [ ] Benchmark computation time: thin vs thick orbits
- [ ] Memory usage analysis for different orbit integration modes
- [ ] Accuracy vs performance trade-off analysis
- [ ] Create automated performance regression testing

## Phase 4: Bounce Integrals and Canonical Frequencies (Advanced Physics)

### 4.1 Bounce Integral Calculations
- [ ] Write failing test for bounce integrals in `test/test_bounce_integrals.f90`
- [ ] Port POTATO's bounce integral methodology to `src/bounce_integrals.f90`
- [ ] Implement general bounce integral computation: `∮ f(x) dl`
- [ ] Add support for arbitrary functionals along orbit trajectory
- [ ] Create resonant integral handling with adaptive sampling
- [ ] Test integral accuracy for known analytical cases

### 4.2 Canonical Frequency Enhancement  
- [ ] Write failing test for canonical frequencies in `test/test_canonical_frequencies.f90`
- [ ] Implement exact canonical frequency calculation from orbit data
- [ ] Add `omega_theta = 2π * dJ_theta/dT_bounce` calculation
- [ ] Add `omega_phi = 2π * dJ_phi/dT_bounce` calculation
- [ ] Handle resonant frequency calculations with root-finding
- [ ] Test frequency accuracy near resonances

### 4.3 Resonance Analysis Extension
- [ ] Write failing test for thick orbit resonances in `test/test_resonances_thick.f90`
- [ ] Extend `resonance.f90` to handle thick orbit resonances
- [ ] Implement resonance identification using thick orbit frequencies
- [ ] Add resonance width calculation from thick orbit dynamics
- [ ] Support multiple simultaneous resonances per orbit
- [ ] Test resonance location accuracy vs thin orbit approximation

## Phase 5: Configuration and User Interface (Integration)

### 5.1 Configuration Interface
- [ ] Write failing test for configuration in `test/test_thick_orbit_config.f90`
- [ ] Extend namelist in `examples/base/driftorbit.in` with thick orbit parameters
- [ ] Add `orbit_type` parameter: "thin", "thick", "auto"
- [ ] Add gyroradius calculation method selection
- [ ] Add thick orbit integration tolerances configuration
- [ ] Test configuration parsing and validation

### 5.2 Runtime Mode Switching
- [ ] Write failing test for mode switching in `test/test_orbit_mode_switch.f90`
- [ ] Implement runtime orbit mode selection in `src/neort.f90`
- [ ] Add automatic mode selection based on physics parameters
- [ ] Create seamless interface switching without code duplication
- [ ] Support hybrid calculations: thin for most, thick for resonances
- [ ] Test mode switching consistency and performance

### 5.3 Documentation and Examples
- [ ] Write failing test for example runs in `test/test_thick_orbit_examples.f90`
- [ ] Create thick orbit example in `examples/thick_orbit/`
- [ ] Document thick orbit usage in CLAUDE.md
- [ ] Add thick orbit section to existing documentation
- [ ] Create comparison plotting scripts in `python/`
- [ ] Test example runs execute successfully

## Testing Strategy and Validation Approach

### Automated Test Requirements
1. **Unit Tests**: Each module tested independently with TDD
2. **Integration Tests**: Thick and thin orbit comparison at system level
3. **Physics Tests**: Conservation laws and known analytical solutions
4. **Performance Tests**: Computational efficiency and memory usage
5. **Regression Tests**: Ensure no degradation of existing thin orbit functionality

### Frequency Comparison Validation
```fortran
! Target test structure for automated comparison
subroutine test_bounce_frequency_comparison()
    real(8) :: omega_b_thin, omega_b_thick, relative_error
    real(8), parameter :: TOLERANCE = 1.0d-3  ! 0.1% agreement required
    
    ! Calculate frequencies with both methods
    call calculate_bounce_frequency_thin(eta, omega_b_thin)
    call calculate_bounce_frequency_thick(eta, omega_b_thick)
    
    ! Validate agreement in appropriate regimes
    relative_error = abs(omega_b_thick - omega_b_thin) / omega_b_thin
    call assert_less_than(relative_error, TOLERANCE)
end subroutine
```

### Success Criteria
- [ ] Bounce frequency agreement within 0.1% for small gyroradius limit
- [ ] Precession frequency agreement within 1% across parameter space  
- [ ] Orbit closure to machine precision for both methods
- [ ] Performance degradation < 10x for thick orbit calculations
- [ ] All existing thin orbit tests continue to pass
- [ ] Comprehensive test coverage (>90%) for new thick orbit modules

## Expected Outcomes

### Physics Improvements
1. **Accurate thick orbit physics** via direct POTATO integration
2. **Enhanced resonance calculations** using POTATO's exact guiding-center dynamics  
3. **Unified orbit framework** supporting thin (NEO-RT) and thick (POTATO) calculations
4. **Validated frequency calculations** with automated testing against POTATO reference

### Code Architecture Benefits
1. **Modular design** with interface layer allowing seamless orbit method selection
2. **Backward compatibility** with all existing thin orbit functionality
3. **Performance options** from fast thin to accurate thick orbits via runtime selection
4. **Direct POTATO integration** avoiding code duplication and maintenance burden

### Revised Integration Timeline
- **Phase 1**: ✅ **COMPLETE** - Interface compatibility layer (3 phases)
- **Phase 2**: 2-3 weeks (POTATO source integration and interface layer)
- **Phase 3**: 2-3 weeks (automated testing framework) 
- **Phase 4**: 2-3 weeks (advanced physics validation)
- **Phase 5**: 1-2 weeks (configuration and documentation)
- **Total**: 7-11 weeks for complete POTATO integration

**Key Insight**: Direct POTATO integration is more efficient than porting, leveraging existing tested POTATO functionality while providing clean interface for NEO-RT integration.