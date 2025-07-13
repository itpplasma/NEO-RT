# TODO: Thick Orbit Integration for NEO-RT

## Project Goal
Extend NEO-RT main code to support full guiding-center (thick) orbit calculations like in POTATO, with automated testing comparing bounce and precession frequencies between thin and thick orbit implementations.

## Phase 1: Interface Compatibility Layer (TDD Foundation)

### 1.1 Create Thick Orbit Abstract Interface
- [ ] Write failing test for thick orbit interface in `test/test_thick_orbit_interface.f90`
- [ ] Create abstract orbit type `orbit_type_t` in new `src/orbit_types.f90` module
- [ ] Define common interface for bounce calculations: `calculate_bounce_time()`, `calculate_frequencies()`
- [ ] Implement thin orbit wrapper using existing `bounce()` routine
- [ ] Write basic integration test comparing thin orbit wrapper to direct calls

### 1.2 Magnetic Field Interface Extension  
- [ ] Write failing test for thick orbit magnetic field interface in `test/test_magfie_thick.f90`
- [ ] Create `magfie_thick_mod` module in `src/magfie_thick.f90`
- [ ] Implement finite gyroradius magnetic field evaluation
- [ ] Add gyroradius calculation routines from POTATO's approach
- [ ] Test magnetic field continuity as gyroradius → 0

### 1.3 Build System Integration
- [ ] Write failing test for CMake thick orbit option in `test/test_build_thick_orbits.cmake`
- [ ] Add `USE_THICK_ORBITS` CMake option to `CMakeLists.txt`
- [ ] Create conditional compilation flags for thick orbit modules
- [ ] Integrate POTATO source files into build when thick orbits enabled
- [ ] Test build works with both thick orbit modes (ON/OFF)

## Phase 2: POTATO Integration (Core Functionality)

### 2.1 Orbit Integration Enhancement
- [ ] Write failing test for POTATO-style bounce integration in `test/test_potato_bounce.f90`
- [ ] Port POTATO's `find_bounce()` routine to `src/orbit_thick.f90`
- [ ] Implement `bounce_thick()` function with same interface as `bounce()`
- [ ] Add support for 5D guiding-center equations: (R, φ, Z, p, λ)
- [ ] Implement adaptive ODE integration with exact bounce detection
- [ ] Test bounce time accuracy for simple magnetic configurations

### 2.2 Velocity Field Implementation
- [ ] Write failing test for guiding-center velocity in `test/test_gc_velocity.f90`
- [ ] Port POTATO's `velo()` routines to `src/velocity_gc.f90`
- [ ] Implement complete magnetic drift calculations (grad-B, curvature, E×B)
- [ ] Add relativistic corrections through gamma factor
- [ ] Handle both trapped and passing particle orbits uniformly
- [ ] Test velocity field conservation properties

### 2.3 Frequency Calculation Extension
- [ ] Write failing test for thick orbit frequencies in `test/test_frequencies_thick.f90`
- [ ] Create `freq_thick.f90` module extending `freq.f90` interface
- [ ] Implement `Om_th_thick()` and `Om_ph_thick()` functions
- [ ] Add direct frequency calculation from orbit integration
- [ ] Maintain spline interface compatibility for performance
- [ ] Test frequency calculation accuracy vs POTATO reference

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
1. **Accurate thick orbit physics** for large gyroradius regimes
2. **Enhanced resonance calculations** with exact orbit dynamics  
3. **Unified orbit framework** supporting thin and thick calculations
4. **Validated frequency calculations** with automated testing

### Code Architecture Benefits
1. **Modular design** allowing easy orbit method selection
2. **Backward compatibility** with all existing thin orbit functionality
3. **Performance options** from fast thin to accurate thick orbits
4. **Extensible framework** for future orbit enhancement

### Integration Timeline
- **Phase 1-2**: 4-6 weeks (core functionality)
- **Phase 3**: 2-3 weeks (testing framework) 
- **Phase 4**: 3-4 weeks (advanced physics)
- **Phase 5**: 2 weeks (integration and documentation)
- **Total**: 11-15 weeks for complete thick orbit integration

This plan ensures robust, tested integration of POTATO-style thick orbit functionality while maintaining the reliability and performance of existing NEO-RT thin orbit calculations.