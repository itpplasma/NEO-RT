# TODO: Real POTATO Integration for Thick Orbit Calculations

## Project Goal
Integrate actual POTATO functionality into NEO-RT to enable real thick orbit calculations, replacing the current stub implementation with genuine finite Larmor radius physics.

## Current Status: Interface Framework Complete ✅

### ✅ **DONE: Interface Architecture (Phases 1-5)**
- **Phase 1**: Interface compatibility layer
- **Phase 2**: POTATO integration architecture with runtime dispatch 
- **Phase 3**: Comprehensive automated testing framework
- **Phase 4**: Advanced physics validation and benchmarking
- **Phase 5**: Configuration, documentation, and examples

### ⚠️ **REALITY CHECK: Stub Implementation Only**
**Current thick_orbit_example.x uses:**
- **Artificial physics**: Stub generates fake eta-dependent values
- **Placeholder differences**: Example multiplies stub results by constants  
- **No real POTATO calls**: All "thick orbit" calculations are synthetic

**What we have**: Production-ready interface framework  
**What we need**: Real POTATO function calls for actual thick orbit physics

## 🔧 **REAL POTATO INTEGRATION PLAN**

**📋 POTATO Structure Reference**: See `POTATO/README.md` for complete codebase organization, key functions, and integration points.

### Phase A: CMake Build Integration (TDD Required)

#### A.1 POTATO Library Integration
- [ ] **Write failing test** for POTATO library availability in `test/test_potato_build.f90`
- [ ] Add POTATO subdirectory to main CMakeLists.txt when `USE_THICK_ORBITS=ON`
- [ ] Create `find_package(POTATO)` or direct subdirectory approach
- [ ] Handle POTATO's dependency on BLAS/LAPACK (already in NEO-RT)
- [ ] Handle POTATO's dependency on VODE (conflict resolution with NEO-RT's VODE)
- [ ] Test POTATO library (`potato_base`) builds successfully within NEO-RT
- [ ] Test POTATO modules are available for linking

#### A.2 Dependency Resolution
- [ ] **Write failing test** for VODE conflict resolution in `test/test_vode_conflict.f90`
- [ ] Resolve POTATO's VODE vs NEO-RT's VODE integration
- [ ] Ensure POTATO's `odeint_allroutines.f` doesn't conflict with NEO-RT
- [ ] Create proper linking order: `potato_base` → `vode` → NEO-RT modules
- [ ] Test no symbol conflicts between POTATO and NEO-RT

### Phase B: Magnetic Field Interface Bridge (TDD Required)

#### B.1 Field Interface Analysis
**CRITICAL DISCOVERY**: POTATO's `psif`, `dpsidr`, `dpsidz` are **module variables**, not functions!

- [ ] **Write failing test** for field bridge in `test/test_potato_field_bridge.f90` 
- [ ] Create `src/potato_field_bridge.f90` module
- [ ] Interface NEO-RT's `magfie` calls to POTATO's field variables:
  ```fortran
  ! POTATO: module variables set by spline interpolation
  use field_sub, only: psif, dpsidr, dpsidz
  
  ! NEO-RT: function calls
  call psif_neo_rt(R, Z, psi_result)
  call dpsidr_neo_rt(R, Z, dpsi_dr_result) 
  call dpsidz_neo_rt(R, Z, dpsi_dz_result)
  ```
- [ ] Call POTATO's `field_eq(R, Z)` to set these variables before each access
- [ ] Handle coordinate system conversion if needed
- [ ] Test field values match between NEO-RT and POTATO approaches

#### B.2 POTATO Field Setup Integration
- [ ] **Write failing test** for POTATO field initialization in `test/test_potato_field_init.f90`
- [ ] Initialize POTATO's magnetic field data from NEO-RT's input files
- [ ] Read POTATO's field data (equilibrium file compatibility)
- [ ] Set up POTATO's spline interpolation tables
- [ ] Call POTATO's field initialization routines before orbit calculations
- [ ] Test POTATO field setup doesn't interfere with NEO-RT's magfie

### Phase C: Core POTATO Function Integration (TDD Required)

#### C.1 Replace Stub with Real `find_bounce`
- [ ] **Write failing test** for real POTATO bounce in `test/test_potato_find_bounce_real.f90`
- [ ] Replace `potato_stub.f90` calls with real POTATO `sub_potato` module:
  ```fortran
  use sub_potato, only: find_bounce  ! Real POTATO function
  ! Remove: use potato_stub, only: potato_find_bounce  ! Stub
  ```
- [ ] Handle POTATO's `find_bounce` interface correctly:
  ```fortran
  ! POTATO signature:
  subroutine find_bounce(next,velo_ext,dtau_in,z_eqm,taub,delphi,extraset)
  !   next - number of extra integrals
  !   velo_ext - external velocity routine 
  !   dtau_in - max time step
  !   z_eqm(5) - phase space variables [R,Z,phi,v_par,v_perp]
  !   taub - bounce time (output)
  !   delphi - toroidal shift per bounce (output)  
  !   extraset(next) - extra integrals (inout)
  ```
- [ ] Convert NEO-RT parameters (v, eta) to POTATO phase space (z_eqm)
- [ ] Handle POTATO's velocity routine requirement (`velo_ext`)
- [ ] Test real bounce times vs stub values show expected physics

#### C.2 POTATO Orbit Integration Setup
- [ ] **Write failing test** for POTATO orbit setup in `test/test_potato_orbit_setup.f90`
- [ ] Initialize POTATO's global state variables needed for orbit integration:
  ```fortran
  use orbit_dim_mod   ! POTATO's orbit dimensions and settings
  use parmot_mod      ! POTATO's motion parameters
  use global_invariants ! POTATO's conserved quantities
  ```
- [ ] Set POTATO's orbit integration parameters (tolerances, etc.)
- [ ] Handle POTATO's file I/O if orbit writing is needed
- [ ] Create proper velocity routine for POTATO's `velo_ext` parameter
- [ ] Test POTATO's orbit integration initialization is complete

#### C.3 Coordinate System Conversion
- [ ] **Write failing test** for coordinate conversion in `test/test_coordinate_conversion.f90`
- [ ] Convert NEO-RT's (v, eta) → POTATO's phase space z_eqm(5):
  ```fortran
  ! NEO-RT input: velocity v, pitch parameter eta
  ! POTATO needs: z_eqm = [R, Z, phi, v_parallel, v_perpendicular]
  ! Conversion logic:
  v_parallel = v * sqrt(1 - eta)     ! From pitch angle
  v_perpendicular = v * sqrt(eta)    ! From pitch angle  
  R = R_initial                      ! Starting position
  Z = Z_initial                      ! Starting position
  phi = 0.0d0                       ! Starting toroidal angle
  ```
- [ ] Handle particle starting position selection
- [ ] Convert POTATO's results back to NEO-RT format
- [ ] Test coordinate conversion preserves physical meaning

### Phase D: Physics Validation and Integration Testing (TDD Required)

#### D.1 Real vs Stub Comparison
- [ ] **Write failing test** for real vs stub comparison in `test/test_real_vs_stub_physics.f90`
- [ ] Run same orbit calculation with stub and real POTATO
- [ ] Verify real POTATO shows different (more accurate) physics than stub
- [ ] Document differences in bounce times and frequencies
- [ ] Validate real POTATO frequencies show velocity-dependent effects
- [ ] Test thick orbit effects are visible (non-zero finite Larmor radius effects)

#### D.2 Thick vs Thin Orbit Physics Validation  
- [ ] **Write failing test** for physics validation in `test/test_thick_thin_physics.f90`
- [ ] Compare real thick orbit (POTATO) vs thin orbit (NEO-RT) frequencies
- [ ] Verify thick orbits show finite orbit width effects
- [ ] Test convergence: thick → thin as gyroradius → 0
- [ ] Validate resonance calculations using real thick orbit frequencies
- [ ] Document performance impact of real POTATO vs stub

#### D.3 End-to-End Integration Testing
- [ ] **Write failing test** for full integration in `test/test_full_potato_integration.f90`
- [ ] Test complete orbit calculation workflow with real POTATO
- [ ] Verify NEO-RT's `driftorbit.f90` works with thick orbit option
- [ ] Test torque calculations using thick orbit bounce integrals
- [ ] Validate example programs work with real POTATO (not just stub)
- [ ] Performance benchmark: real POTATO computational cost vs thin orbits

### Phase E: Documentation and Examples Update (TDD Required)

#### E.1 Update Example Programs
- [ ] **Write failing test** for updated examples in `test/test_potato_examples.f90`
- [ ] Update `thick_orbit_example.f90` to use real POTATO (remove placeholders)
- [ ] Remove artificial difference generation (lines 78-80 in example)
- [ ] Add real physics comparison between thick and thin orbits
- [ ] Create plotting scripts showing real thick vs thin orbit differences
- [ ] Test examples demonstrate meaningful physics (not fake differences)

#### E.2 Documentation Updates
- [ ] Update CLAUDE.md with real POTATO integration status
- [ ] Document POTATO build requirements and dependencies
- [ ] Add instructions for magnetic field compatibility requirements
- [ ] Document performance expectations for thick vs thin orbits
- [ ] Create troubleshooting guide for POTATO integration issues

## Expected Implementation Challenges

### Technical Challenges
1. **VODE Conflict**: POTATO has its own VODE, NEO-RT has VODE - need proper linking
2. **Field Interface**: POTATO uses module variables, NEO-RT uses function calls
3. **Coordinate Systems**: Conversion between NEO-RT and POTATO representations
4. **Global State**: POTATO may rely on module-level state variables
5. **File I/O**: POTATO may expect specific input files or formats

### Physics Challenges
1. **Starting Conditions**: Converting (v,eta) to proper POTATO phase space
2. **Velocity Routine**: POTATO requires external velocity function (`velo_ext`)
3. **Field Setup**: Ensuring POTATO's field matches NEO-RT's magnetic field
4. **Time Scales**: Proper time normalization between POTATO and NEO-RT
5. **Conserved Quantities**: Validating energy and momentum conservation

### Performance Challenges
1. **Computational Cost**: Real POTATO will be much slower than stub
2. **Memory Usage**: POTATO's orbit integration may require significant memory
3. **Convergence**: POTATO's orbit closure may require iterative solving
4. **Scaling**: Performance with large parameter scans

## Success Criteria for Real Integration

### Functional Requirements
- [ ] Real POTATO `find_bounce` successfully called from NEO-RT
- [ ] Bounce times show physical differences from thin orbit calculations
- [ ] Frequencies demonstrate finite Larmor radius effects
- [ ] No crashes or numerical instabilities
- [ ] All existing NEO-RT functionality remains intact

### Physics Requirements  
- [ ] Thick orbit bounce times > thin orbit bounce times (typically)
- [ ] Frequencies show velocity-dependent deviations from thin orbit scaling
- [ ] Resonance locations shift due to thick orbit effects
- [ ] Energy and momentum conservation during orbit integration
- [ ] Proper convergence to thin orbit limit for small gyroradius

### Performance Requirements
- [ ] Real POTATO integration completes in reasonable time (< 10x slower than thin)
- [ ] Memory usage remains manageable for production calculations
- [ ] Numerical stability maintained across parameter ranges
- [ ] Test suite passes with real POTATO (not just stub)

## Key Implementation Insights

**Technical Discoveries**:
1. **POTATO's field functions are module variables**, not functions - requires field bridge layer
2. **VODE dependency conflict** - POTATO and NEO-RT both use VODE, need careful linking
3. **Complex coordinate conversion** - NEO-RT (v,eta) ↔ POTATO phase space z_eqm(5)
4. **Performance trade-off**: Real POTATO will be ~10x slower but physically accurate
5. **No spline optimization possible** - thick orbits require direct integration per particle
6. **Interface framework ready** - all architecture in place for systematic POTATO integration

**⚠️ IMPORTANT LIMITATION: Spline Scaling Invalid for Thick Orbits**

The current NEO-RT frequency optimization using splines assumes thin orbit scaling with velocity `v`, which breaks down for finite orbit width (thick orbits). Key impacts:

1. **Velocity Scaling Breakdown**: Current `freq.f90` splines assume frequencies scale predictably with particle velocity, but thick orbits have complex velocity-dependent drift physics
2. **Finite Orbit Width Effects**: Thick orbits sample different magnetic field regions, invalidating simple scaling relationships
3. **Performance vs Accuracy Trade-off**: Must use direct POTATO bounce integration for each (v,eta) point instead of spline interpolation
4. **No Frequency Splines for Thick Orbits**: Unlike thin orbits, thick orbit frequencies cannot be pre-computed and interpolated

**Consequences**:
- Thick orbit calculations will be computationally more expensive than thin orbits
- Need direct integration for every orbit in thick orbit mode
- Frequency calculations must call POTATO's `find_bounce()` for each evaluation
- Cannot leverage existing `freq.f90` spline optimization infrastructure for thick orbits

---

**CRITICAL**: This is **real physics integration** - not interface framework. Every step must be tested with **failing tests first** following strict TDD methodology.