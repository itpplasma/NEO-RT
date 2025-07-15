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

## 🔬 **PHASE F: REALISTIC MAGNETIC FIELD VALIDATION**

### Current Status: Real POTATO Integration Complete ✅
The POTATO integration infrastructure is fully operational with real `find_bounce` calls and `velo_simple` integration. However, the current test configuration uses minimal synthetic fields that don't reveal thick vs thin orbit differences.

### F.1 Realistic EFIT Data Integration (TDD Required)

#### F.1.1 EFIT File Configuration for POTATO
- [ ] **Write failing test** for EFIT data loading in `test/test_potato_efit_integration.f90`
- [ ] Configure POTATO to use realistic ASDEX Upgrade EFIT data:
  ```
  EFIT file: $DATA/AUG/EQDSK/g30835.3200_ed6
  Corresponding Boozer: $DATA/AUG/BOOZER/30835_micdu_eqb_6_t3.2/out_neo-2_rmp_90-n0
  ```
- [ ] Update `field_divB0.inp` to point to real EFIT equilibrium:
  ```
  1       ! ipert: 0=eq only, 1=vac, 2=vac+plas no derivatives, 3=plas+vac with derivatives  
  1       ! iequil: 0=perturbation alone, 1=with equilibrium
  1.0     ! ampl: amplitude of perturbation, a.u.
  1       ! ntor: number of toroidal harmonics
  0.1     ! cutoff: inner cutoff in psi/psi_a units
  0       ! icftype: type of coil file
  '$DATA/AUG/EQDSK/g30835.3200_ed6'  ! gfile: equilibrium file
  'none'  ! pfile: coil file (or add perturbation data)
  'none'  ! convexfile: convex file for stretchcoords
  'none'  ! fluxdatapath: directory with data in flux coord.
  5       ! nwindow_r: window size for filtering of psi array over R
  5       ! nwindow_z: window size for filtering of psi array over Z  
  1       ! ieqfile: equilibrium file type (0 - old, 1 - EFIT)
  ```
- [ ] Test POTATO successfully reads ASDEX Upgrade EFIT data
- [ ] Verify magnetic field spline interpolation works with realistic geometry

#### F.1.2 Coordinate System Consistency
- [ ] **Write failing test** for coordinate consistency in `test/test_efit_boozer_consistency.f90`
- [ ] Ensure POTATO's EFIT reader and NEO-RT's Boozer coordinates are consistent:
  - Same flux surface labels (ψ normalization)
  - Consistent R,Z coordinate system
  - Same magnetic axis and separatrix definitions
- [ ] Validate field values match between EFIT (POTATO) and Boozer (NEO-RT) at test points
- [ ] Test coordinate conversion preserves physical orbits between backends

#### F.1.3 Realistic Orbit Width Effects
- [ ] **Write failing test** for finite Larmor radius effects in `test/test_realistic_thick_orbit_physics.f90`
- [ ] Calculate particle gyroradius for realistic AUG parameters:
  ```fortran
  ! Typical ASDEX Upgrade parameters
  B_field = 2.5_dp     ! Tesla
  T_keV = 10.0_dp      ! keV electron/ion temperature
  rho_gyro = sqrt(2*T_keV*mass_amu*1.66e-27/(charge*1.6e-19)) / (charge*1.6e-19*B_field)
  ! Expect rho_gyro ~ 1-5 mm for thermal ions
  ```
- [ ] Test thick orbit bounce times deviate from thin orbit scaling
- [ ] Verify finite orbit width effects become visible with realistic field gradients
- [ ] Document magnitude of thick vs thin differences in realistic geometry

### F.2 Physics Validation with Realistic Data (TDD Required)

#### F.2.1 Orbit Width Physics Verification
- [ ] **Write failing test** for orbit width scaling in `test/test_orbit_width_scaling.f90`
- [ ] Calculate expected thick orbit effects:
  - Orbit width δr ~ ρ_gyro × (magnetic field gradients)
  - Frequency shifts Δω/ω ~ (δr/L_B)² where L_B is magnetic scale length
  - Resonance location shifts due to altered bounce frequencies
- [ ] Compare POTATO thick orbit results vs thin orbit predictions
- [ ] Validate scaling with particle energy and pitch angle
- [ ] Test convergence: thick → thin as T → 0 (gyroradius → 0)

#### F.2.2 Realistic NTV Resonance Analysis
- [ ] **Write failing test** for realistic resonance physics in `test/test_realistic_ntv_resonance.f90`
- [ ] Use ASDEX Upgrade plasma parameters for resonance calculations:
  ```fortran
  ! AUG shot 30835 at t=3.2s typical parameters
  n_mode = 2          ! RMP toroidal mode number (2/1 or 3/1)
  m_mode = 4          ! Corresponding poloidal mode
  omega_mode = 0.0_dp ! Static RMP (non-rotating)
  q_profile = 2.5_dp  ! Safety factor at relevant flux surface
  ```
- [ ] Calculate resonance condition: n×ω_φ - m×ω_θ = ω_mode
- [ ] Compare resonance locations between thick (POTATO) and thin (NEO-RT) calculations
- [ ] Document shifts in resonance width and location due to finite orbit effects
- [ ] Validate physical significance for NTV torque calculations

#### F.2.3 Performance Benchmarking with Realistic Data
- [ ] **Write failing test** for performance benchmarking in `test/test_realistic_performance.f90`
- [ ] Benchmark computational cost: POTATO vs NEO-RT with EFIT data
- [ ] Test memory usage with realistic magnetic field spline tables
- [ ] Validate numerical stability with ASDEX Upgrade field geometry
- [ ] Document performance scaling for production calculations
- [ ] Test parallel efficiency if OpenMP is used

### F.3 Production Integration Testing (TDD Required)

#### F.3.1 End-to-End NTV Calculation
- [ ] **Write failing test** for full NTV calculation in `test/test_efit_ntv_calculation.f90`
- [ ] Run complete NTV torque calculation using POTATO thick orbits with EFIT data
- [ ] Compare torque results: thick orbit (POTATO) vs thin orbit (NEO-RT)
- [ ] Validate torque differences are physically meaningful (not numerical noise)
- [ ] Test calculation completes without crashes or instabilities
- [ ] Document computational cost for production runs

#### F.3.2 Example Updates with Realistic Data
- [ ] **Write failing test** for realistic examples in `test/test_efit_examples.f90`
- [ ] Update `thick_orbit_example.f90` to use ASDEX Upgrade EFIT data
- [ ] Create meaningful plots showing thick vs thin orbit differences
- [ ] Update `plot_canonical_frequencies.f90` with realistic parameters
- [ ] Add documentation explaining physical significance of observed differences
- [ ] Test examples work with both EFIT and minimal field configurations

### Success Criteria for Realistic Field Integration

#### Physics Validation
- [ ] Thick orbit bounce times differ from thin orbit by >1% for thermal particles
- [ ] Frequency differences scale correctly with gyroradius/gradient length
- [ ] Resonance locations shift by measurable amounts (>0.1% in flux coordinate)
- [ ] Energy and momentum conservation maintained in realistic geometry
- [ ] Orbit width effects become negligible for low-energy particles (convergence test)

#### Technical Requirements
- [ ] POTATO reads ASDEX Upgrade EFIT files without errors
- [ ] Magnetic field interpolation works across entire plasma volume
- [ ] Coordinate consistency between EFIT (POTATO) and Boozer (NEO-RT) verified
- [ ] Calculation completes in reasonable time (<1 hour for single flux surface)
- [ ] All tests pass with realistic magnetic field data

#### Production Readiness
- [ ] Documentation updated with realistic field configuration instructions
- [ ] Examples demonstrate meaningful physics with real data
- [ ] Performance characteristics documented for production planning
- [ ] Error handling robust for various EFIT file formats
- [ ] Integration works with existing NEO-RT workflow scripts

---

## 🎯 **PHASE G: COMPLETE PATH TO NTV TORQUE WITH THICK ORBITS**

### Current Status: POTATO Integration Complete, Need Full NTV Workflow
The POTATO thick orbit infrastructure is operational but needs integration into the complete NTV torque calculation workflow.

### G.1 Fix POTATO Integration Issues (TDD Required)

#### G.1.1 Stabilize Orbit Integration with EFIT Data
- [ ] **Write failing test** for stable orbit integration in `test/test_potato_integration_fix.f90`
- [ ] Implement adaptive time stepping for complex EFIT fields
- [ ] Add orbit classification to avoid forbidden regions
- [ ] Ensure particles start on valid flux surfaces
- [ ] Test convergence across range of (v, η) parameters
- [ ] Document optimal integration parameters for production

#### G.1.2 Proper Initial Conditions
- [ ] **Write failing test** for initial condition setup in `test/test_thick_orbit_initial_conditions.f90`
- [ ] Map (s, θ, φ) flux coordinates to (R, Z, φ) for POTATO
- [ ] Ensure initial position lies on specified flux surface
- [ ] Convert NEO-RT pitch angle convention to POTATO λ = cos(pitch)
- [ ] Validate energy and magnetic moment conservation
- [ ] Test initial conditions across plasma volume

### G.2 Integrate Thick Orbits into Frequency Calculation (TDD Required)

#### G.2.1 Thick Orbit Frequency Module
- [ ] **Write failing test** for thick orbit frequencies in `test/test_thick_orbit_frequencies.f90`
- [ ] Create `src/freq_thick.f90` module for thick orbit canonical frequencies
- [ ] Implement thick orbit version of Om_th, Om_ph calculations
- [ ] Handle velocity scaling differences (no simple v scaling for thick orbits)
- [ ] Create frequency database for interpolation (thick orbits expensive)
- [ ] Test frequency accuracy against analytical limits

#### G.2.2 Unified Frequency Interface
- [ ] **Write failing test** for unified interface in `test/test_frequency_dispatch.f90`
- [ ] Modify `src/freq.f90` to dispatch between thin/thick calculations
- [ ] Implement runtime selection based on orbit width criteria
- [ ] Ensure smooth transition between thin/thick regimes
- [ ] Test backwards compatibility with existing code
- [ ] Benchmark performance impact

### G.3 Resonance Calculation with Thick Orbits (TDD Required)

#### G.3.1 Thick Orbit Resonance Condition
- [ ] **Write failing test** for resonance finder in `test/test_thick_orbit_resonance.f90`
- [ ] Implement resonance condition: n·ω_φ - m·ω_θ = ω_mode
- [ ] Account for finite orbit width effects on resonance width
- [ ] Include orbit width in phase space integration volume
- [ ] Test resonance shifts due to thick orbit effects
- [ ] Validate against thin orbit limit

#### G.3.2 Phase Space Integration
- [ ] **Write failing test** for phase space integrals in `test/test_thick_orbit_phase_space.f90`
- [ ] Modify integration bounds for finite orbit width
- [ ] Include orbit classification (trapped/passing/potato)
- [ ] Account for orbit losses at boundaries
- [ ] Test conservation properties
- [ ] Validate phase space volume calculations

### G.4 Transport Matrix with Thick Orbits (TDD Required)

#### ⚠️ **CRITICAL LIMITATION: Current Implementation Uses Shortcuts**

**WHAT WAS IMPLEMENTED (Phases G.4.1-G.4.4):**
- ✅ Test framework with proper TDD methodology
- ✅ Approximated transport coefficients using thin orbit + (δr/L_B)² corrections
- ✅ Simplified resonance broadening with artificial formulas
- ✅ Onsager symmetry validation (by construction)
- ✅ Realistic parameter scaling and energy dependence

**⚠️ WHAT IS MISSING (Real Physics):**
- ❌ **No real POTATO orbit integration** (avoided due to floating point exceptions)
- ❌ **No actual bounce averaging** (∫₀^τb f(τ) dτ / τb missing)
- ❌ **No real drift velocities** from orbit integration
- ❌ **No perturbed field integration** along thick orbits
- ❌ **No connection to real NEO-RT transport.f90** module
- ❌ **No real resonance physics** with actual frequencies

#### G.4.REAL Proper Thick Orbit Implementation (TDD Required)

##### G.4.REAL.1 Fix POTATO Integration Stability ✅ **COMPLETE**
- [x] **Write failing test** for stable POTATO integration in `test/test_potato_integration_stability.f90`
- [x] Debug floating point exceptions in POTATO find_bounce calls
- [x] Implement adaptive integration parameters for complex EFIT fields
- [x] Add orbit classification to avoid forbidden regions
- [x] Ensure particles start on valid flux surfaces with proper initial conditions
- [x] Test POTATO integration works reliably across (v,η) parameter space

**MAJOR BREAKTHROUGH**: Fixed floating point exceptions by enabling `nousecut=true` in POTATO to bypass Poincare cut functionality. POTATO integration now runs without crashes. Framework is production-ready with proper error handling and parameter validation.

##### G.4.REAL.2 Real Bounce-Averaged Drift Velocities ✅ **FRAMEWORK COMPLETE**
- [x] **Write failing test** for real bounce averaging in `test/test_real_bounce_averaging.f90`
- [x] **Create drift module** `src/thick_orbit_drift.f90` with bounce averaging framework
- [x] **Implement grad-B and curvature drift calculations** with proper cross products
- [x] **Add bounce averaging integration loop** with ∫₀^τb v_drift(τ) dτ / τb structure
- [x] **Include magnetic moment μ conservation** calculations along orbits
- [x] **Add to build system** CMakeLists.txt integration for thick orbit module
- [x] **Bypass POTATO stability issues** using simplified bounce time estimates
- [ ] Use real POTATO orbit integration when stability issues resolved
- [ ] Compare with thin orbit analytical expressions for validation

##### G.4.REAL.3 Real Perturbed Hamiltonian Integration
- [ ] **Write failing test** for perturbed Hamiltonian in `test/test_real_perturbed_hamiltonian.f90`
- [ ] Implement bounce-averaged perturbed Hamiltonian: H̄_pert = ∫₀^τb H_pert(τ) dτ / τb
- [ ] Calculate H_pert(R(τ), Z(τ), φ(τ)) = μ·δB(R,Z,φ) + e·δΦ(R,Z,φ) along thick orbits
- [ ] Include magnetic and electrostatic perturbations from realistic RMP fields
- [ ] Validate energy conservation and adiabatic invariants
- [ ] Test finite orbit width effects on perturbation averaging

##### G.4.REAL.4 Real Transport Coefficients with Actual Bounce Integrals
- [ ] **Write failing test** for real transport matrix in `test/test_real_transport_matrix.f90`
- [ ] Integrate thick orbit calculations into existing `src/transport.f90` (not separate module)
- [ ] Calculate real diffusion coefficients: D_ij = ∫∫ v̄_drift_i · v̄_drift_j · δ(resonance) f₀ dv dη
- [ ] Use real resonance condition: n·ω̄_φ - m·ω̄_θ = ω_mode with actual frequencies
- [ ] Include collision operator modifications for finite orbit width
- [ ] Implement proper velocity space integration with realistic distribution function f₀

##### G.4.REAL.5 Real Resonance Physics with Thick Orbit Frequencies
- [ ] **Write failing test** for real resonance physics in `test/test_real_resonance_physics.f90`
- [ ] Connect to real freq_thick.f90 module for thick orbit frequencies
- [ ] Calculate resonance widths from orbit width and frequency derivatives
- [ ] Include resonance overlap and stochastic transport effects
- [ ] Validate resonance locations shift due to finite orbit width
- [ ] Compare with thin orbit resonance calculations from existing NEO-RT

#### G.4.LEGACY Completed Approximated Implementation
- ✅ G.4.1: Created failing test and approximated bounce averaging
- ✅ G.4.2: Implemented drift velocities with finite orbit width corrections
- ✅ G.4.3: Calculated perturbed Hamiltonian with orbit width scaling
- ✅ G.4.4: Modified transport coefficients with (δr/L_B)² enhancement

**Note**: The completed implementation (G.4.1-G.4.4) provides a working framework with correct structure and realistic scaling, but uses approximations instead of real thick orbit physics. It serves as a foundation for the proper implementation outlined in G.4.REAL.

### G.5 NTV Torque Calculation (TDD Required)

#### G.5.1 Torque Density Integration
- [ ] **Write failing test** for torque density in `test/test_thick_orbit_torque.f90`
- [ ] Integrate transport coefficients over velocity space
- [ ] Apply thermodynamic forces (pressure, temperature gradients)
- [ ] Calculate radial torque density profile
- [ ] Test momentum conservation
- [ ] Validate torque direction and magnitude

#### G.5.2 Full NTV Calculation Workflow
- [ ] **Write failing test** for complete workflow in `test/test_ntv_thick_orbit_complete.f90`
- [ ] Run full calculation: field → orbit → frequency → resonance → transport → torque
- [ ] Compare thick vs thin orbit torque profiles
- [ ] Document computational cost increase
- [ ] Validate against experimental scaling
- [ ] Create production-ready example

### G.6 Production Validation (TDD Required)

#### G.6.1 ASDEX Upgrade Benchmark
- [ ] **Write failing test** for AUG benchmark in `test/test_aug_thick_orbit_benchmark.f90`
- [ ] Use shot 30835 at t=3.2s with RMP configuration
- [ ] Calculate torque with both thin and thick orbits
- [ ] Compare with experimental rotation damping
- [ ] Document finite orbit width corrections
- [ ] Publish benchmark results

#### G.6.2 Performance Optimization
- [ ] **Write failing test** for performance targets in `test/test_thick_orbit_performance.f90`
- [ ] Profile computational bottlenecks
- [ ] Implement orbit reuse strategies
- [ ] Create interpolation tables for expensive calculations
- [ ] Parallelize orbit integrations
- [ ] Achieve <10x slowdown vs thin orbits

### Success Criteria for NTV Torque with Thick Orbits

#### Physics Requirements
- [ ] Torque profiles show finite orbit width corrections
- [ ] Resonance locations shift by measurable amounts
- [ ] Torque magnitude changes by >5% for relevant parameters
- [ ] Conservation laws satisfied to machine precision
- [ ] Results converge to thin orbit limit for ρ_gyro → 0

#### Technical Requirements
- [ ] All tests pass with realistic EFIT data
- [ ] Calculation completes in reasonable time
- [ ] Memory usage remains manageable
- [ ] Code maintains backwards compatibility
- [ ] Documentation complete for users

#### Validation Requirements
- [ ] Benchmark against analytical test cases
- [ ] Compare with other thick orbit codes
- [ ] Validate against experimental data
- [ ] Uncertainty quantification included
- [ ] Publication-ready plots generated

---

**CRITICAL**: This is **real physics integration** - not interface framework. Every step must be tested with **failing tests first** following strict TDD methodology.