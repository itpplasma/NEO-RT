# POTATO - Particle Orbit Tracing and Transport Operations

POTATO is a specialized code for calculating guiding-center orbit dynamics in tokamak magnetic fields, providing exact thick orbit calculations that include finite Larmor radius effects.

## Project Structure

### Core Source Code (`SRC/`)

#### **Orbit Integration and Physics**
- `sub_potato.f90` - **Main orbit integration routines**
  - `find_bounce()` - **Primary function**: Integrates orbit over one bounce time
  - Handles orbit closure, bounce time calculation, and toroidal shift (delphi)
  - **Interface**: `subroutine find_bounce(next,velo_ext,dtau_in,z_eqm,taub,delphi,extraset)`
- `velo.f90` - Velocity equations for guiding-center dynamics
- `period_mod.f90` - Periodic orbit handling and classification

#### **Magnetic Field Interface**
- `field_eq_mod.f90` - **Equilibrium magnetic field module**
  - Defines field variables: `psif`, `dpsidr`, `dpsidz` (module variables, not functions)
- `field_divB0.f90` - **Field evaluation routines** 
  - Contains `field_eq()` subroutine that calls spline interpolation
  - Sets field module variables via 2D spline evaluation
- `field_mod.f90` - General magnetic field utilities
- `field_c_mod.f90` - Cylindrical coordinate field operations
- `magfie_cyl.f90` - Magnetic field in cylindrical coordinates

#### **Spline Interpolation**
- `spline5_RZ.f90` - **2D spline interpolation for magnetic field**
  - `spline()` subroutine: Evaluates ψ(R,Z) and derivatives
  - **Critical for POTATO integration**: Sets field module variables
- `spl_three_to_five.f90` - Spline utilities and transformations
- `plag_coeff.f90` - Lagrange coefficient calculations for interpolation

#### **Mathematical and Numerical Tools**
- `odeint_allroutines.f` - **ODE integration routines (potential VODE conflict)**
- `math_constants.f90` - Mathematical constants and definitions
- `libneo_kinds.f90` - Precision definitions and data types
- `find_all_roots.f90` - Root finding algorithms for orbit closure

#### **Coordinate Systems and Geometry**
- `theta_rz_mod.f90` - Coordinate transformations (θ,R,Z)
- `extract_fluxcoord_mod.f90` - Flux coordinate extraction
- `bdivfree_mod.f90` - Divergence-free field operations
- `inthecore_mod.f90` - Core region identification

#### **Data Analysis and Output**
- `sample_matrix.f90` - Matrix sampling for data collection
- `sample_matrix_out.f90` - Output formatting for sample matrices
- `binsrc.f90` - Binary search utilities
- `sorting.f90` - Data sorting algorithms

#### **Equilibrium and Profiles**
- `eqmagprofs.f90` - Equilibrium magnetic profiles
- `equimoments.f90` - Equilibrium moments calculation
- `amn_mod.f90` - Amplitude mode analysis
- `bmod_pert.f90` - Magnetic perturbation handling

#### **Input/Output and Configuration**
- `input_files.f90` - File I/O handling
- `profile_input.f90` - Profile data input (main version)
- `profile_input_fixed.f90` - Fixed profile input variant

#### **Special Analysis Tools**
- `resonant_int.f90` - Resonant integral calculations
- `box_counting.f90` - Fractal dimension analysis
- `alpha_lifetime_mod.f90` - Alpha particle lifetime calculations

#### **Testing and Examples**
- `tt.f90` - Main POTATO executable
- `test_profile.f90` - Profile testing routines

### Build System

#### **CMake Configuration**
- `CMakeLists.txt` - Main CMake configuration
- `SRC/CMakeLists.txt` - Source-specific build rules
- **Libraries**: 
  - `potato_base` - Core POTATO functionality (35 source files)
  - `potato` - Main library with profile input
  - `potato_fixed` - Fixed profile variant

#### **Make Interface**
- `Makefile` - High-level Make interface to CMake
- **Targets**: `build`, `test`, `clean`, `reconfigure`
- **Configurations**: Release (default), Debug, Fast

#### **Dependencies**
- **BLAS/LAPACK** - Linear algebra operations
- **VODE** - ODE integration (potential conflict with NEO-RT's VODE)
- **Fortran 2008** - Modern Fortran features required

### Executables

#### **Main Programs**
- `potato.x` - Primary POTATO executable (from `tt.f90`)
- `test_profile.x` - Profile testing executable
- `test_profile_fixed.x` - Fixed profile testing variant

### Auxiliary Tools

#### **BOOZER_TO_EFIT/** - Coordinate Conversion
- Utilities for converting between Boozer coordinates and EFIT format
- **Purpose**: Magnetic equilibrium preprocessing
- Contains plotting tools (`PLOT/plot.jl`) and source code

#### **TEX/** - Documentation and Results
- LaTeX documentation (`equilmaxw.tex`)
- **FIGURES/**: Extensive collection of analysis plots
  - Bounce time plots: `taub*.png`
  - Toroidal shift plots: `delphi*.png`  
  - Resonance analysis: `res_*.png`
  - Orbit visualization: `o*.png`

#### **RUN/** - Execution Environment
- Runtime configuration files
- Input examples and batch scripts

## Key Integration Points for NEO-RT

### **Critical Functions for NEO-RT Integration**

1. **`find_bounce()` in `sub_potato.f90`**
   - **Purpose**: Primary orbit integration function
   - **Returns**: `taub` (bounce time), `delphi` (toroidal shift per bounce)
   - **Requirements**: Proper velocity routine (`velo_ext`), phase space setup

2. **Field Interface in `field_eq_mod.f90` + `field_divB0.f90`**
   - **Challenge**: POTATO uses module variables, NEO-RT uses function calls
   - **Bridge needed**: Convert NEO-RT field calls to POTATO field variable access
   - **Setup required**: Initialize POTATO's spline interpolation from NEO-RT data

3. **Coordinate Conversion**
   - **Input**: NEO-RT (v, eta) → POTATO phase space z_eqm(5) = [R,Z,φ,v_∥,v_⊥]
   - **Output**: POTATO results → NEO-RT bounce averaging format

### **Integration Challenges**

1. **VODE Conflict**: Both POTATO and NEO-RT use VODE for ODE integration
2. **Module Variables vs Functions**: Different field interface paradigms
3. **Global State**: POTATO may require module-level initialization
4. **Performance**: Real thick orbit calculations will be ~10x slower than thin orbits

### **Build Integration Requirements**

1. **Add POTATO as subdirectory** when `USE_THICK_ORBITS=ON`
2. **Link against**: `potato_base` library
3. **Handle dependencies**: BLAS/LAPACK (already in NEO-RT), VODE (conflict resolution)
4. **Module availability**: Ensure POTATO modules accessible from NEO-RT source

## Physics Capabilities

### **Thick Orbit Effects**
- **Finite Larmor radius**: Exact guiding-center dynamics include orbit width effects
- **Drift physics**: Full drift motion in 3D magnetic geometry
- **Resonance analysis**: Accurate resonance identification using exact frequencies

### **Computational Features**
- **Adaptive integration**: Variable time stepping for orbit closure
- **Bounce time calculation**: Exact integration over complete bounce period
- **Canonical frequencies**: Toroidal shift (delphi) provides ω_φ for resonance analysis
- **Conservation laws**: Energy and momentum conservation during orbit integration

## Development Status for NEO-RT Integration

### **✅ Ready for Integration**
- Complete POTATO codebase with tested orbit integration
- Established build system (CMake + Make)
- Clear interface functions (`find_bounce()`)
- Documented field interface structure

### **🔧 Integration Work Required**
- CMake integration with NEO-RT build system
- Magnetic field interface bridge (module variables ↔ function calls)
- Coordinate conversion routines (NEO-RT ↔ POTATO formats)
- VODE dependency conflict resolution
- Performance optimization for production calculations

**Purpose**: POTATO provides the real physics foundation for thick orbit calculations in NEO-RT, replacing the current stub implementation with genuine finite Larmor radius effects.