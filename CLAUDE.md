# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System and Commands

NEO-RT uses a hybrid build system with Make and CMake:

### Dependencies
Install system dependencies:
```bash
make deps          # Auto-detects OS and installs deps (Debian/Ubuntu)
make deps-debian   # For Debian/Ubuntu systems
```

### Building
```bash
make               # Build everything (default: Release mode)
make CONFIG=Debug  # Build in Debug mode
make CONFIG=Fast   # Build in Fast mode with aggressive optimizations

# Enable thick orbit support (POTATO integration)
make CONFIG=Debug USE_THICK_ORBITS=ON
```

### Testing
```bash
make test          # Run all tests (both ctest and pytest)
make ctest         # Run C++/Fortran tests via CTest
make pytest        # Run Python tests
```

### Other Commands
```bash
make clean         # Remove build directory
make reconfigure   # Force reconfigure CMake
```

## Project Architecture

NEO-RT calculates neoclassical toroidal viscosity (NTV) torque in resonant transport regimes using a Hamiltonian approach. 

### **THICK ORBIT BRANCH GOAL**
**This branch (potato) implements the same NTV torque functionality as standard NEO-RT thin orbit approximation, but using thick guiding-center orbits from POTATO.**

This means: frequencies, resonances, transport coefficients, and NTV torque calculations must work cleanly with thick orbits. No shortcuts, no approximations - full physics implementation.

### Code Structure

### Core Components

1. **Main Entry Point**: `src/main.f90` → `src/neort.f90`
   - Initializes magnetic field, profiles, and orchestrates the calculation

2. **Core Physics Modules**:
   - `src/driftorbit.f90`: Drift orbit calculations and resonant transport
   - `src/magfie.f90`: Magnetic field interface and calculations
   - `src/profiles.f90`: Plasma profiles and thermodynamic forces
   - `src/freq.f90`: Frequency calculations (bounce, transit, etc.)
   - `src/orbit.f90`: Particle orbit integration
   - `src/resonance.f90`: Resonance identification and analysis
   - `src/transport.f90`: Transport coefficient calculations
   - `src/nonlin.f90`: Nonlinear physics calculations

3. **Utility Modules**:
   - `src/util.f90`: General utilities and mathematical functions
   - `src/attenuation_factor.f90`: Attenuation factor calculations
   - `src/collis_nbi.f90`: Collision and NBI-related functions

### Key Dependencies
- Uses external libraries: spline, vode, BLAS/LAPACK, SuiteSparse, NetCDF
- Optionally integrates with NEO-2 (controlled by USE_STANDALONE cmake option)
- **POTATO Integration**: Thick orbit support via `USE_THICK_ORBITS=ON` (includes POTATO subdirectory)

### Input Files Structure
NEO-RT requires specific input files in the working directory:
- `<runname>.in`: Main namelist configuration (see `examples/base/driftorbit.in`)
- `in_file`: Boozer coordinate file for axisymmetric magnetic field
- `in_file_pert`: Boozer coordinate file for magnetic perturbations
- `plasma.in`: Plasma thermodynamic profiles (required for torque calculations)
- `profile.in`: Rotation profiles (required for nonlinear calculations)

### Output Structure
- `<runname>_torque.out`: Torque density data
- `<runname>_magfie_param.out`: Magnetic field parameters
- Various other diagnostic outputs depending on configuration

## Development Workflow

### Running Single Calculations
```bash
# Build first
make

# Run with input file (without .in extension)
./build/neo_rt.x runname
```

### Batch Processing
Use Python scripts for flux surface scans:
```bash
# Batch runs across flux surfaces
python3 python/run_driftorbit.py

# Collect results from multiple runs
python3 python/collect_data_from_individual_runs.py
```

### Testing Strategy
- Fortran unit tests in `test/` directory using CMake/CTest
- Python integration tests in `test/ripple_plateau/` 
- Example configurations in `examples/` directory
- POTATO integration tests: `test/test_potato_*.f90` (8 modules for build, field bridge, physics validation)

#### Running Specific Tests
```bash
# Run single test
ctest -R test_potato_build

# Run all POTATO tests  
ctest -R test_potato

# Run specific example
./build/thick_orbit_example.x
./build/plot_canonical_frequencies.x
```

### Code Organization
- Main source in `src/`
- Tests in `test/`
- Examples and plotting utilities in `examples/`
- Python utilities in `python/`
- POTATO sub-project for magnetic field preprocessing
- Documentation in `doc/`

The codebase follows Fortran 90+ standards with extensive use of modules for physics calculations and utilities.

## Coding Standards

### SOLID Principles (Adapted for Fortran) - MANDATORY

**⚠️ CRITICAL: ALL SOLID PRINCIPLES MUST BE FOLLOWED ⚠️**

**S - Single Responsibility**: Each routine has one clear purpose, max 30 lines - **ENFORCED**
**O - Open/Closed**: Extend through inheritance/composition, not modification - **REQUIRED**
**L - Liskov Substitution**: Derived types must work wherever base types do - **MANDATORY**
**I - Interface Segregation**: Keep interfaces focused and minimal - **ENFORCED**
**D - Dependency Inversion**: Depend on abstractions (abstract types), not concrete implementations - **REQUIRED**

### DRY and KISS Principles - MANDATORY

**⚠️ CRITICAL: DRY AND KISS ARE STRICTLY ENFORCED ⚠️**

**DRY - Don't Repeat Yourself**: Extract common functionality into shared modules - **REQUIRED**
- Create common modules for shared logic
- Use procedure pointers for generic operations  
- Centralize constants and magic numbers in one place

**KISS - Keep It Simple, Stupid**: Favor simplicity over cleverness - **MANDATORY**
- Write clear, readable code over "clever" optimizations
- Use straightforward algorithms unless performance demands complexity
- Prefer explicit over implicit behavior
- Choose clear variable names over short abbreviations

### Test-Driven Development (MANDATORY)

**⚠️ CRITICAL: TDD IS MANDATORY FOR ALL FEATURES AND REFACTORING ⚠️**
**⚠️ WRITE TESTS FIRST - NO EXCEPTIONS ⚠️**
**⚠️ DO NOT WRITE ANY CODE WITHOUT A FAILING TEST FIRST ⚠️**

**MANDATORY TDD WORKFLOW - NEVER DEVIATE:**

1. **WRITE FAILING TEST FIRST** in `test/test_*.f90` - **ALWAYS START HERE**
2. **RUN `make test`** to confirm the test fails (RED)
3. **Write minimal code** to make test pass (GREEN)
4. **Refactor** while keeping tests green (REFACTOR)
5. **Repeat RED-GREEN-REFACTOR** for next feature

**FORBIDDEN:**
- Writing implementation code before tests
- Changing code without a test covering the change
- Assuming existing code works without tests
- Skipping tests "just this once"

**TDD is not optional** - it is the foundation of all development in this codebase.
**TESTS FIRST, ALWAYS. NO CODE WITHOUT TESTS.**

### Code Organization (MANDATORY RULES)

**⚠️ CRITICAL: THESE RULES ARE NON-NEGOTIABLE ⚠️**

**Routine Size**: Max 30 lines, single responsibility - **NO EXCEPTIONS**
**Naming**: Use descriptive verbs (`calculate_bounds` not `calc`) - **REQUIRED**
**Placement**: Helper routines after caller, shared utilities at module end - **ENFORCED**
**Comments**: Only for complex algorithms, let code self-document - **MANDATORY**

### State Management - CRITICALLY IMPORTANT

**⚠️ MUTABLE GLOBAL STATE IS THE SOURCE OF ALL EVIL AND MUST BE AVOIDED AT ALL COST ⚠️**

**MANDATORY PRINCIPLES:**
- **NO GLOBAL MUTABLE STATE** - All state must be explicitly passed as parameters
- **IMMUTABLE BY DEFAULT** - Prefer immutable data structures and pure functions
- **EXPLICIT STATE MANAGEMENT** - Always save/restore state when temporarily modifying context
- **STATELESS OPERATIONS** - Functions should not rely on hidden global state
- **CLEAR OWNERSHIP** - Each piece of state must have a clear owner and scope

```fortran
! GOOD: Explicit state save/restore
subroutine good_draw_text(ctx, text)
    real(8) :: saved_width
    saved_width = ctx%current_line_width  ! Save state
    call ctx%set_line_width(0.5d0)       ! Modify
    ! ... draw text ...
    call ctx%set_line_width(saved_width) ! Restore - MANDATORY
end subroutine

! GOOD: Pure functions with explicit parameters
pure function calculate_position(x, y, offset) result(new_pos)
    ! No hidden state dependencies
end function
```

### Constants and Magic Numbers - STRICTLY ENFORCED

**⚠️ MAGIC NUMBER CONSTANTS ARE FORBIDDEN ⚠️**

**MANDATORY PRINCIPLES:**
- **NO MAGIC NUMBERS** - If a number has meaning, it MUST be a named constant
- **DESCRIPTIVE NAMES** - Constant names must clearly indicate their purpose
- **CENTRALIZED CONSTANTS** - Group related constants in parameter declarations
- **DOCUMENTED PURPOSE** - Each constant should have a clear comment explaining its meaning

```fortran
! GOOD: Named constants with clear meaning
real(8), parameter :: DEFAULT_TOLERANCE = 1.0d-12    ! Convergence tolerance
integer, parameter :: MAX_ITERATIONS = 1000          ! Maximum solver iterations
real(8), parameter :: SAFETY_FACTOR = 0.8d0         ! CFL safety factor

! Usage
if (residual < DEFAULT_TOLERANCE) then
    converged = .true.
endif
```

### Fortran-Specific Standards - STRICTLY ENFORCED

**⚠️ CRITICAL: THESE RULES HAVE NO EXCEPTIONS ⚠️**

- Always explicitly import with `use only`. No wildcard imports allowed. - **MANDATORY**
- Use `implicit none` in all modules and programs - **REQUIRED**
- **ALL variable declarations MUST come before any executable code in routines** - **MANDATORY**
  - Variables, parameters, and type declarations first
  - Then executable statements and assignments
  - Fortran requires this strict ordering
- Use double precision (`real(8)`) for all floating-point calculations - **REQUIRED**
- **cd COMMAND IS FORBIDDEN** - Never use `cd` in bash commands. Use absolute paths instead. - **MANDATORY**

## Thick Orbit Integration (POTATO)

### ⚠️ POTATO Integration Status - **STUB IMPLEMENTATION ONLY**
**After comprehensive code audit, the thick orbit integration is NOT complete:**

- **Real POTATO Integration**: **STUB ONLY** - All `find_bounce` calls are commented out or use stubs
- **Field Bridge Layer**: Framework exists but real integration is disabled (`src/potato_field_bridge.f90` line 361)
- **Spline Interpolation**: Framework exists but not connected to real POTATO physics
- **Runtime Architecture**: Infrastructure exists but falls back to thin orbit approximations
- **Test Suite**: Tests exist but expect stub behavior, not real POTATO physics
- **Examples**: All examples use stub implementations, not real thick orbit calculations

### POTATO Build Configuration
When `USE_THICK_ORBITS=ON` is enabled in CMake:
- POTATO subdirectory is included and built as `potato_base` library
- Field evaluation functions (`psif`, `dpsidr`, `dpsidz`) added to `field_eq_mod.f90`
- VODE conflict resolved by disabling `box_counting.f90` in POTATO build
- `field_divB0.inp` configuration file required for POTATO field initialization

### Key POTATO Integration Components - **CURRENT STATUS**
- **Real Integration**: **STUB ONLY** - `real_find_bounce_calculation()` calls stub, not actual POTATO
- **Field Interface**: Framework exists but not connected to real POTATO physics
- **Coordinate Conversion**: Basic conversion exists but uses hardcoded values
- **Physics Results**: **NO REAL PHYSICS** - All calculations use simplified estimates or thin orbit approximations

### Critical Implementation Issues
- **`src/potato_stub.f90`**: Entire stub module returns fake values (**BLOCKING**)
- **`src/potato_wrapper.f90`**: Real POTATO calls commented out (line 74) (**BLOCKING**)
- **`src/potato_field_bridge.f90`**: Real `find_bounce` call commented out (line 361) (**BLOCKING**)
- **`src/thick_orbit_drift.f90`**: Uses simplified estimates, not real POTATO bounce times (**BLOCKING**)
- **`src/transport_thick.f90`**: Falls back to thin orbit approximation (**BLOCKING**)
- **EFIT Integration**: Not implemented for realistic magnetic fields (**BLOCKING**)

### Performance and Physics Limitations
- **Current Status**: No performance data available - no real thick orbit calculations implemented
- **Expected**: POTATO integration ~10x slower than thin orbits due to direct ODE integration
- **Spline Optimization**: Cannot be used for thick orbits (requires particle-by-particle integration)
- **Memory Usage**: Unknown - no real POTATO integration to measure