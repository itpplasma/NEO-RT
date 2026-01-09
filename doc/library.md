---
title: Library Interface
---

# Using NEO-RT as a Library

NEO-RT can be used as a library for embedding NTV transport calculations into other codes. The `neort_lib` module provides a clean API for initialization, profile setup, and thread-safe computation across flux surfaces.

## API Overview

The library follows a three-phase workflow:

1. **Initialization** - Load configuration and magnetic field data (once at startup)
2. **Spline preparation** - Set up plasma profile interpolation (once, or per time step if profiles evolve)
3. **Computation** - Calculate transport at flux surfaces (thread-safe, parallelizable)

## Two Operation Modes

### Spline Mode

Use when computing transport across multiple flux surfaces with interpolated profiles:

- Requires `plasma.in` and `profile.in` files (or equivalent array data)
- Profiles (`M_t`, `vth`, etc.) interpolated from splines at each `s` value
- Typically used with `comptorque = .true.`
- Use cases: time evolution coupling, flux surface scans

### Config-Only Mode

Use for single-point calculations where all values come from the config file or
`config_t` struct:

- No `plasma.in` or `profile.in` needed
- Physics values (`s`, `M_t`, `vth`, `qi`, `mi`) read from config file or `config_t`
- Typically used with `comptorque = .false.`
- Use cases: ripple tests, quick single-point calculations

## Module Import

```fortran
use neort_lib
```

This exports:

| Symbol | Description |
| --- | --- |
| `config_t` | Configuration struct (mirrors the `.in` params namelist) |
| `transport_data_t` | Result data type |
| `neort_init` | Load config and magnetic field data (file or struct) |
| `neort_prepare_splines` | Build profile splines (files or arrays) |
| `neort_compute_at_s` | Compute at flux surface (spline mode) |
| `neort_compute_no_splines` | Compute using config values only |

## Subroutine Reference

### neort_init

```fortran
subroutine neort_init(config_file_base, boozer_file, boozer_pert_file)
    character(len=*), intent(in) :: config_file_base   ! Config file without .in extension
    character(len=*), intent(in) :: boozer_file        ! Boozer coordinate file
    character(len=*), intent(in), optional :: boozer_pert_file  ! Perturbation file
end subroutine

subroutine neort_init(config, boozer_file, boozer_pert_file)
    type(config_t), intent(in) :: config               ! In-memory config parameters
    character(len=*), intent(in) :: boozer_file        ! Boozer coordinate file
    character(len=*), intent(in), optional :: boozer_pert_file  ! Perturbation file
end subroutine
```

Overloaded initializer. The file-based form reads the control file, then loads
Boozer magnetic field data (and optional perturbations). The struct-based form
uses `config_t` instead of a file. Call **once** at startup before any parallel
region.

### neort_prepare_splines (file-based)

```fortran
subroutine neort_prepare_splines(plasma_file, profile_file)
    character(len=*), intent(in) :: plasma_file   ! e.g., "plasma.in"
    character(len=*), intent(in) :: profile_file  ! e.g., "profile.in"
end subroutine
```

Reads plasma and rotation profile data from files and builds spline interpolants. Call from the **main thread** before parallel computation.

### neort_prepare_splines (array-based)

```fortran
subroutine neort_prepare_splines(nplasma, am1, am2, Z1, Z2, plasma_data, profile_data)
    integer, intent(in) :: nplasma                  ! Number of plasma data points
    real(8), intent(in) :: am1, am2                 ! Ion mass numbers (species 1, 2)
    real(8), intent(in) :: Z1, Z2                   ! Ion charges (species 1, 2)
    real(8), intent(in) :: plasma_data(:, :)        ! Shape (nplasma, 6): s, n1, n2, T1, T2, E_r
    real(8), intent(in) :: profile_data(:, :)       ! Shape (n, 2): s, M_t
end subroutine
```

Alternative to file-based initialization. Pass plasma profiles directly as arrays. Useful when profiles are computed externally.

#### Notes on overwrites

- The spline preparation routines override `M_t`, `vth`, and the particle charges/masses (`qi`, `mi`) with values derived from the supplied plasma/profile data.
- Any values set previously through `config_t` or namelist input are replaced once spline mode is engaged, so set them back only if you later switch to config-only runs.

### neort_compute_at_s

```fortran
subroutine neort_compute_at_s(s_val, transport_data_out)
    real(8), intent(in) :: s_val                               ! Normalized toroidal flux
    type(transport_data_t), intent(out) :: transport_data_out  ! Results
end subroutine
```

Computes transport coefficients and torque at the specified flux surface. **Thread-safe**: can be called in parallel over different `s` values. Requires prior call to `neort_prepare_splines`.

### neort_compute_no_splines

```fortran
subroutine neort_compute_no_splines(transport_data_out)
    type(transport_data_t), intent(out) :: transport_data_out
end subroutine
```

Computes transport using `s`, `M_t`, `vth`, etc. from the config file or `config_t`. Use when `plasma.in` and `profile.in` are not available.

## Data Types

### config_t

```fortran
type :: config_t
    ! Same fields as the params namelist in <base>.in
    real(8) :: s, qs, ms
    real(8) :: epsmn, mph
    real(8) :: bfac, efac
    real(8) :: M_t, vth
    integer :: m0, inp_swi, vsteps, log_level
    logical :: comptorque, magdrift, nopassing, noshear, pertfile, nonlin
end type
```

### transport_data_t

```fortran
type :: transport_data_t
    type(transport_summary_t) :: summary                     ! Summed transport coefficients
    type(torque_summary_t) :: torque                         ! Summed torque density
    type(transport_harmonic_t), allocatable :: harmonics(:)  ! Per-harmonic data
end type
```

### torque_summary_t

```fortran
type :: torque_summary_t
    logical :: has_torque  ! Torque calculation enabled
    real(8) :: s           ! Normalized toroidal flux
    real(8) :: dVds        ! dV/ds (flux surface volume derivative)
    real(8) :: M_t         ! Toroidal Mach number
    real(8) :: Tco         ! Torque from co-passing particles
    real(8) :: Tctr        ! Torque from counter-passing particles
    real(8) :: Tt          ! Torque from trapped particles
end type
```

### transport_summary_t

```fortran
type :: transport_summary_t
    real(8) :: M_t          ! Toroidal Mach number
    real(8) :: Dco(2)       ! Transport coefficients (D11, D12) for co-passing
    real(8) :: Dctr(2)      ! Transport coefficients (D11, D12) for counter-passing
    real(8) :: Dt(2)        ! Transport coefficients (D11, D12) for trapped
end type
```

## Usage Examples

### Parallel Flux Surface Scan

```fortran
program batch_calculation
    use neort_lib
    implicit none

    type(transport_data_t), allocatable :: results(:)
    real(8), allocatable :: s_array(:)
    integer :: i, n

    n = 100
    allocate(s_array(n), results(n))

    do i = 1, n
        s_array(i) = real(i, 8) / (n + 1)
    end do

    ! Phase 1: Initialize
    call neort_init("driftorbit", "in_file", "in_file_pert")

    ! Phase 2: Prepare splines
    call neort_prepare_splines("plasma.in", "profile.in")

    ! Phase 3: Parallel computation
    !$omp parallel do schedule(dynamic)
    do i = 1, n
        call neort_compute_at_s(s_array(i), results(i))
    end do

    ! Process results...
end program
```

### Time Evolution with Profile Updates

```fortran
subroutine run_time_evolution()
    use neort_lib
    implicit none

    type(transport_data_t), allocatable :: transport_data(:)
    real(8), allocatable :: s_tor(:), plasma_data(:,:), profile_data(:,:)
    integer :: time_step, s_idx, n_s, n_time

    call neort_init("config", "boozer_file", "boozer_pert_file")

    do time_step = 1, n_time
        ! Update profiles from external state
        call prepare_plasma_data(plasma_data)
        call prepare_profile_data(profile_data)

        ! Re-prepare splines with updated profiles
    call neort_prepare_splines(size(plasma_data,1), am1, am2, Z1, Z2, &
                                   plasma_data, profile_data)

        !$omp parallel do schedule(dynamic)
        do s_idx = 1, n_s
            call neort_compute_at_s(s_tor(s_idx), transport_data(s_idx))
        end do
    end do
end subroutine
```

### Single-Point Calculation

```fortran
program single_point
    use neort_lib
    implicit none

    type(transport_data_t) :: result

    call neort_init("ripple_test", "in_file")
    call neort_compute_no_splines(result)

    print *, "D11:", result%summary%Dt(1)
    print *, "D12:", result%summary%Dt(2)
end program
```

## CMake Integration

Link against the `neort` library target:

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_project Fortran)

add_subdirectory(path/to/neo-rt neo-rt)

add_executable(my_program main.f90)
target_link_libraries(my_program PRIVATE neort)
```

## Thread Safety

- `neort_init` and `neort_prepare_splines` must be called from the **main thread** before any parallel region
- `neort_compute_at_s` is **thread-safe** and can be called concurrently
- Shared data (splines, magnetic field) is read-only after initialization
