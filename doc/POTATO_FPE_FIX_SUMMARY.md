# POTATO Floating Point Exception Fix Summary

## Problem
The POTATO orbit integration was experiencing floating point exceptions (FPE) when calling `find_bounce` through the real integration path. The error occurred in the `rkck_` (Runge-Kutta Cash-Karp) routine at line 208 of `odeint_allroutines.f`.

## Root Causes Identified

1. **Division by zero in velocity calculations** - The `velo` routine had multiple potential division by zero issues:
   - Division by `gamma` (relativistic factor)
   - Division by `bmod` (magnetic field magnitude)
   - Division by `sqrtg` (Jacobian)
   - Division by `hpstar` (effective field denominator)
   - Division by `p` (momentum)
   - Division by `vpa` (parallel velocity)

2. **Uninitialized ODE integration arrays** - The POTATO odeint module arrays were not properly initialized

3. **Aggressive time stepping** - The adaptive time step was too large for some edge cases

4. **Missing input validation** - No checks for valid phase space coordinates before integration

## Fixes Implemented

### 1. Safe Velocity Function (`velo_simple` in `potato_field_bridge.f90`)
- Added safety parameters to prevent division by zero:
  ```fortran
  SMALL_P = 1.0d-10
  SMALL_BMOD = 1.0d-10
  SMALL_SQRTG = 1.0d-10
  SMALL_GAMMA = 1.0d-10
  SMALL_HPSTAR = 1.0d-10
  SMALL_VPA = 1.0d-10
  ```
- Enforced minimum values for all denominators
- Added bounds checking for pitch angle cosine (`alambd`)
- Protected all division operations with safety checks

### 2. Conservative Adaptive Time Stepping
- Changed from aggressive to conservative time stepping:
  - Base time step: `1.0d-5` (was `1.0d-4`)
  - Maximum time step: `1.0d-4` (was `1.0d-3`)
  - Minimum time step: `1.0d-8` (was `1.0d-6`)
- Added pitch-angle dependent scaling:
  - Near turning points (λ < 0.1): use 10% of normal step
  - Deeply trapped particles (λ < 0.3): use 30% of normal step
  - Passing particles: use normal step

### 3. Proper Module Initialization
- Initialize all ODE integration arrays:
  ```fortran
  allocate(ak2(ndim_max), ak3(ndim_max), ak4(ndim_max), ...)
  ak2 = 0.0d0; ak3 = 0.0d0; ...
  ```
- Initialize electric potential polynomial coefficients
- Set conservative gyroradius `ro0 = 1.0d-3`

### 4. Input Validation
- Check momentum is positive: `z_eqm(4) > 0`
- Check pitch cosine is valid: `|z_eqm(5)| ≤ 1`
- Add diagnostic output for debugging

### 5. Test Suite (`test_potato_floating_point_fix.f90`)
Created comprehensive tests for edge cases:
- Very small/large velocities
- Extreme pitch angles (near 0 and π/2)
- Zero velocity components
- Small magnetic field regions

## Usage

To use the fixed POTATO integration:

```fortran
use potato_field_bridge, only: real_find_bounce_calculation
real(8) :: v, eta, taub, delphi
logical :: success

! Initialize field first
call initialize_potato_field(success)

! Calculate bounce time with FPE protection
call real_find_bounce_calculation(v, eta, taub, delphi, success)

if (.not. success) then
    ! Handle error gracefully
end if
```

## Testing

Run the floating point fix test:
```bash
make USE_THICK_ORBITS=ON
ctest -R test_potato_floating_point_fix
```

## Future Improvements

1. Implement proper Poincaré cut initialization (currently using simplified version)
2. Add more sophisticated error recovery in ODE integration
3. Implement variable tolerance based on particle parameters
4. Add field line following tests for extreme geometries