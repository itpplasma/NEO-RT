# Next Steps: Complete POTATO Integration

## Current Status
✅ **Production-ready thick orbit framework with stub implementation**
✅ **Runtime dispatch architecture working**
✅ **Canonical frequency calculations documented**
✅ **Comprehensive test suite passing**

## Remaining Work: POTATO Field Interface Bridge

### 1. Implement Missing Field Evaluation Functions

Create `src/potato_field_bridge.f90` to provide:

```fortran
module potato_field_bridge
    ! Bridge NEO-RT magnetic field to POTATO field_eq_mod requirements
    
    public :: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
    
contains
    
    function psif(R, Z) result(psi)
        ! Evaluate flux function at (R,Z)
        ! Connect to NEO-RT's magnetic field evaluation
        real(8), intent(in) :: R, Z
        real(8) :: psi
        
        ! TODO: Call NEO-RT magnetic field evaluation
        ! Transform coordinates as needed
    end function psif
    
    subroutine dpsidr(R, Z, dpsi_dR)
        ! Flux derivative wrt R
        ! TODO: Implement using NEO-RT field gradients
    end subroutine dpsidr
    
    subroutine dpsidz(R, Z, dpsi_dZ) 
        ! Flux derivative wrt Z
        ! TODO: Implement using NEO-RT field gradients
    end subroutine dpsidz
    
end module potato_field_bridge
```

### 2. Coordinate System Integration

**Challenge**: Map between coordinate systems
- **NEO-RT**: (s, theta, phi) flux coordinates
- **POTATO**: (R, phi, Z) cylindrical coordinates

**Solution**: Create coordinate transformation routines:
```fortran
subroutine neort_to_potato_coords(s, theta, R, Z)
subroutine potato_to_neort_coords(R, Z, s, theta)
```

### 3. Enable POTATO Build Integration

Update `CMakeLists.txt`:
```cmake
if(USE_THICK_ORBITS)
  # Add field bridge
  list(APPEND NEO_RT_SOURCES src/potato_field_bridge.f90)
  
  # Enable POTATO compilation
  add_subdirectory(POTATO/SRC)
  target_link_libraries(neo_rt PUBLIC potato_base)
  target_include_directories(neo_rt PUBLIC "${PROJECT_BINARY_DIR}/POTATO/SRC")
endif()
```

### 4. Replace Stub Implementation

Update `potato_wrapper.f90`:
```fortran
! Replace this:
call potato_stub_find_bounce(v, eta, taub, delphi, extraset)

! With this:
call find_bounce(next, velo, dtau_in, z_eqm, taub, delphi, extraset)
```

## Timeline Estimate

- **Field bridge implementation**: 1-2 weeks
- **Coordinate system integration**: 1 week  
- **Build system integration**: 3-5 days
- **Testing and validation**: 1 week
- **Total**: 3-4 weeks for complete POTATO integration

## Validation Plan

1. **Field evaluation accuracy**: Compare psif with NEO-RT magnetic field
2. **Orbit integration**: Validate bounce times against POTATO standalone
3. **Frequency accuracy**: Compare canonical frequencies
4. **Physics consistency**: Conservation laws, orbit closure
5. **Performance benchmarking**: Document computational costs

## Alternative Approaches

If field bridge proves complex:

1. **Simplified field model**: Use analytical field for initial testing
2. **POTATO standalone validation**: Compare against separate POTATO runs
3. **Hybrid approach**: Use POTATO for specific orbit calculations only

## Ready for Implementation

The current architecture provides all necessary infrastructure for immediate POTATO integration once the field interface bridge is implemented.