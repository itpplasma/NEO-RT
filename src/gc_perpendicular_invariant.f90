module neort_gc_perpendicular_invariant
    !! Explicit conversions for the perpendicular guiding-center invariant.
    !!
    !! The canonical dimensional quantity used by the cylindrical backend is
    !! the physical magnetic moment
    !!
    !!     mu_phys = m v_perp**2 / (2 B),
    !!
    !! with units mass*velocity**2/magnetic_field (energy/gauss in the
    !! repository's CGS convention).  It is independent of the sign of the
    !! particle charge.  No normalized or code-specific quantity is stored in
    !! place of mu_phys.
    !!
    !! Buchholz et al. (2022), Eq. 4, uses
    !!
    !!     J_K = c mu_phys / abs(q),
    !!     J_K omega_c = v_perp**2 / 2,
    !!     omega_c = abs(q) B / (m c).
    !!
    !! POTATO's documented quantity is instead
    !!
    !!     J_potato = p**2 (1-xi**2) / B
    !!               = v_perp**2 / (v0**2 B),
    !!
    !! where p=v/v0 and xi=v_parallel/v.  The v0 and B factors are therefore
    !! explicit at every conversion boundary.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    character(len=*), parameter, public :: GC_PERP_MU_PHYS_UNITS = &
        'mass*velocity^2/magnetic_field'
    character(len=*), parameter, public :: GC_PERP_BUCHHOLZ_JK_UNITS = &
        'c*mu_phys/abs(charge)'
    character(len=*), parameter, public :: GC_PERP_POTATO_JPERP_UNITS = &
        'p^2*(1-xi^2)/B'

    public :: gc_buchholz_jk_from_mu_phys
    public :: gc_buchholz_jk_from_vperp
    public :: gc_mu_phys_from_buchholz_jk
    public :: gc_mu_phys_from_potato_jperp
    public :: gc_mu_phys_from_vperp
    public :: gc_perpendicular_specific_energy_from_buchholz_jk
    public :: gc_perpendicular_specific_energy_from_mu_phys
    public :: gc_potato_jperp_from_mu_phys
    public :: gc_potato_jperp_from_vperp
    public :: gc_vperp_squared_from_buchholz_jk
    public :: gc_vperp_squared_from_mu_phys
    public :: gc_vperp_squared_from_potato_jperp

contains

    pure function gc_mu_phys_from_vperp(v_perp, mass, bmod) result(mu_phys)
        real(dp), intent(in) :: v_perp, mass, bmod
        real(dp) :: mu_phys

        mu_phys = mass*v_perp**2/(2.0_dp*bmod)
    end function gc_mu_phys_from_vperp

    pure function gc_vperp_squared_from_mu_phys(mu_phys, mass, bmod) &
            result(v_perp_squared)
        real(dp), intent(in) :: mu_phys, mass, bmod
        real(dp) :: v_perp_squared

        v_perp_squared = 2.0_dp*mu_phys*bmod/mass
    end function gc_vperp_squared_from_mu_phys

    pure function gc_buchholz_jk_from_mu_phys(mu_phys, charge, c_light) &
            result(j_k)
        real(dp), intent(in) :: mu_phys, charge, c_light
        real(dp) :: j_k

        j_k = c_light*mu_phys/abs(charge)
    end function gc_buchholz_jk_from_mu_phys

    pure function gc_mu_phys_from_buchholz_jk(j_k, charge, c_light) &
            result(mu_phys)
        real(dp), intent(in) :: j_k, charge, c_light
        real(dp) :: mu_phys

        mu_phys = abs(charge)*j_k/c_light
    end function gc_mu_phys_from_buchholz_jk

    pure function gc_buchholz_jk_from_vperp(v_perp, mass, bmod, charge, &
            c_light) result(j_k)
        real(dp), intent(in) :: v_perp, mass, bmod, charge, c_light
        real(dp) :: j_k

        j_k = gc_buchholz_jk_from_mu_phys( &
            gc_mu_phys_from_vperp(v_perp, mass, bmod), charge, c_light)
    end function gc_buchholz_jk_from_vperp

    pure function gc_vperp_squared_from_buchholz_jk(j_k, charge, c_light, &
            mass, bmod) result(v_perp_squared)
        real(dp), intent(in) :: j_k, charge, c_light, mass, bmod
        real(dp) :: v_perp_squared

        v_perp_squared = 2.0_dp*gc_mu_phys_from_buchholz_jk( &
            j_k, charge, c_light)*bmod/mass
    end function gc_vperp_squared_from_buchholz_jk

    pure function gc_perpendicular_specific_energy_from_mu_phys(mu_phys, &
            bmod, mass) result(specific_energy)
        real(dp), intent(in) :: mu_phys, bmod, mass
        real(dp) :: specific_energy

        specific_energy = mu_phys*bmod/mass
    end function gc_perpendicular_specific_energy_from_mu_phys

    pure function gc_perpendicular_specific_energy_from_buchholz_jk(j_k, &
            omega_c) result(specific_energy)
        real(dp), intent(in) :: j_k, omega_c
        real(dp) :: specific_energy

        specific_energy = j_k*omega_c
    end function gc_perpendicular_specific_energy_from_buchholz_jk

    pure function gc_potato_jperp_from_vperp(v_perp, v0, bmod) &
            result(j_potato)
        real(dp), intent(in) :: v_perp, v0, bmod
        real(dp) :: j_potato

        j_potato = v_perp**2/(v0**2*bmod)
    end function gc_potato_jperp_from_vperp

    pure function gc_vperp_squared_from_potato_jperp(j_potato, v0, bmod) &
            result(v_perp_squared)
        real(dp), intent(in) :: j_potato, v0, bmod
        real(dp) :: v_perp_squared

        v_perp_squared = j_potato*v0**2*bmod
    end function gc_vperp_squared_from_potato_jperp

    pure function gc_potato_jperp_from_mu_phys(mu_phys, mass, bmod, v0) &
            result(j_potato)
        real(dp), intent(in) :: mu_phys, mass, bmod, v0
        real(dp) :: j_potato

        j_potato = gc_potato_jperp_from_vperp( &
            sqrt(gc_vperp_squared_from_mu_phys(mu_phys, mass, bmod)), v0, bmod)
    end function gc_potato_jperp_from_mu_phys

    pure function gc_mu_phys_from_potato_jperp(j_potato, mass, bmod, v0) &
            result(mu_phys)
        real(dp), intent(in) :: j_potato, mass, bmod, v0
        real(dp) :: mu_phys

        mu_phys = gc_mu_phys_from_vperp( &
            sqrt(gc_vperp_squared_from_potato_jperp(j_potato, v0, bmod)), &
            mass, bmod)
    end function gc_mu_phys_from_potato_jperp

end module neort_gc_perpendicular_invariant
