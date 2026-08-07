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
    !!     J_K = m c mu_phys / abs(q),
    !!     J_K omega_c = mu_phys B = m v_perp**2 / 2,
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
    use neort_buchholz_action_symbolic, only: &
        evaluate_neort_buchholz_action
    use neort_buchholz_energy_symbolic, only: &
        evaluate_neort_buchholz_energy
    use neort_buchholz_specific_energy_symbolic, only: &
        evaluate_neort_buchholz_specific_energy
    use neort_physical_mu_symbolic, only: evaluate_neort_physical_mu
    use neort_potato_mu_symbolic, only: evaluate_neort_potato_mu
    use neort_potato_velocity_symbolic, only: evaluate_neort_potato_velocity

    implicit none
    private

    character(len=*), parameter, public :: GC_PERP_MU_PHYS_UNITS = &
        'mass*velocity^2/magnetic_field'
    character(len=*), parameter, public :: GC_PERP_BUCHHOLZ_JK_UNITS = &
        'mass*c*mu_phys/abs(charge)'
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
        real(dp) :: unused_vperp_squared, unused_specific_energy

        call evaluate_neort_physical_mu(v_perp, mass, bmod, 0.0_dp, &
            mu_phys, unused_vperp_squared, unused_specific_energy)
    end function gc_mu_phys_from_vperp

    pure function gc_vperp_squared_from_mu_phys(mu_phys, mass, bmod) &
            result(v_perp_squared)
        real(dp), intent(in) :: mu_phys, mass, bmod
        real(dp) :: v_perp_squared
        real(dp) :: unused_mu, unused_specific_energy

        call evaluate_neort_physical_mu(0.0_dp, mass, bmod, mu_phys, &
            unused_mu, v_perp_squared, unused_specific_energy)
    end function gc_vperp_squared_from_mu_phys

    pure function gc_buchholz_jk_from_mu_phys(mu_phys, mass, charge, c_light) &
            result(j_k)
        real(dp), intent(in) :: mu_phys, mass, charge, c_light
        real(dp) :: j_k
        real(dp) :: unused_mu

        call evaluate_neort_buchholz_action(mu_phys, 0.0_dp, mass, charge, &
            c_light, j_k, unused_mu)
    end function gc_buchholz_jk_from_mu_phys

    pure function gc_mu_phys_from_buchholz_jk(j_k, mass, charge, c_light) &
            result(mu_phys)
        real(dp), intent(in) :: j_k, mass, charge, c_light
        real(dp) :: mu_phys
        real(dp) :: unused_j_k

        call evaluate_neort_buchholz_action(0.0_dp, j_k, mass, charge, &
            c_light, unused_j_k, mu_phys)
    end function gc_mu_phys_from_buchholz_jk

    pure function gc_buchholz_jk_from_vperp(v_perp, mass, bmod, charge, &
            c_light) result(j_k)
        real(dp), intent(in) :: v_perp, mass, bmod, charge, c_light
        real(dp) :: j_k

        j_k = gc_buchholz_jk_from_mu_phys( &
            gc_mu_phys_from_vperp(v_perp, mass, bmod), mass, charge, c_light)
    end function gc_buchholz_jk_from_vperp

    pure function gc_vperp_squared_from_buchholz_jk(j_k, charge, c_light, &
            mass, bmod) result(v_perp_squared)
        real(dp), intent(in) :: j_k, charge, c_light, mass, bmod
        real(dp) :: v_perp_squared
        real(dp) :: unused_omega_c, unused_specific_energy

        call evaluate_neort_buchholz_energy(j_k, mass, charge, c_light, &
            bmod, unused_omega_c, v_perp_squared, unused_specific_energy)
    end function gc_vperp_squared_from_buchholz_jk

    pure function gc_perpendicular_specific_energy_from_mu_phys(mu_phys, &
            bmod, mass) result(specific_energy)
        real(dp), intent(in) :: mu_phys, bmod, mass
        real(dp) :: specific_energy
        real(dp) :: unused_mu, unused_vperp_squared

        call evaluate_neort_physical_mu(0.0_dp, mass, bmod, mu_phys, &
            unused_mu, unused_vperp_squared, specific_energy)
    end function gc_perpendicular_specific_energy_from_mu_phys

    pure function gc_perpendicular_specific_energy_from_buchholz_jk(j_k, &
            omega_c, mass) result(specific_energy)
        real(dp), intent(in) :: j_k, omega_c, mass
        real(dp) :: specific_energy

        call evaluate_neort_buchholz_specific_energy(j_k, omega_c, mass, &
            specific_energy)
    end function gc_perpendicular_specific_energy_from_buchholz_jk

    pure function gc_potato_jperp_from_vperp(v_perp, v0, bmod) &
            result(j_potato)
        real(dp), intent(in) :: v_perp, v0, bmod
        real(dp) :: j_potato
        real(dp) :: unused_vperp_squared

        call evaluate_neort_potato_velocity(v_perp, 0.0_dp, v0, bmod, &
            j_potato, unused_vperp_squared)
    end function gc_potato_jperp_from_vperp

    pure function gc_vperp_squared_from_potato_jperp(j_potato, v0, bmod) &
            result(v_perp_squared)
        real(dp), intent(in) :: j_potato, v0, bmod
        real(dp) :: v_perp_squared
        real(dp) :: unused_j_potato

        call evaluate_neort_potato_velocity(0.0_dp, j_potato, v0, bmod, &
            unused_j_potato, v_perp_squared)
    end function gc_vperp_squared_from_potato_jperp

    pure function gc_potato_jperp_from_mu_phys(mu_phys, mass, bmod, v0) &
            result(j_potato)
        real(dp), intent(in) :: mu_phys, mass, bmod, v0
        real(dp) :: j_potato
        real(dp) :: unused_mu

        associate (unused_bmod => bmod)
        end associate
        call evaluate_neort_potato_mu(mu_phys, 0.0_dp, mass, v0, &
            j_potato, unused_mu)
    end function gc_potato_jperp_from_mu_phys

    pure function gc_mu_phys_from_potato_jperp(j_potato, mass, bmod, v0) &
            result(mu_phys)
        real(dp), intent(in) :: j_potato, mass, bmod, v0
        real(dp) :: mu_phys
        real(dp) :: unused_j_potato

        associate (unused_bmod => bmod)
        end associate
        call evaluate_neort_potato_mu(0.0_dp, j_potato, mass, v0, &
            unused_j_potato, mu_phys)
    end function gc_mu_phys_from_potato_jperp

end module neort_gc_perpendicular_invariant
