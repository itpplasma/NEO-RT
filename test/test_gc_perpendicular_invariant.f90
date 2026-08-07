program test_gc_perpendicular_invariant
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: &
        GC_CYL_SUCCESS, gc_buchholz_jk_from_state, &
        gc_cylindrical_field_sample_t, gc_cylindrical_state_t, &
        energy_from_state, gc_potato_jperp_from_state, &
        make_gc_cylindrical_field_sample
    use neort_gc_perpendicular_invariant, only: &
        gc_buchholz_jk_from_mu_phys, gc_buchholz_jk_from_vperp, &
        gc_mu_phys_from_buchholz_jk, gc_mu_phys_from_potato_jperp, &
        gc_mu_phys_from_vperp, &
        gc_perpendicular_specific_energy_from_buchholz_jk, &
        gc_perpendicular_specific_energy_from_mu_phys, &
        gc_potato_jperp_from_mu_phys, gc_potato_jperp_from_vperp, &
        gc_vperp_squared_from_buchholz_jk, gc_vperp_squared_from_mu_phys, &
        gc_vperp_squared_from_potato_jperp
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: mass = 3.7_dp
    real(dp), parameter :: charge = -2.1_dp
    real(dp), parameter :: c_light = 2.8e3_dp
    real(dp), parameter :: bmod = 4.2_dp
    real(dp), parameter :: v_perp = 1.7_dp
    real(dp), parameter :: v_parallel = -0.8_dp
    real(dp), parameter :: v0 = 2.3_dp
    real(dp), parameter :: potential = 0.13_dp
    real(dp), parameter :: radius = 1.8_dp

    type(gc_cylindrical_field_sample_t) :: field
    type(gc_cylindrical_state_t) :: state
    real(dp) :: db_dq(3, 3)
    real(dp) :: mu_ref, h_ref, j_k_ref, omega_c, e_ref, jhat
    real(dp) :: potato_ref, speed, xi, j_k_reversed, mu_reversed
    real(dp) :: potato_from_mu, mu_from_potato
    integer :: status

    mu_ref = mass*v_perp**2/(2.0_dp*bmod)
    h_ref = 0.5_dp*mass*(v_parallel**2 + v_perp**2) &
        +charge*potential
    j_k_ref = mass**2*c_light*v_perp**2/(2.0_dp*abs(charge)*bmod)
    omega_c = abs(charge)*bmod/(mass*c_light)
    e_ref = 0.5_dp*mass*v0**2
    jhat = abs(charge)*j_k_ref/(mass*c_light*e_ref)
    speed = sqrt(v_parallel**2 + v_perp**2)
    xi = v_parallel/speed
    potato_ref = (speed/v0)**2*(1.0_dp - xi**2)/bmod

    db_dq = 0.0_dp
    call make_gc_cylindrical_field_sample(radius, [0.0_dp, 0.0_dp, bmod], &
        db_dq, 0.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], field, status)
    call require(status == GC_CYL_SUCCESS, 'field sample construction failed')

    state = gc_cylindrical_state_t()
    state%R = radius
    state%p_parallel = mass*v_parallel
    state%mu = mu_ref
    call require_close('H decomposition', energy_from_state(state, field, &
        potential, charge, mass), h_ref, 1.0e-13_dp)

    call require_close('mu from v_perp', gc_mu_phys_from_vperp(v_perp, mass, &
        bmod), mu_ref, 1.0e-13_dp)
    call require_close('v_perp^2 from mu', gc_vperp_squared_from_mu_phys(&
        mu_ref, mass, bmod), v_perp**2, 1.0e-13_dp)
    call require_close('Buchholz J_K from v_perp', &
        gc_buchholz_jk_from_vperp(v_perp, mass, bmod, charge, c_light), &
        j_k_ref, 1.0e-13_dp)
    call require_close('Buchholz J_K from mu', &
        gc_buchholz_jk_from_mu_phys(mu_ref, mass, charge, c_light), j_k_ref, &
        1.0e-13_dp)
    call require_close('mu from Buchholz J_K', &
        gc_mu_phys_from_buchholz_jk(j_k_ref, mass, charge, c_light), mu_ref, &
        1.0e-13_dp)
    call require_close('v_perp^2 from Buchholz J_K', &
        gc_vperp_squared_from_buchholz_jk(j_k_ref, charge, c_light, mass, &
        bmod), v_perp**2, 1.0e-13_dp)
    call require_close('J_K omega_c energy', j_k_ref*omega_c, &
        mass*v_perp**2/2.0_dp, 1.0e-13_dp)
    call require_close('J_K omega_c energy identity', j_k_ref*omega_c, &
        0.5_dp*mass*v_perp**2, 1.0e-13_dp)
    call require_close('J_K omega_c divided by mass', j_k_ref*omega_c/mass, &
        0.5_dp*v_perp**2, 1.0e-13_dp)
    call require_close('specific energy from mu', &
        gc_perpendicular_specific_energy_from_mu_phys(mu_ref, bmod, mass), &
        v_perp**2/2.0_dp, 1.0e-13_dp)
    call require_close('specific energy from J_K', &
        gc_perpendicular_specific_energy_from_buchholz_jk(j_k_ref, omega_c, mass), &
        v_perp**2/2.0_dp, 1.0e-13_dp)
    call require_close('Eq. 18 normalized J_K', jhat, &
        (v_perp/v0)**2/bmod, 1.0e-13_dp)
    call require_close('specific H decomposition', h_ref/mass, &
        v_parallel**2/2.0_dp + (j_k_ref*omega_c)/mass &
        +charge*potential/mass, 1.0e-13_dp)

    call require_close('state to Buchholz J_K', &
        gc_buchholz_jk_from_state(state, mass, charge, c_light), j_k_ref, &
        1.0e-13_dp)

    call require_close('POTATO J_perp from v_perp', &
        gc_potato_jperp_from_vperp(v_perp, v0, bmod), potato_ref, &
        1.0e-13_dp)
    potato_from_mu = gc_potato_jperp_from_mu_phys(mu_ref, mass, bmod, v0)
    call require_close('POTATO J_perp from mu', potato_from_mu, potato_ref, &
        1.0e-13_dp)
    call require_close('state to POTATO J_perp', &
        gc_potato_jperp_from_state(state, mass, bmod, v0), potato_ref, &
        1.0e-13_dp)
    call require_close('v_perp^2 from POTATO J_perp', &
        gc_vperp_squared_from_potato_jperp(potato_ref, v0, bmod), &
        v_perp**2, 1.0e-13_dp)
    mu_from_potato = gc_mu_phys_from_potato_jperp(potato_ref, mass, bmod, v0)
    call require_close('POTATO round-trip mu', mu_from_potato, mu_ref, &
        1.0e-13_dp)

    j_k_reversed = gc_buchholz_jk_from_mu_phys(mu_ref, mass, -charge, c_light)
    mu_reversed = gc_mu_phys_from_buchholz_jk(j_k_ref, mass, -charge, c_light)
    call require_close('charge reversal leaves J_K magnitude', j_k_reversed, &
        j_k_ref, 1.0e-13_dp)
    call require_close('charge reversal inverse mu', mu_reversed, mu_ref, &
        1.0e-13_dp)

    write (*, '(A)') 'test_gc_perpendicular_invariant OK'
    call pass_test

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop trim(message)
    end subroutine require

    subroutine require_close(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(reference)))
        if (abs(actual - reference) > tolerance*scale) then
            write (*, '(A,2(1X,ES24.16),1X,ES12.4)') trim(label), actual, &
                reference, tolerance
            error stop 'perpendicular-invariant oracle failed'
        end if
    end subroutine require_close

end program test_gc_perpendicular_invariant
