program test_chartmap_transport
    ! Behavioral test: resonant transport coefficients from the Boozer chartmap
    ! input (inp_swi=10) against the .bc input (inp_swi=8) for the same field.
    !
    ! The chartmap fixture test/fixtures/circ_chartmap.nc was generated from
    ! examples/circ.bc. Both runs use the analytic perturbation harmonic
    ! (pertfile=.false.), so they differ only in the equilibrium representation.
    !
    ! D11 is normalized with dsdreff = 2*sqrt(s)/a, i.e. proportional to a**2.
    ! The chartmap minor radius is the equal-area edge radius, which differs from
    ! the .bc header value by about 1e-3, so the .bc oracle is rescaled by
    ! (a_chartmap/a_bc)**2 before comparison. The remaining difference comes from
    ! the rho spline and theta Fourier interpolation of the chartmap tables.
    !
    ! Regressions caught: a = 0 from the chartmap reader made D identically zero;
    ! linear theta interpolation of Bmod biased D11 by about 6 percent.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: a, magfie_thread_init
    use logger, only: set_log_level
    use neort_config, only: config_t
    use neort_datatypes, only: transport_data_t
    use neort_lib, only: neort_init, neort_compute_no_splines

    implicit none

    real(dp), parameter :: tol_rel = 2.0e-3_dp
    type(transport_data_t) :: bc_result, cm_result
    real(dp) :: a_bc, a_cm, d11_bc(3), d11_cm(3), d11_bc_total, d11_cm_total
    character(len=512) :: chartmap_path, bc_path
    integer :: stat, k

    call get_environment_variable("CHARTMAP_FILE", chartmap_path, status=stat)
    if (stat /= 0 .or. len_trim(chartmap_path) == 0) then
        error stop "test_chartmap_transport: CHARTMAP_FILE required"
    end if
    call get_environment_variable("BC_FILE", bc_path, status=stat)
    if (stat /= 0 .or. len_trim(bc_path) == 0) then
        error stop "test_chartmap_transport: BC_FILE required"
    end if

    ! Chartmap first: module state such as the minor radius must not be
    ! inherited from a previous .bc read.
    call run_case(10, trim(chartmap_path), cm_result)
    a_cm = a
    ! Reset the per-thread field buffers, which the chartmap path leaves
    ! unallocated, before switching to the .bc representation.
    call magfie_thread_init()
    call run_case(8, trim(bc_path), bc_result)
    a_bc = a

    d11_bc = [bc_result%summary%Dco(1), bc_result%summary%Dctr(1), &
              bc_result%summary%Dt(1)] * (a_cm / a_bc)**2
    d11_cm = [cm_result%summary%Dco(1), cm_result%summary%Dctr(1), &
              cm_result%summary%Dt(1)]
    d11_bc_total = sum(d11_bc)
    d11_cm_total = sum(d11_cm)

    print *, "a [cm]       bc =", a_bc, " chartmap =", a_cm
    print *, "D11 (co, ctr, t), bc rescaled:", d11_bc
    print *, "D11 (co, ctr, t), chartmap   :", d11_cm

    if (.not. (d11_cm_total > 0.0_dp)) then
        print *, "ERROR: chartmap D11 must be positive, got", d11_cm_total
        error stop "test_chartmap_transport: chartmap D11 not positive"
    end if
    if (.not. (d11_bc_total > 0.0_dp)) then
        error stop "test_chartmap_transport: .bc D11 not positive"
    end if

    do k = 1, 3
        call assert_close("D11 component", d11_cm(k), d11_bc(k))
    end do
    call assert_close("D11 total", d11_cm_total, d11_bc_total)

    print *, "test_chartmap_transport OK: chartmap D11 matches .bc to", tol_rel

contains

    subroutine run_case(input_switch, boozer_file, result)
        integer, intent(in) :: input_switch
        character(len=*), intent(in) :: boozer_file
        type(transport_data_t), intent(out) :: result
        type(config_t) :: config

        config%s = 0.5_dp
        config%M_t = 0.036_dp
        config%qs = 1.0_dp
        config%ms = 2.014_dp
        config%vth = 4.0e7_dp
        config%epsmn = 1.0e-3_dp
        config%m0 = 0
        config%mph = 3
        config%magdrift = .true.
        config%noshear = .true.
        config%pertfile = .false.
        config%inp_swi = input_switch
        config%vsteps = 256
        config%log_level = -1

        call neort_init(config, boozer_file)
        call set_log_level(-1)
        call neort_compute_no_splines(result)
    end subroutine run_case

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: rel_err

        rel_err = abs(actual - expected) / abs(expected)
        if (rel_err > tol_rel) then
            print *, "ERROR: ", label, " chartmap =", actual, " bc =", expected, &
                " rel_err =", rel_err
            error stop "test_chartmap_transport: chartmap differs from .bc"
        end if
    end subroutine assert_close

end program test_chartmap_transport
