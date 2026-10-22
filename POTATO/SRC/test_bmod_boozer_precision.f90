program test_bmod_boozer_precision
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use potato_input_mod, only: read_potato_input, n_tor
    use boozer_perturbation_mod, only: boozer_geometry_file
    use libneo_coordinates_base, only: coordinate_system_t
    use libneo_coordinates_chartmap, only: make_chartmap_coordinate_system
    use do_magfie_pert_mod, only: do_magfie_pert_amp, init_magfie_pert_at_s
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    class(coordinate_system_t), allocatable :: chart
    real(dp) :: u(3), xyz(3), r, phi, worst
    complex(dp) :: got, amplitude, expected
    integer :: j

    call read_potato_input('potato.in')
    call make_chartmap_coordinate_system(chart, trim(boozer_geometry_file))
    ! These are forward-map points, so a strict inverse exists independently of
    ! the cold solver's step-size exit. The sqrt(s) chart is steep near the axis.
    call chart%evaluate_cart([0.01_dp, 0.0_dp, 0.13_dp], xyz)
    call bmod_pert(sqrt(sum(xyz(1:2)**2)), xyz(3), got)
    worst = 0.0_dp
    !$omp parallel private(u,xyz,r,phi,got,amplitude,expected)
    call init_magfie_pert_at_s()
    !$omp do reduction(max:worst)
    do j = 0, 511
        u = [0.001_dp + 0.019_dp*real(mod(j, 31), dp)/30.0_dp, &
             2*pi*real(mod(j, 64), dp)/64.0_dp, 0.13_dp]
        call chart%evaluate_cart(u, xyz)
        r = sqrt(sum(xyz(1:2)**2))
        phi = atan2(xyz(2), xyz(1))
        call bmod_pert(r, xyz(3), got)
        call do_magfie_pert_amp([u(1), 0.0_dp, u(2)], amplitude)
        expected = amplitude*exp(cmplx(0.0_dp, n_tor*(u(3)-phi), dp))
        worst = max(worst, abs(got-expected)/max(1.0_dp, abs(expected)))
    end do
    !$omp end do
    !$omp end parallel
    print *, 'maximum field error after strict inverse:', worst
    if (worst > 1.0e-9_dp) error stop 'strict Boozer inverse field mismatch'
end program test_bmod_boozer_precision
