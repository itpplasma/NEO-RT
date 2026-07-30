module boozer_jxb_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap_geometry, only: boozer_chartmap_geometry_t, &
        evaluate_boozer_metric
    use boozer_jxb, only: boozer_local_toroidal_torque
    use boozer_vector_io, only: boozer_vector_harmonics_t, &
        compute_boozer_weighted_current
    implicit none
    private

    public :: compute_boozer_jxb_profile

contains

    subroutine compute_boozer_jxb_profile( &
            geometry, field, mu0, torque_density)
        type(boozer_chartmap_geometry_t), intent(in) :: geometry
        type(boozer_vector_harmonics_t), intent(in) :: field
        real(dp), intent(in) :: mu0
        real(dp), allocatable, intent(out) :: torque_density(:)
        complex(dp), allocatable :: weighted_current(:, :, :)
        real(dp), allocatable :: metric(:, :, :), jacobian(:)
        complex(dp) :: b_covariant(3), current(3), weighted_magnetic(3)
        complex(dp) :: phase
        real(dp) :: local_density
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: mode_index, radial_index, theta_index

        call compute_boozer_weighted_current(field, mu0, weighted_current)
        allocate (torque_density(field%radial_count))
        allocate (metric(3, 3, geometry%theta_count))
        allocate (jacobian(geometry%theta_count))
        torque_density = 0.0_dp
        do radial_index = 1, field%radial_count
            call evaluate_boozer_metric( &
                geometry, field%s(radial_index), metric, jacobian)
            do theta_index = 1, geometry%theta_count
                b_covariant = cmplx(0.0_dp, 0.0_dp, dp)
                current = cmplx(0.0_dp, 0.0_dp, dp)
                do mode_index = 1, field%mode_count
                    phase = exp(cmplx( &
                        0.0_dp, &
                        real(field%modes(mode_index), dp)* &
                        geometry%theta(theta_index), dp))
                    b_covariant = b_covariant + &
                        field%b_covariant(:, radial_index, mode_index)*phase
                    current = current + &
                        weighted_current(:, radial_index, mode_index)*phase
                end do
                call boozer_local_toroidal_torque( &
                    metric(:, :, theta_index), jacobian(theta_index), &
                    b_covariant, current, weighted_magnetic, local_density)
                torque_density(radial_index) = &
                    torque_density(radial_index) + local_density
            end do
        end do
        torque_density = 4.0_dp*pi**2*torque_density/ &
            real(geometry%theta_count, dp)
    end subroutine compute_boozer_jxb_profile

end module boozer_jxb_profile
