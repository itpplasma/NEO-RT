module boozer_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use generated_boozer_jxb_kernel, only: boozer_jxb_kernel
    implicit none
    private

    public :: boozer_local_toroidal_torque

contains

    pure subroutine boozer_local_toroidal_torque( &
            covariant_metric, jacobian, b_covariant, weighted_current, &
            weighted_magnetic, torque_density)
        real(dp), intent(in) :: covariant_metric(3, 3), jacobian
        complex(dp), intent(in) :: b_covariant(3), weighted_current(3)
        complex(dp), intent(out) :: weighted_magnetic(3)
        real(dp), intent(out) :: torque_density
        real(dp) :: cs_re, cs_im, ctheta_re, ctheta_im

        call boozer_jxb_kernel( &
            1.0_dp/jacobian, &
            covariant_metric(1, 1), covariant_metric(1, 2), &
            covariant_metric(1, 3), covariant_metric(2, 2), &
            covariant_metric(2, 3), covariant_metric(3, 3), &
            real(b_covariant(1), dp), aimag(b_covariant(1)), &
            real(b_covariant(2), dp), aimag(b_covariant(2)), &
            real(b_covariant(3), dp), aimag(b_covariant(3)), &
            real(weighted_current(1), dp), aimag(weighted_current(1)), &
            real(weighted_current(3), dp), aimag(weighted_current(3)), &
            cs_re, cs_im, ctheta_re, ctheta_im, torque_density)
        weighted_magnetic = cmplx(0.0_dp, 0.0_dp, dp)
        weighted_magnetic(1) = cmplx(cs_re, cs_im, dp)
        weighted_magnetic(3) = cmplx(ctheta_re, ctheta_im, dp)
    end subroutine boozer_local_toroidal_torque

end module boozer_jxb
