module boozer_ampere
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use generated_boozer_ampere_kernel, only: boozer_ampere_kernel
    implicit none
    private

    public :: boozer_harmonic_current

contains

    pure subroutine boozer_harmonic_current( &
            signed_jacobian, mu0, m, n, b_covariant, db_covariant_ds, &
            j_contravariant)
        ! Ampere's law for B_hat(s)*exp(i*(n*phi + m*theta)).
        !
        ! Coordinates are ordered (s, phi, theta).  B is covariant, J is
        ! contravariant, and signed_jacobian must use that exact coordinate
        ! orientation.  Inputs must be in one consistent unit system; SI uses
        ! mu0 = 4*pi*1e-7 H/m.
        real(dp), intent(in) :: signed_jacobian, mu0
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: b_covariant(3), db_covariant_ds(3)
        complex(dp), intent(out) :: j_contravariant(3)
        real(dp) :: js_re, js_im, jphi_re, jphi_im, jtheta_re, jtheta_im
        real(dp) :: scale

        if (signed_jacobian == 0.0_dp) error stop "zero Boozer Jacobian"
        if (mu0 == 0.0_dp) error stop "zero magnetic permeability"
        scale = 1.0_dp/(mu0*signed_jacobian)
        call boozer_ampere_kernel( &
            scale, real(m, dp), real(n, dp), &
            real(b_covariant(1), dp), aimag(b_covariant(1)), &
            real(b_covariant(2), dp), aimag(b_covariant(2)), &
            real(b_covariant(3), dp), aimag(b_covariant(3)), &
            real(db_covariant_ds(2), dp), aimag(db_covariant_ds(2)), &
            real(db_covariant_ds(3), dp), aimag(db_covariant_ds(3)), &
            js_re, js_im, jphi_re, jphi_im, jtheta_re, jtheta_im)
        j_contravariant(1) = cmplx(js_re, js_im, dp)
        j_contravariant(2) = cmplx(jphi_re, jphi_im, dp)
        j_contravariant(3) = cmplx(jtheta_re, jtheta_im, dp)
    end subroutine boozer_harmonic_current

end module boozer_ampere
