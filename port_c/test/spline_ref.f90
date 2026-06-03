! Dumps reference spline coefficients and evaluations from the real itpplasma
! spline module, for cross-checking the C port. Output: machine-readable lines.
program spline_ref
    use spline, only: spline_coeff, spline_val_0
    implicit none
    integer, parameter :: np = 9
    real(8) :: x(np), y(np), coeff(np - 1, 5), val(3)
    real(8) :: xe
    integer :: i, k

    do i = 1, np
        x(i) = real(i - 1, 8) * 0.37d0
        y(i) = sin(1.3d0 * x(i)) + 0.5d0 * x(i)
    end do

    coeff = spline_coeff(x, y)

    do i = 1, np - 1
        write (*, '("COEFF",1X,I0,5(1X,ES24.16E3))') i - 1, (coeff(i, k), k=1, 5)
    end do

    do i = 0, 40
        xe = -0.2d0 + real(i, 8) * (x(np) + 0.4d0) / 40.0d0
        val = spline_val_0(coeff, xe)
        write (*, '("VAL",1X,ES24.16E3,3(1X,ES24.16E3))') xe, val(1), val(2), val(3)
    end do
end program spline_ref
