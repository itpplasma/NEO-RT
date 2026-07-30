module response_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: cylindrical_toroidal_torque
    public :: mars_surface_torque

contains

    pure function cylindrical_toroidal_torque(radius, weight, jr, jz, br, bz) result(torque)
        real(dp), intent(in) :: radius(:), weight(:)
        complex(dp), intent(in) :: jr(:), jz(:), br(:), bz(:)
        real(dp) :: torque

        torque = 0.5_dp*sum(weight*radius*real(jz*conjg(br) - jr*conjg(bz), dp))
    end function cylindrical_toroidal_torque

    pure function mars_surface_torque(j1, j2, b1, b2) result(torque)
        complex(dp), intent(in) :: j1(:), j2(:), b1(:), b2(:)
        real(dp) :: torque

        torque = 0.5_dp*sum(real(conjg(j1)*b2 - conjg(j2)*b1, dp))
    end function mars_surface_torque

end module response_jxb
