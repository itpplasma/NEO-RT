module electromagnetic_torque
    ! Algebraic electromagnetic-torque contraction for native complex harmonics.
    !
    ! J and B are stored non-axisymmetric Fourier amplitudes. The factor 1/2
    ! is the real-field average of one complex harmonic; callers must use one
    ! consistent Fourier sign convention for both vectors. Results retain the
    ! input J*B normalization and are not SI torque by themselves. Geometry,
    ! B0**2/mu0 scaling, radial-coordinate conversion, filtering, and edge
    ! policies belong to an adapter outside this coordinate-independent kernel.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: cross_contraction
    public :: staggered_cross_contraction

contains

    pure elemental function cross_contraction(j1, b2, j2, b1) result(value)
        complex(dp), intent(in) :: j1, b2, j2, b1
        real(dp) :: value

        value = 0.5_dp*real(conjg(j1)*b2 - conjg(j2)*b1, dp)
    end function cross_contraction

    function staggered_cross_contraction(j1_half, b2_half, j2_full, b1_full) &
            result(density)
        ! Return the unsmoothed native radial integration density on half cells.
        ! Full-mesh J2*B1 products are averaged at adjacent cell endpoints;
        ! half-mesh J1*B2 products are used at their native locations.
        complex(dp), intent(in) :: j1_half(:, :), b2_half(:, :)
        complex(dp), intent(in) :: j2_full(:, :), b1_full(:, :)
        real(dp) :: density(size(j1_half, 1))
        integer :: cell_count, mode

        call require_same_shape(j1_half, b2_half, "half-mesh B2")
        call require_same_shape(j2_full, b1_full, "full-mesh B1")
        if (size(j2_full, 2) /= size(j1_half, 2)) then
            error stop "full/half meshes have different mode counts"
        end if
        cell_count = size(j1_half, 1)
        if (size(j2_full, 1) /= cell_count + 1) then
            error stop "full mesh must have one more radial point than half mesh"
        end if
        density = 0.0_dp
        do mode = 1, size(j1_half, 2)
            density = density + 0.5_dp*real( &
                conjg(j1_half(:, mode))*b2_half(:, mode) - 0.5_dp*( &
                conjg(j2_full(:cell_count, mode))*b1_full(:cell_count, mode) + &
                conjg(j2_full(2:, mode))*b1_full(2:, mode)), dp)
        end do
    end function staggered_cross_contraction

    subroutine require_same_shape(first, second, label)
        complex(dp), intent(in) :: first(:, :), second(:, :)
        character(*), intent(in) :: label

        if (size(first, 1) /= size(second, 1)) then
            error stop label//" has a different radial size"
        end if
        if (size(first, 2) /= size(second, 2)) then
            error stop label//" has a different mode count"
        end if
    end subroutine require_same_shape

end module electromagnetic_torque
