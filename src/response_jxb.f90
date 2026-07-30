module response_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: cylindrical_toroidal_torque
    public :: integrate_mars_profile
    public :: mars_half_mesh_torque
    public :: mars_surface_torque
    public :: reconstruct_mars_profile

contains

    function cylindrical_toroidal_torque(radius, weight, jr, jz, br, bz) result(torque)
        real(dp), intent(in) :: radius(:), weight(:)
        complex(dp), intent(in) :: jr(:), jz(:), br(:), bz(:)
        real(dp) :: torque

        call require_size("cylindrical weights", size(radius), size(weight))
        call require_size("cylindrical jr", size(radius), size(jr))
        call require_size("cylindrical jz", size(radius), size(jz))
        call require_size("cylindrical br", size(radius), size(br))
        call require_size("cylindrical bz", size(radius), size(bz))
        torque = 0.5_dp*sum(weight*radius*real(jz*conjg(br) - jr*conjg(bz), dp))
    end function cylindrical_toroidal_torque

    function mars_surface_torque(j1, j2, b1, b2) result(torque)
        complex(dp), intent(in) :: j1(:), j2(:), b1(:), b2(:)
        real(dp) :: torque

        call require_size("MARS j2", size(j1), size(j2))
        call require_size("MARS b1", size(j1), size(b1))
        call require_size("MARS b2", size(j1), size(b2))
        torque = 0.5_dp*sum(real(conjg(j1)*b2 - conjg(j2)*b1, dp))
    end function mars_surface_torque

    function mars_half_mesh_torque( &
            j1_half, b2_half, j2_left, b1_left, j2_right, b1_right) result(torque)
        complex(dp), intent(in) :: j1_half(:), b2_half(:)
        complex(dp), intent(in) :: j2_left(:), b1_left(:), j2_right(:), b1_right(:)
        real(dp) :: torque

        call require_size("MARS half-mesh b2", size(j1_half), size(b2_half))
        call require_size("MARS left j2", size(j1_half), size(j2_left))
        call require_size("MARS left b1", size(j1_half), size(b1_left))
        call require_size("MARS right j2", size(j1_half), size(j2_right))
        call require_size("MARS right b1", size(j1_half), size(b1_right))
        torque = 0.5_dp*sum(real(conjg(j1_half)*b2_half - &
            0.5_dp*(conjg(j2_left)*b1_left + conjg(j2_right)*b1_right), dp))
    end function mars_half_mesh_torque

    function integrate_mars_profile(edges, density, scale) result(torque)
        real(dp), intent(in) :: edges(:), density(:), scale
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: torque
        integer :: number_of_cells

        number_of_cells = size(density)
        if (size(edges) /= number_of_cells + 1) then
            error stop "MARS profile requires one more edge than density value"
        end if
        torque = 4.0_dp*pi**2*scale*sum( &
            density*(edges(2:number_of_cells + 1) - edges(:number_of_cells)))
    end function integrate_mars_profile

    subroutine reconstruct_mars_profile( &
            j1, j2, b1, b2, half_mesh, first_torque_cell, smoothing_passes, &
            edge_cutoff, raw, processed)
        complex(dp), intent(in) :: j1(:, :), j2(:, :), b1(:, :), b2(:, :)
        real(dp), intent(in) :: half_mesh(:), edge_cutoff
        integer, intent(in) :: first_torque_cell, smoothing_passes
        real(dp), allocatable, intent(out) :: raw(:), processed(:)
        real(dp), allocatable :: previous(:)
        integer :: cell, number_of_cells, pass

        number_of_cells = size(half_mesh)
        call check_mars_profile_shapes(j1, j2, b1, b2, number_of_cells)
        if (first_torque_cell < 1) error stop "Invalid first MARS torque cell"
        if (first_torque_cell > number_of_cells) error stop "Invalid first MARS torque cell"
        if (smoothing_passes < 0) error stop "Invalid MARS smoothing pass count"
        allocate (raw(number_of_cells), processed(number_of_cells), previous(number_of_cells))
        raw = 0.0_dp
        do cell = first_torque_cell, number_of_cells
            raw(cell) = mars_half_mesh_torque( &
                j1(cell, :), b2(cell, :), j2(cell, :), b1(cell, :), &
                j2(cell + 1, :), b1(cell + 1, :))
        end do
        processed = raw
        do pass = 1, smoothing_passes
            previous = processed
            processed(2:number_of_cells - 1) = 0.1_dp*previous(:number_of_cells - 2) + &
                0.8_dp*previous(2:number_of_cells - 1) + 0.1_dp*previous(3:)
        end do
        processed(:first_torque_cell - 1) = 0.0_dp
        where (half_mesh > edge_cutoff) processed = 0.0_dp
    end subroutine reconstruct_mars_profile

    subroutine check_mars_profile_shapes(j1, j2, b1, b2, number_of_cells)
        complex(dp), intent(in) :: j1(:, :), j2(:, :), b1(:, :), b2(:, :)
        integer, intent(in) :: number_of_cells
        integer :: number_of_modes

        number_of_modes = size(j1, 2)
        if (size(j1, 1) < number_of_cells + 1) error stop "MARS j1 radial size"
        if (size(j2, 1) < number_of_cells + 1) error stop "MARS j2 radial size"
        if (size(b1, 1) < number_of_cells + 1) error stop "MARS b1 radial size"
        if (size(b2, 1) < number_of_cells + 1) error stop "MARS b2 radial size"
        call require_size("MARS j2 modes", number_of_modes, size(j2, 2))
        call require_size("MARS b1 modes", number_of_modes, size(b1, 2))
        call require_size("MARS b2 modes", number_of_modes, size(b2, 2))
    end subroutine check_mars_profile_shapes

    subroutine require_size(label, expected, actual)
        character(*), intent(in) :: label
        integer, intent(in) :: expected, actual

        if (actual /= expected) error stop label//" has inconsistent size"
    end subroutine require_size

end module response_jxb
