program test_mars_response
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mars_response_io, only: mars_vector_harmonics_t, read_mars_vector
    use response_jxb, only: reconstruct_mars_profile

    implicit none

    character(1024) :: fixture_dir
    type(mars_vector_harmonics_t) :: magnetic, current
    real(dp), allocatable :: raw(:), processed(:)
    real(dp) :: half_mesh(3)

    call get_command_argument(1, fixture_dir)
    call read_mars_vector(trim(fixture_dir)//"/mars_bplasma_small.out", magnetic)
    call read_mars_vector(trim(fixture_dir)//"/mars_jplasma_small.out", current)

    if (magnetic%radial_count /= 4) error stop "magnetic radial count"
    if (magnetic%mode_count /= 1) error stop "magnetic mode count"
    if (magnetic%toroidal_mode /= -3) error stop "magnetic toroidal mode"
    if (magnetic%modes(1) /= 2) error stop "magnetic poloidal mode"
    if (current%modes(1) /= magnetic%modes(1)) error stop "mode mismatch"

    half_mesh = [0.2_dp, 0.5_dp, 0.8_dp]
    call reconstruct_mars_profile( &
        current%c1, current%c2, magnetic%c1, magnetic%c2, half_mesh, &
        first_torque_cell=2, smoothing_passes=0, edge_cutoff=1.0_dp, &
        raw=raw, processed=processed)

    call assert_vector(raw, [0.0_dp, 2.0_dp, 3.0_dp], "raw profile")
    call assert_vector(processed, raw, "processed profile")

contains

    subroutine assert_vector(actual, expected, label)
        real(dp), intent(in) :: actual(:), expected(:)
        character(*), intent(in) :: label

        if (size(actual) /= size(expected)) error stop label//" size"
        if (maxval(abs(actual - expected)) > 1.0e-13_dp) error stop label
    end subroutine assert_vector

end program test_mars_response
