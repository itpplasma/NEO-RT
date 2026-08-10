program test_neo2_reader
    use neo2_reader, only: dp, neo2_data_t, read_neo2_data
    implicit none

    type(neo2_data_t) :: data
    character(len=1024) :: fixture
    integer :: i

    call get_command_argument(1, fixture)
    call read_neo2_data(trim(fixture), data)

    if (any(abs(data%stor - [0.0_dp, 0.5_dp, 1.0_dp]) > 1.0e-12_dp)) &
        error stop "toroidal flux differs"
    if (any(abs(data%q - 2.0_dp) > 1.0e-12_dp)) error stop "safety factor differs"
    if (any(abs(data%spol - [0.0_dp, 0.5_dp, 1.0_dp]) > 1.0e-12_dp)) &
        error stop "poloidal flux differs"
    if (any(shape(data%n_spec) /= [3, 2])) error stop "rank-two shape differs"
    if (any(abs(data%n_spec - reshape([(real(i, dp), i=1, 6)], [3, 2])) > &
        1.0e-12_dp)) error stop "rank-two values differ"
    if (any(abs(data%Er - [-2.0_dp, 0.0_dp, 2.0_dp]) > 1.0e-12_dp)) &
        error stop "electric field differs"
end program test_neo2_reader
