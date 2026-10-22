program reader_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32
    use potato_input_mod, only: read_potato_input
    use do_magfie_pert_mod, only: do_magfie_pert_amp, init_magfie_pert_at_s
    implicit none
    character(len=512) :: path
    integer :: input, output, j
    integer(int32) :: n
    real(dp), allocatable :: query(:, :)
    complex(dp), allocatable :: got(:), expected(:)
    complex(dp) :: amp
    real(dp) :: x(3)

    call read_potato_input('potato.in')
    call get_command_argument(1, path)
    open(newunit=input, file=trim(path), form='unformatted', access='stream', &
         status='old')
    read(input) n
    allocate(query(n, 6), got(n), expected(n))
    ! Columns: R[cm], Z[cm], s, theta_B, Re(phase), Im(phase).
    read(input) query
    close(input)
    call bmod_pert(query(1, 1), query(1, 2), got(1))
    !$omp parallel private(x,amp)
    call init_magfie_pert_at_s()
    !$omp do
    do j = 1, n
        call bmod_pert(query(j, 1), query(j, 2), got(j))
        x = [query(j, 3), 0.0_dp, query(j, 4)]
        call do_magfie_pert_amp(x, amp)
        expected(j) = amp*cmplx(query(j, 5), query(j, 6), dp)
    end do
    !$omp end do
    !$omp end parallel
    call get_command_argument(2, path)
    open(newunit=output, file=trim(path), form='unformatted', access='stream', &
         status='replace')
    write(output) got, expected
    close(output)
end program reader_readback
