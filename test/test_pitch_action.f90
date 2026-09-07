program test_pitch_action
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use diag_pitch_action, only: pitch_point_t, read_pitch_points, &
        resonance_coefficients
    implicit none
    type(pitch_point_t), allocatable :: points(:)
    real(dp) :: coeff(3)
    character(len=32) :: mode
    character(len=128) :: path
    integer :: unit

    call get_command_argument(1, mode)
    if (len_trim(mode) == 0) mode = "valid"
    path = "pitch_points_"//trim(mode)//".txt"
    open (newunit=unit, file=trim(path), status="replace")
    select case (trim(mode))
    case ("valid")
        write (unit, '(A)') "# s branch mth eta ux"
        write (unit, '(A)') "0.5 3 -1 2.0D-5 1.0"
        write (unit, '(A)') ""
        write (unit, '(A)') "0.5 2 2 1.0D-5 2.3"
    case ("extra")
        write (unit, '(A)') "0.5 3 -1 2.0D-5 1.0 extra"
    case ("missing")
        write (unit, '(A)') "0.5 3 -1 2.0D-5"
    case ("branch")
        write (unit, '(A)') "0.5 4 -1 2.0D-5 1.0"
    case ("nonfinite")
        write (unit, '(A)') "0.5 3 -1 NaN 1.0"
    case ("endpoint")
        write (unit, '(A)') "1.0 3 -1 2.0D-5 1.0"
    case ("empty")
        write (unit, '(A)') "# no data"
    case default
        error stop "invalid test mode"
    end select
    close (unit)
    call read_pitch_points(trim(path), points)
    if (mode /= "valid") stop 0 ! WILL_FAIL tests must reject in the reader.
    if (size(points) /= 2) error stop "wrong point count"
    if (points(1)%harmonic /= -1) error stop "signed harmonic changed"
    if (points(2)%branch /= 2) error stop "branch changed"
    if (abs(points(2)%ux - 2.3_dp) > 1.0e-14_dp) error stop "speed changed"
    if (abs(points(1)%eta - 2.0e-5_dp) > 1.0e-18_dp) error stop "pitch changed"
    coeff = resonance_coefficients(-1, 3, .false., -0.5_dp, -7.0_dp, &
        5.0_dp, -2.0_dp)
    if (coeff(1) /= -6.0_dp) error stop "quadratic coefficient"
    if (coeff(2) /= -5.0_dp) error stop "trapped linear coefficient"
    if (coeff(3) /= -21.0_dp) error stop "electric coefficient"
    coeff = resonance_coefficients(-1, 3, .true., -0.5_dp, -7.0_dp, &
        5.0_dp, -2.0_dp)
    if (coeff(2) /= -35.0_dp) error stop "passing transit shift"
    coeff = resonance_coefficients(-1, 3, .true., -0.5_dp, -7.0_dp, &
        -5.0_dp, -2.0_dp)
    if (coeff(2) /= 35.0_dp) error stop "counter-passing orientation"
end program test_pitch_action
