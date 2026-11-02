program test_resonance_scan_profile
    ! The production profile input is a two-column (s, M_t) table.  The
    ! resonance-scan diagnostic must read one record per surface and preserve
    ! both values; a three-column reader consumes adjacent records incorrectly.
    use iso_fortran_env, only: dp => real64
    use diag_resonance_scan, only: read_profile_table
    implicit none

    real(dp), allocatable :: profile(:, :)
    integer :: unit

    open (newunit=unit, file="profile_two_columns.in", status="replace", action="write")
    write (unit, *) 0.125_dp, 0.25_dp
    write (unit, *) 0.875_dp, -0.50_dp
    close (unit)

    call read_profile_table("profile_two_columns.in", profile)
    if (size(profile, 1) /= 2 .or. size(profile, 2) /= 2) then
        error stop "profile reader did not preserve two rows and two columns"
    end if
    if (abs(profile(1, 1) - 0.125_dp) > 1.0e-14_dp .or. &
        abs(profile(1, 2) - 0.25_dp) > 1.0e-14_dp .or. &
        abs(profile(2, 1) - 0.875_dp) > 1.0e-14_dp .or. &
        abs(profile(2, 2) + 0.50_dp) > 1.0e-14_dp) then
        error stop "profile reader changed a two-column record"
    end if

    print *, "test_resonance_scan_profile PASSED"
end program test_resonance_scan_profile
