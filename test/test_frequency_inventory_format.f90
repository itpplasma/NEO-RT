program test_frequency_inventory_format
    !! Independent manufactured oracle for the whitespace-delimited row
    !! contracts.  It does not import the diagnostic implementation.
    use iso_fortran_env, only: dp => real64
    implicit none

    character(len=8) :: class_label
    integer :: model, sign_value, status, orbit_status, mph, mth
    real(dp) :: eta, etatp, ratio, omega_b, omega_phi, omega_magnetic
    real(dp) :: omega_electric, residual, derivative
    character(len=1024) :: line

    line = '0 passing -1 1.0 2.0 0.5 3.0 4.0 NaN NaN 7 9'
    read (line, *) model, class_label, sign_value, eta, etatp, ratio, omega_b, &
        omega_phi, omega_magnetic, omega_electric, status, orbit_status
    if (model /= 0 .or. trim(class_label) /= 'passing' .or. sign_value /= -1) then
        error stop 'frequency row prefix format failure'
    end if
    if (abs(eta - 1.0_dp) > 1.0e-12_dp .or. abs(ratio - 0.5_dp) > 1.0e-12_dp) then
        error stop 'frequency row numeric format failure'
    end if
    if (status /= 7 .or. orbit_status /= 9) error stop 'status column format failure'

    line = '2 trapped 1 3 -2 1.0 2.0 0.5 0.0 NaN 3 4'
    read (line, *) model, class_label, sign_value, mph, mth, eta, etatp, ratio, &
        residual, derivative, status, orbit_status
    if (model /= 2 .or. trim(class_label) /= 'trapped' .or. sign_value /= 1) then
        error stop 'resonance row prefix format failure'
    end if
    if (mph /= 3 .or. mth /= -2 .or. abs(ratio - 0.5_dp) > 1.0e-12_dp) then
        error stop 'resonance row numeric format failure'
    end if
    if (status /= 3 .or. orbit_status /= 4) then
        error stop 'resonance status column format failure'
    end if
end program test_frequency_inventory_format
