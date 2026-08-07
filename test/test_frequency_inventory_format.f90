program test_frequency_inventory_format
    !! Independent manufactured oracle for the public labels and the
    !! whitespace-delimited frequency-row column contract.
    use iso_fortran_env, only: dp => real64
    use diag_frequency_inventory, only: frequency_inventory_model_label, &
        frequency_inventory_class_label
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_TRAPPED
    implicit none

    character(len=8) :: label
    integer :: model, sign_value, status, orbit_status
    real(dp) :: eta, etatp, ratio, omega_b, omega_phi, omega_magnetic
    real(dp) :: omega_electric
    character(len=1024) :: line

    label = frequency_inventory_model_label(0)
    if (trim(label) /= 'boozer0') error stop 'model 0 label format failure'
    label = frequency_inventory_model_label(2)
    if (trim(label) /= 'real2') error stop 'model 2 label format failure'
    label = frequency_inventory_class_label(GC_ORBIT_PASSING)
    if (trim(label) /= 'passing') error stop 'passing label format failure'
    label = frequency_inventory_class_label(GC_ORBIT_TRAPPED)
    if (trim(label) /= 'trapped') error stop 'trapped label format failure'

    line = '2 passing -1 1.0 2.0 0.5 3.0 4.0 NaN NaN 7 9'
    read (line, *) model, label, sign_value, eta, etatp, ratio, omega_b, &
        omega_phi, omega_magnetic, omega_electric, status, orbit_status
    if (model /= 2 .or. trim(label) /= 'passing' .or. sign_value /= -1) then
        error stop 'frequency row prefix format failure'
    end if
    if (abs(eta - 1.0_dp) > 1.0e-12_dp .or. abs(ratio - 0.5_dp) > 1.0e-12_dp) then
        error stop 'frequency row numeric format failure'
    end if
    if (status /= 7 .or. orbit_status /= 9) error stop 'status column format failure'
end program test_frequency_inventory_format
