module potato_profile_data_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use phielec_of_psi_mod, only: polydens, polytemp, polyphi
    use spline, only: spline_coeff, spline_val_0
    implicit none
    private
    public :: load_potato_profiles, tabulated_profiles, evaluate_profile_table

    logical :: tabulated_profiles = .false.
    real(dp), allocatable :: coefficients(:, :, :)
    real(dp) :: lower_s, upper_s

contains

    subroutine load_potato_profiles(filename, energy_ev, charge, energy_scale, &
                                   field_scale, e_charge, ev)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: energy_ev, charge, energy_scale, field_scale, e_charge, ev
        character(len=256) :: header
        integer :: unit, ios

        if (energy_ev <= 0 .or. energy_scale <= 0 .or. field_scale == 0) &
            error stop 'invalid POTATO profile normalization'
        tabulated_profiles = .false.
        if (allocated(coefficients)) deallocate(coefficients)
        open(newunit=unit, file=filename, status='old', action='read')
        read(unit, '(A)', iostat=ios) header
        if (ios /= 0) error stop 'cannot read POTATO profile header'
        if (trim(header) == 'POTATO_PROFILE_TABLE_V1') then
            call load_profile_table(unit, energy_ev, charge, energy_scale, &
                                    field_scale, e_charge, ev)
            tabulated_profiles = .true.
        else
            call load_profile_polynomial(unit, energy_ev, charge, energy_scale, &
                                         field_scale, e_charge, ev)
        end if
        close(unit)
    end subroutine load_potato_profiles

    subroutine load_profile_table(unit, energy_ev, charge, energy_scale, &
                                  field_scale, e_charge, ev)
        integer, intent(in) :: unit
        real(dp), intent(in) :: energy_ev, charge, energy_scale, field_scale, e_charge, ev
        integer :: count, k, ios
        real(dp), allocatable :: table(:, :)

        read(unit, *, iostat=ios) count
        if (ios /= 0) error stop 'cannot read POTATO profile row count'
        if (count < 4) error stop 'POTATO profile table needs at least four rows'
        allocate(table(count, 4))
        do k = 1, count
            read(unit, *, iostat=ios) table(k, :)
            if (ios /= 0) error stop 'truncated POTATO profile table'
        end do
        call validate_profile_table(table)
        lower_s = table(1, 1)
        upper_s = table(count, 1)
        table(:, 3) = table(:, 3)/energy_scale/energy_ev
        table(:, 4) = table(:, 4)/field_scale*charge*e_charge/(energy_ev*ev)
        allocate(coefficients(count-1, 5, 3))
        do k = 1, 3
            coefficients(:, :, k) = spline_coeff(table(:, 1), table(:, k+1))
        end do
    end subroutine load_profile_table

    subroutine validate_profile_table(table)
        real(dp), intent(in) :: table(:, :)
        integer :: count

        count = size(table, 1)
        if (.not. all(ieee_is_finite(table))) error stop 'nonfinite POTATO profile table'
        if (any(table(2:, 1) <= table(:count-1, 1))) &
            error stop 'POTATO profile coordinates must increase strictly'
        if (any(table(:, 2:3) <= 0)) &
            error stop 'POTATO density and temperature must be positive'
        if (table(1, 1) > 0 .or. table(count, 1) < 1) &
            error stop 'POTATO profile table must cover s_pol=0..1'
    end subroutine validate_profile_table

    subroutine load_profile_polynomial(unit, energy_ev, charge, energy_scale, &
                                       field_scale, e_charge, ev)
        integer, intent(in) :: unit
        real(dp), intent(in) :: energy_ev, charge, energy_scale, field_scale, e_charge, ev

        ! Preserve the existing six-record polynomial file and normalization.
        read(unit, *)
        read(unit, *) polydens
        read(unit, *) polytemp
        read(unit, *) polytemp
        read(unit, *) polyphi
        polytemp = polytemp/energy_scale/energy_ev
        polyphi = polyphi/field_scale*charge*e_charge/(energy_ev*ev)
    end subroutine load_profile_polynomial

    subroutine evaluate_profile_table(s_pol, column, value, derivative)
        real(dp), intent(in) :: s_pol
        integer, intent(in) :: column
        real(dp), intent(out) :: value, derivative
        real(dp) :: evaluated(3)

        if (.not. tabulated_profiles) error stop 'POTATO tabulated profiles not loaded'
        if (.not. ieee_is_finite(s_pol)) error stop 'nonfinite POTATO profile coordinate'
        ! No silent extrapolation into an undefined edge/SOL profile. Only permit
        ! the roundoff margin in normalizing the dimensional poloidal flux.
        if (s_pol < lower_s-1d-12 .or. s_pol > upper_s+1d-12) &
            error stop 'POTATO orbit is outside tabulated profile support'
        if (column < 1 .or. column > 3) error stop 'invalid POTATO profile column'
        evaluated = spline_val_0(coefficients(:, :, column), s_pol)
        value = evaluated(1)
        derivative = evaluated(2)
        if (column <= 2 .and. value <= 0) &
            error stop 'POTATO interpolated density or temperature is nonpositive'
    end subroutine evaluate_profile_table
end module potato_profile_data_mod
