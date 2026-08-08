module neort_gc_eqdsk_cut_interval
    !! Outward enclosure of the Fortsym-generated EQDSK Eq. 13 numerator.
    !!
    !! The caller owns spatial subdivision and supplies one tensor-quintic
    !! cell.  This module marshals libneo coefficients into generated interval
    !! kernels and takes a hull over every F(psi) spline branch reached by the
    !! box.  It contains no handwritten Eq. 13 or field algebra.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: btf, hfpol, icp, ipoint, nrad, nzet, rad, zet, &
        psi_sep, rtf, splfpol, splpsi, use_fpol
    use neort_eqdsk_cut_numerator_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_numerator_interval
    use neort_eqdsk_cut_mean_value_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_mean_value_interval
    use neort_eqdsk_cut_axis_curvature_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_curvature_interval
    use neort_eqdsk_cut_r_flux_chart_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_chart_interval
    use neort_eqdsk_quintic_cell_jet_interval_symbolic, only: &
        evaluate_neort_eqdsk_quintic_cell_jet_interval
    use neort_eqdsk_quintic_profile_jet_interval_symbolic, only: &
        evaluate_neort_eqdsk_quintic_profile_jet_interval
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t, &
        operator(-), operator(/)
    implicit none
    private

    integer, parameter, public :: EQDSK_CUT_INTERVAL_SUCCESS = 0
    integer, parameter, public :: EQDSK_CUT_INTERVAL_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_CUT_INTERVAL_UNINITIALIZED = 2
    integer, parameter, public :: EQDSK_CUT_INTERVAL_CELL_MISMATCH = 3
    integer, parameter, public :: EQDSK_CUT_INTERVAL_NONFINITE = 4

    type, public :: eqdsk_cut_interval_result_t
        integer :: cell_R = 0
        integer :: cell_Z = 0
        type(gc_outward_interval_t) :: psi
        type(gc_outward_interval_t) :: psi_hat
        type(gc_outward_interval_t) :: F
        type(gc_outward_interval_t) :: dF_dpsi_hat
        type(gc_outward_interval_t) :: numerator
        type(gc_outward_interval_t) :: numerator_R
        type(gc_outward_interval_t) :: numerator_Z
        type(gc_outward_interval_t) :: positive_denominator
        type(gc_outward_interval_t) :: dZ_dR
        type(gc_outward_interval_t) :: dpsihat_dR
        type(gc_outward_interval_t) :: axis_flux_curvature
        type(gc_outward_interval_t) :: axis_hessian_determinant
        logical :: denominator_positive_certified = .false.
        logical :: r_chart_certified = .false.
        integer :: profile_cells_covered = 0
        logical :: vacuum_branch_covered = .false.
    end type eqdsk_cut_interval_result_t

    public :: evaluate_eqdsk_cut_interval_box

contains

    subroutine evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, R_lo, R_hi, &
            Z_lo, Z_hi, result, status)
        integer, intent(in) :: cell_R, cell_Z
        real(dp), intent(in) :: R_lo, R_hi, Z_lo, Z_hi
        type(eqdsk_cut_interval_result_t), intent(out) :: result
        integer, intent(out) :: status

        type(gc_outward_interval_t) :: coefficient(6,6), jet(10)
        type(gc_outward_interval_t) :: midpoint_jet(10)
        type(gc_outward_interval_t) :: radius, delta_R, delta_Z, separatrix
        type(gc_outward_interval_t) :: midpoint_profile(4), profile(4)
        type(gc_outward_interval_t) :: midpoint_numerator
        type(gc_outward_interval_t) :: midpoint_numerator_R
        type(gc_outward_interval_t) :: midpoint_numerator_Z
        type(gc_outward_interval_t) :: midpoint_denominator
        type(gc_outward_interval_t) :: mean_value_numerator
        type(gc_outward_interval_t) :: centered_R, centered_Z
        real(dp) :: R_mid, Z_mid
        integer :: cell_pointer, i, j, midpoint_cells, midpoint_status
        logical :: midpoint_vacuum

        result = eqdsk_cut_interval_result_t()
        status = EQDSK_CUT_INTERVAL_INVALID_INPUT
        if (.not. all(ieee_is_finite([R_lo, R_hi, Z_lo, Z_hi]))) return
        if (R_lo <= 0.0_dp .or. R_hi < R_lo .or. Z_hi < Z_lo) return
        status = EQDSK_CUT_INTERVAL_UNINITIALIZED
        if (nrad < 2 .or. nzet < 2 .or. icp < 1) return
        if (.not. allocated(rad) .or. .not. allocated(zet)) return
        if (.not. allocated(splpsi) .or. .not. allocated(ipoint)) return
        if (size(rad) /= nrad .or. size(zet) /= nzet) return
        if (size(splpsi,1) /= 6 .or. size(splpsi,2) /= 6) return
        if (size(splpsi,3) < icp) return
        if (size(ipoint,1) /= nrad .or. size(ipoint,2) /= nzet) return
        if (.not. all(ieee_is_finite([rad, zet, psi_sep, btf, rtf]))) return
        ! Profile-cell ordering assumes psi_hat=psi/psi_sep increases from
        ! zero to one.  Signed-flux support would require a distinct generated
        ! normalization contract rather than silently accepting reversal.
        if (psi_sep <= tiny(1.0_dp)) return

        status = EQDSK_CUT_INTERVAL_CELL_MISMATCH
        if (cell_R < 1 .or. cell_R >= nrad) return
        if (cell_Z < 1 .or. cell_Z >= nzet) return
        if (R_lo < rad(cell_R) .or. R_hi > rad(cell_R+1)) return
        if (Z_lo < zet(cell_Z) .or. Z_hi > zet(cell_Z+1)) return
        cell_pointer = ipoint(cell_R, cell_Z)
        if (cell_pointer < 1 .or. cell_pointer > size(splpsi,3)) return
        if (.not. all(ieee_is_finite(splpsi(:,:,cell_pointer)))) then
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if

        do i = 1, 6
            do j = 1, 6
                coefficient(i,j) = point(splpsi(i,j,cell_pointer))
            end do
        end do
        radius = gc_outward_interval(R_lo, R_hi)
        delta_R = radius-point(rad(cell_R))
        delta_Z = gc_outward_interval(Z_lo, Z_hi)-point(zet(cell_Z))
        call evaluate_cell_jet(delta_R, delta_Z, coefficient, jet)
        if (.not. valid_intervals(jet)) then
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if

        result%cell_R = cell_R
        result%cell_Z = cell_Z
        result%psi = jet(1)
        separatrix = point(psi_sep)
        result%psi_hat = result%psi/separatrix
        call evaluate_profile_hull(result%psi_hat, profile, &
            result%profile_cells_covered, result%vacuum_branch_covered, status)
        if (status /= EQDSK_CUT_INTERVAL_SUCCESS) return
        result%F = profile(1)
        result%dF_dpsi_hat = profile(2)

        call evaluate_neort_eqdsk_cut_numerator_interval(radius, jet(2), &
            jet(3), jet(4), jet(5), jet(6), jet(7), jet(8), jet(9), jet(10), &
            profile(1), profile(2), separatrix, result%numerator, &
            result%numerator_R, result%numerator_Z, &
            result%positive_denominator)
        if (.not. valid_intervals([result%numerator, result%numerator_R, &
                result%numerator_Z, result%positive_denominator])) then
            result = eqdsk_cut_interval_result_t()
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if

        ! Tighten the direct natural interval with the Fortsym-generated
        ! mean-value enclosure.  The midpoint value and the full-box partial
        ! derivative enclosures are independently outward-rounded; their
        ! intersection therefore remains a certified enclosure of N(R,Z).
        R_mid = R_lo+0.5_dp*(R_hi-R_lo)
        Z_mid = Z_lo+0.5_dp*(Z_hi-Z_lo)
        call evaluate_cell_jet(point(R_mid-rad(cell_R)), &
            point(Z_mid-zet(cell_Z)), coefficient, midpoint_jet)
        if (.not. valid_intervals(midpoint_jet)) then
            result = eqdsk_cut_interval_result_t()
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if
        call evaluate_profile_hull(midpoint_jet(1)/separatrix, &
            midpoint_profile, midpoint_cells, midpoint_vacuum, midpoint_status)
        if (midpoint_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
            result = eqdsk_cut_interval_result_t()
            status = midpoint_status
            return
        end if
        call evaluate_neort_eqdsk_cut_numerator_interval(point(R_mid), &
            midpoint_jet(2), midpoint_jet(3), midpoint_jet(4), &
            midpoint_jet(5), midpoint_jet(6), midpoint_jet(7), &
            midpoint_jet(8), midpoint_jet(9), midpoint_jet(10), &
            midpoint_profile(1), midpoint_profile(2), separatrix, &
            midpoint_numerator, midpoint_numerator_R, midpoint_numerator_Z, &
            midpoint_denominator)
        centered_R = radius-point(R_mid)
        centered_Z = gc_outward_interval(Z_lo, Z_hi)-point(Z_mid)
        call evaluate_neort_eqdsk_cut_mean_value_interval(midpoint_numerator, &
            result%numerator_R, result%numerator_Z, centered_R, centered_Z, &
            mean_value_numerator)
        if (.not. valid_intervals([midpoint_numerator, &
                midpoint_numerator_R, midpoint_numerator_Z, &
                midpoint_denominator, mean_value_numerator])) then
            result = eqdsk_cut_interval_result_t()
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if
        result%numerator%lo = max(result%numerator%lo, mean_value_numerator%lo)
        result%numerator%hi = min(result%numerator%hi, mean_value_numerator%hi)
        if (.not. gc_outward_interval_is_valid(result%numerator)) then
            result = eqdsk_cut_interval_result_t()
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if
        result%denominator_positive_certified = &
            result%positive_denominator%lo > 0.0_dp
        if (result%numerator_Z%lo > 0.0_dp .or. &
                result%numerator_Z%hi < 0.0_dp) then
            call evaluate_neort_eqdsk_cut_r_flux_chart_interval( &
                result%numerator_R, result%numerator_Z, jet(2), jet(3), &
                separatrix, result%dZ_dR, result%dpsihat_dR)
            if (.not. valid_intervals([result%dZ_dR, &
                    result%dpsihat_dR])) then
                result = eqdsk_cut_interval_result_t()
                status = EQDSK_CUT_INTERVAL_NONFINITE
                return
            end if
            call evaluate_neort_eqdsk_cut_axis_curvature_interval( &
                result%dZ_dR, jet(4), jet(5), jet(6), separatrix, &
                result%axis_flux_curvature, &
                result%axis_hessian_determinant)
            if (.not. valid_intervals([result%axis_flux_curvature, &
                    result%axis_hessian_determinant])) then
                result = eqdsk_cut_interval_result_t()
                status = EQDSK_CUT_INTERVAL_NONFINITE
                return
            end if
            result%r_chart_certified = .true.
        end if
        status = EQDSK_CUT_INTERVAL_SUCCESS
    end subroutine evaluate_eqdsk_cut_interval_box

    subroutine evaluate_profile_hull(psi_hat, profile, cells_covered, &
            vacuum_covered, status)
        type(gc_outward_interval_t), intent(in) :: psi_hat
        type(gc_outward_interval_t), intent(out) :: profile(4)
        integer, intent(out) :: cells_covered
        logical, intent(out) :: vacuum_covered
        integer, intent(out) :: status

        type(gc_outward_interval_t) :: coefficient(0:5), branch(4)
        type(gc_outward_interval_t) :: delta, vacuum_b, vacuum_r
        real(dp) :: branch_lo, branch_hi, overlap_lo, overlap_hi, origin
        integer :: cell, i
        logical :: initialized

        profile = point(0.0_dp)
        cells_covered = 0
        vacuum_covered = .false.
        status = EQDSK_CUT_INTERVAL_UNINITIALIZED
        if (.not. gc_outward_interval_is_valid(psi_hat)) then
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if
        vacuum_b = point(btf)
        vacuum_r = point(rtf)
        if (.not. use_fpol) then
            coefficient = point(0.0_dp)
            delta = point(0.0_dp)
            call evaluate_profile_jet(delta, coefficient, vacuum_b, vacuum_r, &
                profile)
            vacuum_covered = .true.
            status = EQDSK_CUT_INTERVAL_SUCCESS
            return
        end if
        if (.not. allocated(splfpol)) return
        if (size(splfpol,1) < 6 .or. size(splfpol,2) < nrad) return
        if (.not. ieee_is_finite(hfpol) .or. hfpol <= 0.0_dp) return
        if (.not. all(ieee_is_finite(splfpol(0:5,1:nrad)))) then
            status = EQDSK_CUT_INTERVAL_NONFINITE
            return
        end if

        initialized = .false.
        do cell = 1, nrad
            origin = real(cell-1,dp)*hfpol
            if (cell == 1) then
                branch_lo = min(psi_hat%lo, 0.0_dp)
            else
                branch_lo = origin
            end if
            branch_hi = origin+hfpol
            if (cell == nrad) branch_hi = 1.0_dp
            overlap_lo = max(psi_hat%lo, branch_lo)
            overlap_hi = min(psi_hat%hi, branch_hi)
            if (overlap_hi < overlap_lo) cycle
            do i = 0, 5
                coefficient(i) = point(splfpol(i,cell))
            end do
            delta = gc_outward_interval(overlap_lo, overlap_hi)-point(origin)
            call evaluate_profile_jet(delta, coefficient, vacuum_b, vacuum_r, &
                branch)
            if (.not. valid_intervals(branch)) then
                status = EQDSK_CUT_INTERVAL_NONFINITE
                return
            end if
            call include_profile_branch(profile, branch, initialized)
            cells_covered = cells_covered+1
        end do

        if (psi_hat%hi > 1.0_dp) then
            branch(1) = point(splfpol(0,nrad))
            branch(2:4) = point(0.0_dp)
            call include_profile_branch(profile, branch, initialized)
            vacuum_covered = .true.
        end if
        if (.not. initialized) return
        status = EQDSK_CUT_INTERVAL_SUCCESS
    end subroutine evaluate_profile_hull

    subroutine include_profile_branch(profile, branch, initialized)
        type(gc_outward_interval_t), intent(inout) :: profile(4)
        type(gc_outward_interval_t), intent(in) :: branch(4)
        logical, intent(inout) :: initialized
        integer :: i

        if (.not. initialized) then
            profile = branch
            initialized = .true.
            return
        end if
        do i = 1, 4
            profile(i)%lo = min(profile(i)%lo, branch(i)%lo)
            profile(i)%hi = max(profile(i)%hi, branch(i)%hi)
        end do
    end subroutine include_profile_branch

    subroutine evaluate_cell_jet(delta_R, delta_Z, coefficient, jet)
        type(gc_outward_interval_t), intent(in) :: delta_R, delta_Z
        type(gc_outward_interval_t), intent(in) :: coefficient(6,6)
        type(gc_outward_interval_t), intent(out) :: jet(10)

        call evaluate_neort_eqdsk_quintic_cell_jet_interval(delta_R, delta_Z, &
            coefficient(1,1), coefficient(1,2), coefficient(1,3), &
            coefficient(1,4), coefficient(1,5), coefficient(1,6), &
            coefficient(2,1), coefficient(2,2), coefficient(2,3), &
            coefficient(2,4), coefficient(2,5), coefficient(2,6), &
            coefficient(3,1), coefficient(3,2), coefficient(3,3), &
            coefficient(3,4), coefficient(3,5), coefficient(3,6), &
            coefficient(4,1), coefficient(4,2), coefficient(4,3), &
            coefficient(4,4), coefficient(4,5), coefficient(4,6), &
            coefficient(5,1), coefficient(5,2), coefficient(5,3), &
            coefficient(5,4), coefficient(5,5), coefficient(5,6), &
            coefficient(6,1), coefficient(6,2), coefficient(6,3), &
            coefficient(6,4), coefficient(6,5), coefficient(6,6), &
            jet(1), jet(2), jet(3), jet(4), jet(5), jet(6), jet(7), jet(8), &
            jet(9), jet(10))
    end subroutine evaluate_cell_jet

    subroutine evaluate_profile_jet(delta, coefficient, vacuum_b, vacuum_r, &
            profile)
        type(gc_outward_interval_t), intent(in) :: delta, coefficient(0:5)
        type(gc_outward_interval_t), intent(in) :: vacuum_b, vacuum_r
        type(gc_outward_interval_t), intent(out) :: profile(4)

        call evaluate_neort_eqdsk_quintic_profile_jet_interval(delta, &
            coefficient(0), coefficient(1), coefficient(2), coefficient(3), &
            coefficient(4), coefficient(5), vacuum_b, vacuum_r, profile(1), &
            profile(2), profile(3), profile(4))
    end subroutine evaluate_profile_jet

    pure elemental function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    logical function valid_intervals(values)
        type(gc_outward_interval_t), intent(in) :: values(:)
        integer :: i

        valid_intervals = .true.
        do i = 1, size(values)
            valid_intervals = valid_intervals .and. &
                gc_outward_interval_is_valid(values(i))
        end do
    end function valid_intervals

end module neort_gc_eqdsk_cut_interval
