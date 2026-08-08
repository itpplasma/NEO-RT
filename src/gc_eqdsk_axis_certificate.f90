module neort_gc_eqdsk_axis_certificate
    !! Interval certificate for one nondegenerate magnetic axis.
    !!
    !! Existence follows from the strict opposing signs of psi_R on the two R
    !! faces and psi_Z on the two Z faces (Poincare-Miranda).  Uniqueness
    !! follows from the Fortsym-generated positive-definite Hessian test on
    !! every EQDSK cell tile (Sylvester criterion).  The same tiles certify a
    !! regular Eq.13 R-chart and positive full along-cut flux curvature.  No
    !! sampled point is promoted to a proof.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, rad, zet
    use neort_gc_eqdsk_cut_interval, only: &
        EQDSK_CUT_INTERVAL_SUCCESS, eqdsk_cut_interval_result_t, &
        evaluate_eqdsk_cut_interval_box
    implicit none
    private

    integer, parameter, public :: EQDSK_AXIS_CERT_SUCCESS = 0
    integer, parameter, public :: EQDSK_AXIS_CERT_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_AXIS_CERT_GRID_INVALID = 2
    integer, parameter, public :: EQDSK_AXIS_CERT_INTERVAL_FAILURE = 3
    integer, parameter, public :: EQDSK_AXIS_CERT_NOT_BRACKETED = 4
    integer, parameter, public :: EQDSK_AXIS_CERT_NOT_CONVEX = 5
    integer, parameter, public :: EQDSK_AXIS_CERT_CUT_DEGENERATE = 6
    integer, parameter, public :: EQDSK_AXIS_CERT_INVALID_CERTIFICATE = 7
    integer, parameter, public :: EQDSK_AXIS_CERTIFICATE_ID = 130014

    type, public :: eqdsk_axis_certificate_t
        real(dp) :: r_lo = 0.0_dp
        real(dp) :: r_hi = 0.0_dp
        real(dp) :: z_lo = 0.0_dp
        real(dp) :: z_hi = 0.0_dp
        real(dp) :: psi_hat_lo = 0.0_dp
        real(dp) :: psi_hat_hi = 0.0_dp
        real(dp) :: psi_RR_lower = 0.0_dp
        real(dp) :: hessian_determinant_lower = 0.0_dp
        real(dp) :: cut_flux_curvature_lower = 0.0_dp
        real(dp) :: left_psi_R_upper = 0.0_dp
        real(dp) :: right_psi_R_lower = 0.0_dp
        real(dp) :: bottom_psi_Z_upper = 0.0_dp
        real(dp) :: top_psi_Z_lower = 0.0_dp
        real(dp) :: left_dpsihat_dR_upper = 0.0_dp
        real(dp) :: right_dpsihat_dR_lower = 0.0_dp
        integer :: cell_tiles_covered = 0
        integer :: certificate_id = 0
        logical :: unique_axis_certified = .false.
        logical :: nondegenerate_cut_limit_certified = .false.
    end type eqdsk_axis_certificate_t

    public :: build_eqdsk_axis_certificate
    public :: validate_eqdsk_axis_certificate

contains

    subroutine build_eqdsk_axis_certificate(r_lo, r_hi, z_lo, z_hi, &
            certificate, status)
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi
        type(eqdsk_axis_certificate_t), intent(out) :: certificate
        integer, intent(out) :: status

        type(eqdsk_cut_interval_result_t) :: tile, face
        integer :: cell_R, cell_Z, interval_status
        real(dp) :: tile_r_lo, tile_r_hi, tile_z_lo, tile_z_hi
        logical :: initialized

        certificate = eqdsk_axis_certificate_t()
        status = EQDSK_AXIS_CERT_INVALID_INPUT
        if (.not. all(ieee_is_finite([r_lo, r_hi, z_lo, z_hi]))) return
        if (r_lo >= r_hi .or. z_lo >= z_hi .or. r_lo <= 0.0_dp) return
        status = EQDSK_AXIS_CERT_GRID_INVALID
        if (.not. valid_grid()) return
        if (r_lo < rad(1) .or. r_hi > rad(nrad) .or. &
                z_lo < zet(1) .or. z_hi > zet(nzet)) return

        certificate%r_lo = r_lo
        certificate%r_hi = r_hi
        certificate%z_lo = z_lo
        certificate%z_hi = z_hi
        certificate%psi_RR_lower = huge(1.0_dp)
        certificate%hessian_determinant_lower = huge(1.0_dp)
        certificate%cut_flux_curvature_lower = huge(1.0_dp)
        certificate%left_psi_R_upper = -huge(1.0_dp)
        certificate%right_psi_R_lower = huge(1.0_dp)
        certificate%bottom_psi_Z_upper = -huge(1.0_dp)
        certificate%top_psi_Z_lower = huge(1.0_dp)
        certificate%left_dpsihat_dR_upper = -huge(1.0_dp)
        certificate%right_dpsihat_dR_lower = huge(1.0_dp)
        initialized = .false.

        do cell_R = 1, nrad-1
            tile_r_lo = max(r_lo, rad(cell_R))
            tile_r_hi = min(r_hi, rad(cell_R+1))
            if (tile_r_hi <= tile_r_lo) cycle
            do cell_Z = 1, nzet-1
                tile_z_lo = max(z_lo, zet(cell_Z))
                tile_z_hi = min(z_hi, zet(cell_Z+1))
                if (tile_z_hi <= tile_z_lo) cycle
                call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
                    tile_r_lo, tile_r_hi, tile_z_lo, tile_z_hi, tile, &
                    interval_status)
                if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
                    status = EQDSK_AXIS_CERT_INTERVAL_FAILURE
                    return
                end if
                if (.not. tile%r_chart_certified) then
                    status = EQDSK_AXIS_CERT_CUT_DEGENERATE
                    return
                end if
                if (tile%psi_RR%lo <= 0.0_dp .or. &
                        tile%axis_hessian_determinant%lo <= 0.0_dp) then
                    status = EQDSK_AXIS_CERT_NOT_CONVEX
                    return
                end if
                if (tile%d2psihat_dR2%lo <= 0.0_dp) then
                    status = EQDSK_AXIS_CERT_CUT_DEGENERATE
                    return
                end if
                call include_tile(certificate, tile, initialized)
                certificate%cell_tiles_covered = &
                    certificate%cell_tiles_covered+1

                if (tile_r_lo == r_lo) then
                    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
                        r_lo, r_lo, tile_z_lo, tile_z_hi, face, interval_status)
                    if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
                        status = EQDSK_AXIS_CERT_INTERVAL_FAILURE
                        return
                    end if
                    if (.not. face%r_chart_certified) then
                        status = EQDSK_AXIS_CERT_CUT_DEGENERATE
                        return
                    end if
                    certificate%left_psi_R_upper = max( &
                        certificate%left_psi_R_upper, face%psi_R%hi)
                    certificate%left_dpsihat_dR_upper = max( &
                        certificate%left_dpsihat_dR_upper, &
                        face%dpsihat_dR%hi)
                end if
                if (tile_r_hi == r_hi) then
                    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
                        r_hi, r_hi, tile_z_lo, tile_z_hi, face, interval_status)
                    if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
                        status = EQDSK_AXIS_CERT_INTERVAL_FAILURE
                        return
                    end if
                    if (.not. face%r_chart_certified) then
                        status = EQDSK_AXIS_CERT_CUT_DEGENERATE
                        return
                    end if
                    certificate%right_psi_R_lower = min( &
                        certificate%right_psi_R_lower, face%psi_R%lo)
                    certificate%right_dpsihat_dR_lower = min( &
                        certificate%right_dpsihat_dR_lower, &
                        face%dpsihat_dR%lo)
                end if
                if (tile_z_lo == z_lo) then
                    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
                        tile_r_lo, tile_r_hi, z_lo, z_lo, face, interval_status)
                    if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
                        status = EQDSK_AXIS_CERT_INTERVAL_FAILURE
                        return
                    end if
                    certificate%bottom_psi_Z_upper = max( &
                        certificate%bottom_psi_Z_upper, face%psi_Z%hi)
                end if
                if (tile_z_hi == z_hi) then
                    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
                        tile_r_lo, tile_r_hi, z_hi, z_hi, face, interval_status)
                    if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
                        status = EQDSK_AXIS_CERT_INTERVAL_FAILURE
                        return
                    end if
                    certificate%top_psi_Z_lower = min( &
                        certificate%top_psi_Z_lower, face%psi_Z%lo)
                end if
            end do
        end do

        if (.not. initialized .or. certificate%cell_tiles_covered < 1) then
            status = EQDSK_AXIS_CERT_GRID_INVALID
            return
        end if
        if (certificate%left_psi_R_upper >= 0.0_dp .or. &
                certificate%right_psi_R_lower <= 0.0_dp .or. &
                certificate%bottom_psi_Z_upper >= 0.0_dp .or. &
                certificate%top_psi_Z_lower <= 0.0_dp .or. &
                certificate%left_dpsihat_dR_upper >= 0.0_dp .or. &
                certificate%right_dpsihat_dR_lower <= 0.0_dp) then
            status = EQDSK_AXIS_CERT_NOT_BRACKETED
            return
        end if
        certificate%unique_axis_certified = .true.
        certificate%nondegenerate_cut_limit_certified = .true.
        certificate%certificate_id = EQDSK_AXIS_CERTIFICATE_ID
        call validate_eqdsk_axis_certificate(certificate, status)
    end subroutine build_eqdsk_axis_certificate

    subroutine validate_eqdsk_axis_certificate(certificate, status)
        type(eqdsk_axis_certificate_t), intent(in) :: certificate
        integer, intent(out) :: status

        status = EQDSK_AXIS_CERT_INVALID_CERTIFICATE
        if (.not. certificate%unique_axis_certified) return
        if (.not. certificate%nondegenerate_cut_limit_certified) return
        if (certificate%certificate_id /= EQDSK_AXIS_CERTIFICATE_ID) return
        if (certificate%cell_tiles_covered < 1) return
        if (.not. all(ieee_is_finite([certificate%r_lo, certificate%r_hi, &
                certificate%z_lo, certificate%z_hi, certificate%psi_hat_lo, &
                certificate%psi_hat_hi, certificate%psi_RR_lower, &
                certificate%hessian_determinant_lower, &
                certificate%cut_flux_curvature_lower, &
                certificate%left_psi_R_upper, &
                certificate%right_psi_R_lower, &
                certificate%bottom_psi_Z_upper, &
                certificate%top_psi_Z_lower, &
                certificate%left_dpsihat_dR_upper, &
                certificate%right_dpsihat_dR_lower]))) return
        if (certificate%r_lo >= certificate%r_hi .or. &
                certificate%z_lo >= certificate%z_hi) return
        if (certificate%psi_hat_lo > certificate%psi_hat_hi) return
        if (certificate%psi_RR_lower <= 0.0_dp .or. &
                certificate%hessian_determinant_lower <= 0.0_dp .or. &
                certificate%cut_flux_curvature_lower <= 0.0_dp) return
        if (certificate%left_psi_R_upper >= 0.0_dp .or. &
                certificate%right_psi_R_lower <= 0.0_dp .or. &
                certificate%bottom_psi_Z_upper >= 0.0_dp .or. &
                certificate%top_psi_Z_lower <= 0.0_dp .or. &
                certificate%left_dpsihat_dR_upper >= 0.0_dp .or. &
                certificate%right_dpsihat_dR_lower <= 0.0_dp) return
        status = EQDSK_AXIS_CERT_SUCCESS
    end subroutine validate_eqdsk_axis_certificate

    subroutine include_tile(certificate, tile, initialized)
        type(eqdsk_axis_certificate_t), intent(inout) :: certificate
        type(eqdsk_cut_interval_result_t), intent(in) :: tile
        logical, intent(inout) :: initialized

        if (.not. initialized) then
            certificate%psi_hat_lo = tile%psi_hat%lo
            certificate%psi_hat_hi = tile%psi_hat%hi
            initialized = .true.
        else
            certificate%psi_hat_lo = min(certificate%psi_hat_lo, &
                tile%psi_hat%lo)
            certificate%psi_hat_hi = max(certificate%psi_hat_hi, &
                tile%psi_hat%hi)
        end if
        certificate%psi_RR_lower = min(certificate%psi_RR_lower, &
            tile%psi_RR%lo)
        certificate%hessian_determinant_lower = min( &
            certificate%hessian_determinant_lower, &
            tile%axis_hessian_determinant%lo)
        certificate%cut_flux_curvature_lower = min( &
            certificate%cut_flux_curvature_lower, &
            tile%d2psihat_dR2%lo)
    end subroutine include_tile

    logical function valid_grid()
        integer :: i

        valid_grid = .false.
        if (nrad < 2 .or. nzet < 2) return
        if (.not. allocated(rad) .or. .not. allocated(zet)) return
        if (size(rad) /= nrad .or. size(zet) /= nzet) return
        if (.not. all(ieee_is_finite(rad)) .or. &
                .not. all(ieee_is_finite(zet))) return
        do i = 1, nrad-1
            if (rad(i+1) <= rad(i)) return
        end do
        do i = 1, nzet-1
            if (zet(i+1) <= zet(i)) return
        end do
        valid_grid = .true.
    end function valid_grid

end module neort_gc_eqdsk_axis_certificate
