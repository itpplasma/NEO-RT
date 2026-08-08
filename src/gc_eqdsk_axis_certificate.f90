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
    use neort_eqdsk_axis_stationarity_krawczyk_interval_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_krawczyk_interval
    use neort_eqdsk_axis_stationarity_newton_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_newton
    use neort_eqdsk_axis_stationarity_system_interval_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_system_interval
    use neort_eqdsk_axis_stationarity_system_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_system
    use neort_gc_eqdsk_cut_interval, only: &
        EQDSK_CUT_INTERVAL_SUCCESS, eqdsk_cut_interval_result_t, &
        evaluate_eqdsk_cut_interval_box
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t, operator(-)
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
    integer, parameter, public :: EQDSK_AXIS_CERT_POINT_FAILURE = 8
    integer, parameter, public :: EQDSK_AXIS_CERT_KRAWCZYK_FAILURE = 9
    integer, parameter, public :: EQDSK_AXIS_CERTIFICATE_ID = 130014

    type, public :: eqdsk_axis_certificate_t
        real(dp) :: r_lo = 0.0_dp
        real(dp) :: r_hi = 0.0_dp
        real(dp) :: z_lo = 0.0_dp
        real(dp) :: z_hi = 0.0_dp
        real(dp) :: psi_hat_lo = 0.0_dp
        real(dp) :: psi_hat_hi = 0.0_dp
        real(dp) :: psi_R_lower = 0.0_dp
        real(dp) :: psi_R_upper = 0.0_dp
        real(dp) :: psi_Z_lower = 0.0_dp
        real(dp) :: psi_Z_upper = 0.0_dp
        real(dp) :: psi_RR_lower = 0.0_dp
        real(dp) :: psi_RR_upper = 0.0_dp
        real(dp) :: psi_RZ_lower = 0.0_dp
        real(dp) :: psi_RZ_upper = 0.0_dp
        real(dp) :: psi_ZZ_lower = 0.0_dp
        real(dp) :: psi_ZZ_upper = 0.0_dp
        real(dp) :: hessian_determinant_lower = 0.0_dp
        real(dp) :: cut_flux_curvature_lower = 0.0_dp
        real(dp) :: left_psi_R_upper = 0.0_dp
        real(dp) :: right_psi_R_lower = 0.0_dp
        real(dp) :: bottom_psi_Z_upper = 0.0_dp
        real(dp) :: top_psi_Z_lower = 0.0_dp
        real(dp) :: left_dpsihat_dR_upper = 0.0_dp
        real(dp) :: right_dpsihat_dR_lower = 0.0_dp
        real(dp) :: axis_newton_R = 0.0_dp
        real(dp) :: axis_newton_Z = 0.0_dp
        real(dp) :: axis_R_lo = 0.0_dp
        real(dp) :: axis_R_hi = 0.0_dp
        real(dp) :: axis_Z_lo = 0.0_dp
        real(dp) :: axis_Z_hi = 0.0_dp
        real(dp) :: axis_hessian_determinant = 0.0_dp
        integer :: cell_tiles_covered = 0
        integer :: certificate_id = 0
        logical :: unique_axis_certified = .false.
        logical :: nondegenerate_cut_limit_certified = .false.
        logical :: stationarity_enclosure_certified = .false.
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
        type(gc_outward_interval_t) :: interval_system(6)
        type(gc_outward_interval_t) :: delta_R, delta_Z
        type(gc_outward_interval_t) :: krawczyk_R, krawczyk_Z
        real(dp) :: point_system(6), point_R, point_Z, next_R, next_Z
        real(dp) :: determinant, inverse(4)
        integer :: cell_R, cell_Z, interval_status
        integer :: iteration, point_status
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
        certificate%psi_RR_upper = -huge(1.0_dp)
        certificate%psi_RZ_lower = huge(1.0_dp)
        certificate%psi_RZ_upper = -huge(1.0_dp)
        certificate%psi_ZZ_lower = huge(1.0_dp)
        certificate%psi_ZZ_upper = -huge(1.0_dp)
        certificate%psi_R_lower = huge(1.0_dp)
        certificate%psi_R_upper = -huge(1.0_dp)
        certificate%psi_Z_lower = huge(1.0_dp)
        certificate%psi_Z_upper = -huge(1.0_dp)
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

        call evaluate_neort_eqdsk_axis_stationarity_system_interval( &
            gc_outward_interval(certificate%psi_R_lower, &
            certificate%psi_R_upper), &
            gc_outward_interval(certificate%psi_Z_lower, &
            certificate%psi_Z_upper), &
            gc_outward_interval(certificate%psi_RR_lower, &
            certificate%psi_RR_upper), &
            gc_outward_interval(certificate%psi_RZ_lower, &
            certificate%psi_RZ_upper), &
            gc_outward_interval(certificate%psi_ZZ_lower, &
            certificate%psi_ZZ_upper), interval_system(1), interval_system(2), &
            interval_system(3), interval_system(4), interval_system(5), &
            interval_system(6))
        if (.not. all_valid(interval_system)) then
            status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
            return
        end if

        point_R = r_lo+0.5_dp*(r_hi-r_lo)
        point_Z = z_lo+0.5_dp*(z_hi-z_lo)
        do iteration = 1, 8
            call evaluate_axis_point(point_R, point_Z, point_system, &
                point_status)
            if (point_status /= EQDSK_AXIS_CERT_SUCCESS) then
                status = point_status
                return
            end if
            call evaluate_neort_eqdsk_axis_stationarity_newton( &
                point_R, point_Z, point_system(1), point_system(2), &
                point_system(3), point_system(4), point_system(6), &
                determinant, next_R, next_Z, inverse(1), inverse(2), &
                inverse(3), inverse(4))
            if (.not. all(ieee_is_finite([determinant, next_R, next_Z, &
                    inverse])) .or. determinant == 0.0_dp) then
                status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
                return
            end if
            if (next_R <= r_lo .or. next_R >= r_hi .or. &
                    next_Z <= z_lo .or. next_Z >= z_hi) then
                status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
                return
            end if
            point_R = next_R
            point_Z = next_Z
        end do
        call evaluate_axis_point(point_R, point_Z, point_system, point_status)
        if (point_status /= EQDSK_AXIS_CERT_SUCCESS) then
            status = point_status
            return
        end if
        call evaluate_neort_eqdsk_axis_stationarity_newton(point_R, point_Z, &
            point_system(1), point_system(2), point_system(3), point_system(4), &
            point_system(6), determinant, next_R, next_Z, inverse(1), &
            inverse(2), inverse(3), inverse(4))
        if (.not. all(ieee_is_finite([determinant, next_R, next_Z, &
                inverse])) .or. determinant == 0.0_dp) then
            status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
            return
        end if
        delta_R = gc_outward_interval(r_lo, r_hi)-point(point_R)
        delta_Z = gc_outward_interval(z_lo, z_hi)-point(point_Z)
        call evaluate_neort_eqdsk_axis_stationarity_krawczyk_interval( &
            point(point_R), point(point_Z), point(point_system(1)), &
            point(point_system(2)), point(point_system(3)), &
            point(point_system(4)), point(point_system(6)), &
            interval_system(3), interval_system(4), interval_system(6), &
            delta_R, delta_Z, krawczyk_R, krawczyk_Z)
        if (.not. gc_outward_interval_is_valid(krawczyk_R) .or. &
                .not. gc_outward_interval_is_valid(krawczyk_Z)) then
            status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
            return
        end if
        if (krawczyk_R%lo <= r_lo .or. krawczyk_R%hi >= r_hi .or. &
                krawczyk_Z%lo <= z_lo .or. krawczyk_Z%hi >= z_hi) then
            status = EQDSK_AXIS_CERT_KRAWCZYK_FAILURE
            return
        end if
        certificate%axis_newton_R = point_R
        certificate%axis_newton_Z = point_Z
        certificate%axis_R_lo = krawczyk_R%lo
        certificate%axis_R_hi = krawczyk_R%hi
        certificate%axis_Z_lo = krawczyk_Z%lo
        certificate%axis_Z_hi = krawczyk_Z%hi
        certificate%axis_hessian_determinant = determinant
        certificate%stationarity_enclosure_certified = .true.
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
        if (.not. certificate%stationarity_enclosure_certified) return
        if (certificate%certificate_id /= EQDSK_AXIS_CERTIFICATE_ID) return
        if (certificate%cell_tiles_covered < 1) return
        if (.not. all(ieee_is_finite([certificate%r_lo, certificate%r_hi, &
                certificate%z_lo, certificate%z_hi, certificate%psi_hat_lo, &
                certificate%psi_hat_hi, certificate%psi_R_lower, &
                certificate%psi_R_upper, certificate%psi_Z_lower, &
                certificate%psi_Z_upper, certificate%psi_RR_lower, &
                certificate%psi_RR_upper, certificate%psi_RZ_lower, &
                certificate%psi_RZ_upper, certificate%psi_ZZ_lower, &
                certificate%psi_ZZ_upper, &
                certificate%hessian_determinant_lower, &
                certificate%cut_flux_curvature_lower, &
                certificate%left_psi_R_upper, &
                certificate%right_psi_R_lower, &
                certificate%bottom_psi_Z_upper, &
                certificate%top_psi_Z_lower, &
                certificate%left_dpsihat_dR_upper, &
                certificate%right_dpsihat_dR_lower, &
                certificate%axis_newton_R, certificate%axis_newton_Z, &
                certificate%axis_R_lo, certificate%axis_R_hi, &
                certificate%axis_Z_lo, certificate%axis_Z_hi, &
                certificate%axis_hessian_determinant]))) return
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
        if (certificate%psi_R_lower > certificate%psi_R_upper .or. &
                certificate%psi_Z_lower > certificate%psi_Z_upper .or. &
                certificate%psi_RR_lower > certificate%psi_RR_upper .or. &
                certificate%psi_RZ_lower > certificate%psi_RZ_upper .or. &
                certificate%psi_ZZ_lower > certificate%psi_ZZ_upper) return
        if (certificate%axis_hessian_determinant == 0.0_dp) return
        if (certificate%axis_newton_R <= certificate%r_lo .or. &
                certificate%axis_newton_R >= certificate%r_hi .or. &
                certificate%axis_newton_Z <= certificate%z_lo .or. &
                certificate%axis_newton_Z >= certificate%z_hi) return
        if (certificate%axis_R_lo <= certificate%r_lo .or. &
                certificate%axis_R_hi >= certificate%r_hi .or. &
                certificate%axis_Z_lo <= certificate%z_lo .or. &
                certificate%axis_Z_hi >= certificate%z_hi) return
        if (certificate%axis_R_lo > certificate%axis_R_hi .or. &
                certificate%axis_Z_lo > certificate%axis_Z_hi) return
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
        certificate%psi_R_lower = min(certificate%psi_R_lower, tile%psi_R%lo)
        certificate%psi_R_upper = max(certificate%psi_R_upper, tile%psi_R%hi)
        certificate%psi_Z_lower = min(certificate%psi_Z_lower, tile%psi_Z%lo)
        certificate%psi_Z_upper = max(certificate%psi_Z_upper, tile%psi_Z%hi)
        certificate%psi_RR_lower = min(certificate%psi_RR_lower, &
            tile%psi_RR%lo)
        certificate%psi_RR_upper = max(certificate%psi_RR_upper, &
            tile%psi_RR%hi)
        certificate%psi_RZ_lower = min(certificate%psi_RZ_lower, &
            tile%psi_RZ%lo)
        certificate%psi_RZ_upper = max(certificate%psi_RZ_upper, &
            tile%psi_RZ%hi)
        certificate%psi_ZZ_lower = min(certificate%psi_ZZ_lower, &
            tile%psi_ZZ%lo)
        certificate%psi_ZZ_upper = max(certificate%psi_ZZ_upper, &
            tile%psi_ZZ%hi)
        certificate%hessian_determinant_lower = min( &
            certificate%hessian_determinant_lower, &
            tile%axis_hessian_determinant%lo)
        certificate%cut_flux_curvature_lower = min( &
            certificate%cut_flux_curvature_lower, &
            tile%d2psihat_dR2%lo)
    end subroutine include_tile

    subroutine evaluate_axis_point(radius, height, system, status)
        real(dp), intent(in) :: radius, height
        real(dp), intent(out) :: system(6)
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: jet
        integer :: jet_status

        system = 0.0_dp
        call evaluate_eqdsk_cut_jet([radius, height, 0.0_dp], 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_AXIS_CERT_POINT_FAILURE
            return
        end if
        call evaluate_neort_eqdsk_axis_stationarity_system( &
            jet%psi_jet(2), jet%psi_jet(3), jet%psi_jet(4), &
            jet%psi_jet(5), jet%psi_jet(6), system(1), system(2), &
            system(3), system(4), system(5), system(6))
        if (.not. all(ieee_is_finite(system))) then
            status = EQDSK_AXIS_CERT_POINT_FAILURE
            return
        end if
        status = EQDSK_AXIS_CERT_SUCCESS
    end subroutine evaluate_axis_point

    logical function all_valid(values)
        type(gc_outward_interval_t), intent(in) :: values(:)
        integer :: i

        all_valid = .false.
        do i = 1, size(values)
            if (.not. gc_outward_interval_is_valid(values(i))) return
        end do
        all_valid = .true.
    end function all_valid

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

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
