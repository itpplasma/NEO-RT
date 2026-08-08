module neort_gc_eqdsk_cut_endpoint_certificate
    !! Existence-and-uniqueness certificate for a regular Eq.13/flux endpoint.
    !!
    !! The physical residual, Jacobian, Newton inverse, and Krawczyk map are
    !! Fortsym-generated.  This module owns only grid tiling, interval hulling,
    !! generated-kernel orchestration, and theorem acceptance.  A regular
    !! certificate requires K(X) strictly inside X.  The magnetic axis is
    !! excluded because its flux coordinate has a separate square-root limit.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, psi_sep, rad, zet
    use neort_eqdsk_cut_endpoint_krawczyk_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval
    use neort_eqdsk_cut_endpoint_newton_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_newton
    use neort_eqdsk_cut_endpoint_system_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_system_interval
    use neort_eqdsk_cut_endpoint_system_symbolic, only: &
        evaluate_neort_eqdsk_cut_endpoint_system
    use neort_gc_eqdsk_cut_interval, only: &
        EQDSK_CUT_INTERVAL_SUCCESS, eqdsk_cut_interval_result_t, &
        evaluate_eqdsk_cut_interval_box
    use neort_gc_eqdsk_cut_jet, only: &
        EQDSK_CUT_JET_SUCCESS, eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_contains_zero, gc_outward_interval_is_valid, &
        gc_outward_interval_t, operator(-)
    implicit none
    private

    integer, parameter, public :: EQDSK_ENDPOINT_CERT_SUCCESS = 0
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_GRID_INVALID = 2
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_INTERVAL_FAILURE = 3
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_DOMAIN_DEGENERATE = 4
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_POINT_FAILURE = 5
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_SINGULAR = 6
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_NEWTON_ESCAPE = 7
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_NO_ROOT_CANDIDATE = 8
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_KRAWCZYK_FAILURE = 9
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_NOT_INCLUDED = 10
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_INVALID_CERTIFICATE = 11
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_AXIS_REQUIRES_LIMIT = 12
    integer, parameter, public :: EQDSK_ENDPOINT_CERT_NONFINITE = 13
    integer, parameter, public :: EQDSK_ENDPOINT_CERTIFICATE_ID = 130015

    type, public :: eqdsk_cut_endpoint_options_t
        integer :: max_newton_iterations = 8
    end type eqdsk_cut_endpoint_options_t

    type, public :: eqdsk_cut_endpoint_certificate_t
        real(dp) :: r_lo = 0.0_dp
        real(dp) :: r_hi = 0.0_dp
        real(dp) :: z_lo = 0.0_dp
        real(dp) :: z_hi = 0.0_dp
        real(dp) :: target_psihat = 0.0_dp
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: separatrix_flux = 0.0_dp
        real(dp) :: positive_denominator_lower = 0.0_dp
        real(dp) :: newton_point_R = 0.0_dp
        real(dp) :: newton_point_Z = 0.0_dp
        real(dp) :: jacobian_determinant = 0.0_dp
        real(dp) :: cut_residual_lo = 0.0_dp
        real(dp) :: cut_residual_hi = 0.0_dp
        real(dp) :: flux_residual_lo = 0.0_dp
        real(dp) :: flux_residual_hi = 0.0_dp
        real(dp) :: krawczyk_R_lo = 0.0_dp
        real(dp) :: krawczyk_R_hi = 0.0_dp
        real(dp) :: krawczyk_Z_lo = 0.0_dp
        real(dp) :: krawczyk_Z_hi = 0.0_dp
        integer :: cell_tiles_covered = 0
        integer :: newton_iterations = 0
        integer :: certificate_id = 0
        logical :: regular_domain_certified = .false.
        logical :: strict_inclusion_certified = .false.
        logical :: unique_endpoint_certified = .false.
    end type eqdsk_cut_endpoint_certificate_t

    public :: build_eqdsk_cut_endpoint_certificate
    public :: validate_eqdsk_cut_endpoint_certificate

contains

    subroutine build_eqdsk_cut_endpoint_certificate(r_lo, r_hi, z_lo, z_hi, &
            target_psihat, field_scale, seed_R, seed_Z, options, certificate, &
            status)
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi
        real(dp), intent(in) :: target_psihat, field_scale, seed_R, seed_Z
        type(eqdsk_cut_endpoint_options_t), intent(in) :: options
        type(eqdsk_cut_endpoint_certificate_t), intent(out) :: certificate
        integer, intent(out) :: status

        type(gc_outward_interval_t) :: box_system(6)
        type(gc_outward_interval_t) :: delta_R, delta_Z
        type(gc_outward_interval_t) :: krawczyk_R, krawczyk_Z
        real(dp) :: point_system(6), point_R, point_Z
        real(dp) :: determinant, next_R, next_Z, inverse(4)
        real(dp) :: denominator_lower
        integer :: iteration, point_status, tiles

        certificate = eqdsk_cut_endpoint_certificate_t()
        status = EQDSK_ENDPOINT_CERT_INVALID_INPUT
        if (.not. all(ieee_is_finite([r_lo, r_hi, z_lo, z_hi, &
                target_psihat, field_scale, seed_R, seed_Z]))) then
            status = EQDSK_ENDPOINT_CERT_NONFINITE
            return
        end if
        if (target_psihat == 0.0_dp) then
            status = EQDSK_ENDPOINT_CERT_AXIS_REQUIRES_LIMIT
            return
        end if
        if (r_lo <= 0.0_dp .or. r_lo >= r_hi .or. z_lo >= z_hi) return
        if (target_psihat < 0.0_dp .or. target_psihat > 1.0_dp) return
        if (field_scale <= 0.0_dp) return
        if (seed_R <= r_lo .or. seed_R >= r_hi .or. &
                seed_Z <= z_lo .or. seed_Z >= z_hi) return
        if (options%max_newton_iterations < 0) return
        status = validate_grid()
        if (status /= EQDSK_ENDPOINT_CERT_SUCCESS) return
        if (r_lo < rad(1) .or. r_hi > rad(nrad) .or. &
                z_lo < zet(1) .or. z_hi > zet(nzet)) then
            status = EQDSK_ENDPOINT_CERT_INVALID_INPUT
            return
        end if

        call evaluate_endpoint_box_hull(r_lo, r_hi, z_lo, z_hi, &
            target_psihat, box_system, denominator_lower, tiles, status)
        if (status /= EQDSK_ENDPOINT_CERT_SUCCESS) return
        if (.not. gc_outward_interval_contains_zero(box_system(1)) .or. &
                .not. gc_outward_interval_contains_zero(box_system(2))) then
            status = EQDSK_ENDPOINT_CERT_NO_ROOT_CANDIDATE
            return
        end if

        point_R = seed_R
        point_Z = seed_Z
        do iteration = 1, options%max_newton_iterations
            call evaluate_point_system(point_R, point_Z, target_psihat, &
                field_scale, point_system, point_status)
            if (point_status /= EQDSK_ENDPOINT_CERT_SUCCESS) then
                status = point_status
                return
            end if
            call evaluate_neort_eqdsk_cut_endpoint_newton(point_R, point_Z, &
                point_system(1), point_system(2), point_system(3), &
                point_system(4), point_system(5), point_system(6), &
                determinant, next_R, next_Z, inverse(1), inverse(2), &
                inverse(3), inverse(4))
            if (.not. all(ieee_is_finite([determinant, next_R, next_Z, &
                    inverse]))) then
                status = EQDSK_ENDPOINT_CERT_SINGULAR
                return
            end if
            if (determinant == 0.0_dp) then
                status = EQDSK_ENDPOINT_CERT_SINGULAR
                return
            end if
            if (next_R <= r_lo .or. next_R >= r_hi .or. &
                    next_Z <= z_lo .or. next_Z >= z_hi) then
                status = EQDSK_ENDPOINT_CERT_NEWTON_ESCAPE
                return
            end if
            point_R = next_R
            point_Z = next_Z
        end do

        call evaluate_point_system(point_R, point_Z, target_psihat, &
            field_scale, point_system, point_status)
        if (point_status /= EQDSK_ENDPOINT_CERT_SUCCESS) then
            status = point_status
            return
        end if
        call evaluate_neort_eqdsk_cut_endpoint_newton(point_R, point_Z, &
            point_system(1), point_system(2), point_system(3), &
            point_system(4), point_system(5), point_system(6), determinant, &
            next_R, next_Z, inverse(1), inverse(2), inverse(3), inverse(4))
        if (.not. all(ieee_is_finite([determinant, next_R, next_Z, &
                inverse])) .or. determinant == 0.0_dp) then
            status = EQDSK_ENDPOINT_CERT_SINGULAR
            return
        end if

        delta_R = gc_outward_interval(r_lo, r_hi)-point(point_R)
        delta_Z = gc_outward_interval(z_lo, z_hi)-point(point_Z)
        call evaluate_neort_eqdsk_cut_endpoint_krawczyk_interval( &
            point(point_R), point(point_Z), point(point_system(1)), &
            point(point_system(2)), point(point_system(3)), &
            point(point_system(4)), point(point_system(5)), &
            point(point_system(6)), box_system(3), box_system(4), &
            box_system(5), box_system(6), delta_R, delta_Z, krawczyk_R, &
            krawczyk_Z)
        if (.not. gc_outward_interval_is_valid(krawczyk_R) .or. &
                .not. gc_outward_interval_is_valid(krawczyk_Z)) then
            status = EQDSK_ENDPOINT_CERT_KRAWCZYK_FAILURE
            return
        end if
        if (krawczyk_R%lo <= r_lo .or. krawczyk_R%hi >= r_hi .or. &
                krawczyk_Z%lo <= z_lo .or. krawczyk_Z%hi >= z_hi) then
            status = EQDSK_ENDPOINT_CERT_NOT_INCLUDED
            return
        end if

        certificate%r_lo = r_lo
        certificate%r_hi = r_hi
        certificate%z_lo = z_lo
        certificate%z_hi = z_hi
        certificate%target_psihat = target_psihat
        certificate%field_scale = field_scale
        certificate%separatrix_flux = psi_sep
        certificate%positive_denominator_lower = denominator_lower
        certificate%newton_point_R = point_R
        certificate%newton_point_Z = point_Z
        certificate%jacobian_determinant = determinant
        certificate%cut_residual_lo = box_system(1)%lo
        certificate%cut_residual_hi = box_system(1)%hi
        certificate%flux_residual_lo = box_system(2)%lo
        certificate%flux_residual_hi = box_system(2)%hi
        certificate%krawczyk_R_lo = krawczyk_R%lo
        certificate%krawczyk_R_hi = krawczyk_R%hi
        certificate%krawczyk_Z_lo = krawczyk_Z%lo
        certificate%krawczyk_Z_hi = krawczyk_Z%hi
        certificate%cell_tiles_covered = tiles
        certificate%newton_iterations = options%max_newton_iterations
        certificate%regular_domain_certified = .true.
        certificate%strict_inclusion_certified = .true.
        certificate%unique_endpoint_certified = .true.
        certificate%certificate_id = EQDSK_ENDPOINT_CERTIFICATE_ID
        call validate_eqdsk_cut_endpoint_certificate(certificate, status)
    end subroutine build_eqdsk_cut_endpoint_certificate

    subroutine validate_eqdsk_cut_endpoint_certificate(certificate, status)
        type(eqdsk_cut_endpoint_certificate_t), intent(in) :: certificate
        integer, intent(out) :: status

        status = EQDSK_ENDPOINT_CERT_INVALID_CERTIFICATE
        if (.not. certificate%regular_domain_certified) return
        if (.not. certificate%strict_inclusion_certified) return
        if (.not. certificate%unique_endpoint_certified) return
        if (certificate%certificate_id /= EQDSK_ENDPOINT_CERTIFICATE_ID) return
        if (certificate%cell_tiles_covered < 1) return
        if (certificate%newton_iterations < 0) return
        if (.not. all(ieee_is_finite([certificate%r_lo, certificate%r_hi, &
                certificate%z_lo, certificate%z_hi, &
                certificate%target_psihat, certificate%field_scale, &
                certificate%separatrix_flux, &
                certificate%positive_denominator_lower, &
                certificate%newton_point_R, certificate%newton_point_Z, &
                certificate%jacobian_determinant, &
                certificate%cut_residual_lo, certificate%cut_residual_hi, &
                certificate%flux_residual_lo, &
                certificate%flux_residual_hi, certificate%krawczyk_R_lo, &
                certificate%krawczyk_R_hi, certificate%krawczyk_Z_lo, &
                certificate%krawczyk_Z_hi]))) then
            status = EQDSK_ENDPOINT_CERT_NONFINITE
            return
        end if
        if (certificate%r_lo <= 0.0_dp .or. &
                certificate%r_lo >= certificate%r_hi .or. &
                certificate%z_lo >= certificate%z_hi) return
        if (certificate%target_psihat <= 0.0_dp .or. &
                certificate%target_psihat > 1.0_dp) return
        if (certificate%field_scale <= 0.0_dp .or. &
                certificate%separatrix_flux <= tiny(1.0_dp)) return
        if (certificate%positive_denominator_lower <= 0.0_dp) return
        if (certificate%jacobian_determinant == 0.0_dp) return
        if (certificate%cut_residual_lo > 0.0_dp .or. &
                certificate%cut_residual_hi < 0.0_dp) return
        if (certificate%flux_residual_lo > 0.0_dp .or. &
                certificate%flux_residual_hi < 0.0_dp) return
        if (certificate%newton_point_R <= certificate%r_lo .or. &
                certificate%newton_point_R >= certificate%r_hi .or. &
                certificate%newton_point_Z <= certificate%z_lo .or. &
                certificate%newton_point_Z >= certificate%z_hi) return
        if (certificate%krawczyk_R_lo <= certificate%r_lo .or. &
                certificate%krawczyk_R_hi >= certificate%r_hi .or. &
                certificate%krawczyk_Z_lo <= certificate%z_lo .or. &
                certificate%krawczyk_Z_hi >= certificate%z_hi) return
        if (certificate%krawczyk_R_lo > certificate%krawczyk_R_hi .or. &
                certificate%krawczyk_Z_lo > certificate%krawczyk_Z_hi) return
        status = EQDSK_ENDPOINT_CERT_SUCCESS
    end subroutine validate_eqdsk_cut_endpoint_certificate

    subroutine evaluate_endpoint_box_hull(r_lo, r_hi, z_lo, z_hi, &
            target_psihat, system, denominator_lower, tiles, status)
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi, target_psihat
        type(gc_outward_interval_t), intent(out) :: system(6)
        real(dp), intent(out) :: denominator_lower
        integer, intent(out) :: tiles, status

        type(eqdsk_cut_interval_result_t) :: tile
        type(gc_outward_interval_t) :: tile_system(6)
        real(dp) :: tile_r_lo, tile_r_hi, tile_z_lo, tile_z_hi
        integer :: cell_R, cell_Z, interval_status
        logical :: initialized

        system = point(0.0_dp)
        denominator_lower = 0.0_dp
        tiles = 0
        initialized = .false.
        status = EQDSK_ENDPOINT_CERT_INTERVAL_FAILURE
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
                if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) return
                if (.not. tile%denominator_positive_certified) then
                    status = EQDSK_ENDPOINT_CERT_DOMAIN_DEGENERATE
                    return
                end if
                call evaluate_neort_eqdsk_cut_endpoint_system_interval( &
                    tile%numerator, tile%numerator_R, tile%numerator_Z, &
                    tile%psi, tile%psi_R, tile%psi_Z, point(psi_sep), &
                    point(target_psihat), tile_system(1), tile_system(2), &
                    tile_system(3), tile_system(4), tile_system(5), &
                    tile_system(6))
                if (.not. all_valid(tile_system)) then
                    status = EQDSK_ENDPOINT_CERT_NONFINITE
                    return
                end if
                if (.not. initialized) then
                    system = tile_system
                    denominator_lower = tile%positive_denominator%lo
                    initialized = .true.
                else
                    call include_hull(system, tile_system)
                    denominator_lower = min(denominator_lower, &
                        tile%positive_denominator%lo)
                end if
                tiles = tiles+1
            end do
        end do
        if (.not. initialized .or. tiles < 1) return
        if (denominator_lower <= 0.0_dp) then
            status = EQDSK_ENDPOINT_CERT_DOMAIN_DEGENERATE
            return
        end if
        status = EQDSK_ENDPOINT_CERT_SUCCESS
    end subroutine evaluate_endpoint_box_hull

    subroutine evaluate_point_system(radius, height, target_psihat, &
            field_scale, system, status)
        real(dp), intent(in) :: radius, height, target_psihat, field_scale
        real(dp), intent(out) :: system(6)
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: jet
        integer :: jet_status

        system = 0.0_dp
        call evaluate_eqdsk_cut_jet([radius, height, 0.0_dp], field_scale, &
            1, [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_ENDPOINT_CERT_POINT_FAILURE
            return
        end if
        call evaluate_neort_eqdsk_cut_endpoint_system(jet%cut_numerator, &
            jet%d_cut_numerator_d_R, jet%d_cut_numerator_d_Z, &
            jet%psi_jet(1), jet%psi_jet(2), jet%psi_jet(3), psi_sep, &
            target_psihat, system(1), system(2), system(3), system(4), &
            system(5), system(6))
        if (.not. all(ieee_is_finite(system))) then
            status = EQDSK_ENDPOINT_CERT_NONFINITE
            return
        end if
        status = EQDSK_ENDPOINT_CERT_SUCCESS
    end subroutine evaluate_point_system

    subroutine include_hull(hull, values)
        type(gc_outward_interval_t), intent(inout) :: hull(:)
        type(gc_outward_interval_t), intent(in) :: values(:)
        integer :: i

        do i = 1, size(hull)
            hull(i)%lo = min(hull(i)%lo, values(i)%lo)
            hull(i)%hi = max(hull(i)%hi, values(i)%hi)
        end do
    end subroutine include_hull

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

    integer function validate_grid()
        integer :: i

        validate_grid = EQDSK_ENDPOINT_CERT_GRID_INVALID
        if (nrad < 2 .or. nzet < 2) return
        if (.not. allocated(rad) .or. .not. allocated(zet)) return
        if (size(rad) /= nrad .or. size(zet) /= nzet) return
        if (.not. all(ieee_is_finite(rad)) .or. &
                .not. all(ieee_is_finite(zet))) return
        if (.not. ieee_is_finite(psi_sep) .or. psi_sep <= tiny(1.0_dp)) return
        do i = 1, nrad-1
            if (rad(i+1) <= rad(i)) return
        end do
        do i = 1, nzet-1
            if (zet(i+1) <= zet(i)) return
        end do
        validate_grid = EQDSK_ENDPOINT_CERT_SUCCESS
    end function validate_grid

end module neort_gc_eqdsk_cut_endpoint_certificate
