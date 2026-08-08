module neort_gc_eqdsk_flux_profile_map
    !! Certified coordinate map between native toroidal flux and the
    !! normalized poloidal flux used by the direct-EQDSK cut atlas.
    !!
    !! This module owns only profile-segment selection and certificate state.
    !! Every forward, inverse, and rho_tor derivative formula is emitted by
    !! Fortsym.  In particular, s_tor and psi_hat are never identified.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_flux_profile_segment_symbolic, only: &
        evaluate_neort_eqdsk_flux_profile_segment
    use neort_eqdsk_flux_profile_segment_interval_symbolic, only: &
        evaluate_neort_eqdsk_flux_profile_segment_interval
    use neort_eqdsk_flux_profile_rho_chain_symbolic, only: &
        evaluate_neort_eqdsk_flux_profile_rho_chain
    use neort_eqdsk_rho_tor_map_symbolic, only: &
        evaluate_neort_eqdsk_rho_tor_map
    use neort_eqdsk_scaled_flux_normalization_symbolic, only: &
        evaluate_neort_eqdsk_scaled_flux_normalization
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t
    implicit none
    private

    integer, parameter, public :: EQDSK_FLUX_MAP_SUCCESS = 0
    integer, parameter, public :: EQDSK_FLUX_MAP_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_FLUX_MAP_NONFINITE = 2
    integer, parameter, public :: EQDSK_FLUX_MAP_NOT_MONOTONE = 3
    integer, parameter, public :: EQDSK_FLUX_MAP_ENDPOINT_MISMATCH = 4
    integer, parameter, public :: EQDSK_FLUX_MAP_UNINITIALIZED = 5
    integer, parameter, public :: EQDSK_FLUX_MAP_OUT_OF_RANGE = 6
    integer, parameter, public :: EQDSK_FLUX_MAP_INVALID_CERTIFICATE = 7
    integer, parameter, public :: EQDSK_FLUX_MAP_CERTIFICATE_ID = 130015

    type, public :: eqdsk_flux_profile_map_t
        real(dp), allocatable :: s_tor(:)
        real(dp), allocatable :: scaled_psi(:)
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: psi_sep = 0.0_dp
        real(dp) :: psihat_axis = 0.0_dp
        real(dp) :: psihat_edge = 0.0_dp
        real(dp) :: minimum_dpsihat_dstor = 0.0_dp
        integer :: certificate_id = 0
        logical :: endpoints_certified = .false.
        logical :: strict_monotonicity_certified = .false.
        logical :: initialized = .false.
        character(len=64) :: method = ''
    end type eqdsk_flux_profile_map_t

    public :: clear_eqdsk_flux_profile_map
    public :: initialize_eqdsk_flux_profile_map
    public :: validate_eqdsk_flux_profile_map
    public :: map_eqdsk_rho_tor_to_psihat
    public :: map_eqdsk_s_tor_to_psihat
    public :: map_eqdsk_psihat_to_s_tor
    public :: map_eqdsk_psihat_interval_to_s_tor
    public :: map_eqdsk_scaled_psi_to_s_tor

contains

    subroutine initialize_eqdsk_flux_profile_map(s_tor, scaled_psi, &
            field_scale, separatrix_psi, map, status)
        real(dp), intent(in) :: s_tor(:), scaled_psi(:)
        real(dp), intent(in) :: field_scale, separatrix_psi
        type(eqdsk_flux_profile_map_t), intent(inout) :: map
        integer, intent(out) :: status

        real(dp) :: psihat, derivative, inverse_s, tolerance
        integer :: i, count

        call clear_eqdsk_flux_profile_map(map)
        status = EQDSK_FLUX_MAP_INVALID_INPUT
        count = size(s_tor)
        if (count < 2 .or. size(scaled_psi) /= count) return
        if (.not. all(ieee_is_finite([s_tor, scaled_psi, field_scale, &
                separatrix_psi]))) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        if (field_scale <= 0.0_dp .or. separatrix_psi <= 0.0_dp) return
        if (s_tor(1) /= 0.0_dp .or. s_tor(count) /= 1.0_dp) then
            status = EQDSK_FLUX_MAP_ENDPOINT_MISMATCH
            return
        end if
        do i = 1, count-1
            if (s_tor(i+1) <= s_tor(i)) then
                status = EQDSK_FLUX_MAP_NOT_MONOTONE
                return
            end if
        end do

        allocate(map%s_tor(count), map%scaled_psi(count))
        map%s_tor = s_tor
        map%scaled_psi = scaled_psi
        map%field_scale = field_scale
        map%psi_sep = separatrix_psi
        map%minimum_dpsihat_dstor = huge(1.0_dp)
        do i = 1, count-1
            call evaluate_neort_eqdsk_flux_profile_segment(s_tor(i), 0.0_dp, &
                s_tor(i), s_tor(i+1), scaled_psi(i), scaled_psi(i+1), &
                field_scale, separatrix_psi, psihat, derivative, inverse_s)
            if (.not. all(ieee_is_finite([psihat, derivative, inverse_s]))) then
                status = EQDSK_FLUX_MAP_NONFINITE
                call clear_eqdsk_flux_profile_map(map)
                return
            end if
            if (derivative <= 0.0_dp) then
                status = EQDSK_FLUX_MAP_NOT_MONOTONE
                call clear_eqdsk_flux_profile_map(map)
                return
            end if
            map%minimum_dpsihat_dstor = min( &
                map%minimum_dpsihat_dstor, derivative)
        end do

        call evaluate_neort_eqdsk_flux_profile_segment(s_tor(1), 0.0_dp, &
            s_tor(1), s_tor(2), scaled_psi(1), scaled_psi(2), field_scale, &
            separatrix_psi, map%psihat_axis, derivative, inverse_s)
        call evaluate_neort_eqdsk_flux_profile_segment(s_tor(count), 1.0_dp, &
            s_tor(count-1), s_tor(count), scaled_psi(count-1), &
            scaled_psi(count), field_scale, separatrix_psi, map%psihat_edge, &
            derivative, inverse_s)
        if (.not. all(ieee_is_finite([map%psihat_axis, map%psihat_edge]))) then
            status = EQDSK_FLUX_MAP_NONFINITE
            call clear_eqdsk_flux_profile_map(map)
            return
        end if
        tolerance = 256.0_dp*epsilon(1.0_dp)
        if (abs(map%psihat_axis) > tolerance .or. &
                abs(map%psihat_edge-1.0_dp) > tolerance) then
            status = EQDSK_FLUX_MAP_ENDPOINT_MISMATCH
            call clear_eqdsk_flux_profile_map(map)
            return
        end if

        map%endpoints_certified = .true.
        map%strict_monotonicity_certified = .true.
        map%certificate_id = EQDSK_FLUX_MAP_CERTIFICATE_ID
        map%method = 'fortsym-piecewise-affine-s_tor-to-psi_hat'
        map%initialized = .true.
        call validate_eqdsk_flux_profile_map(map, status)
    end subroutine initialize_eqdsk_flux_profile_map

    subroutine clear_eqdsk_flux_profile_map(map)
        type(eqdsk_flux_profile_map_t), intent(inout) :: map

        if (allocated(map%s_tor)) deallocate(map%s_tor)
        if (allocated(map%scaled_psi)) deallocate(map%scaled_psi)
        map%field_scale = 0.0_dp
        map%psi_sep = 0.0_dp
        map%psihat_axis = 0.0_dp
        map%psihat_edge = 0.0_dp
        map%minimum_dpsihat_dstor = 0.0_dp
        map%certificate_id = 0
        map%endpoints_certified = .false.
        map%strict_monotonicity_certified = .false.
        map%initialized = .false.
        map%method = ''
    end subroutine clear_eqdsk_flux_profile_map

    subroutine validate_eqdsk_flux_profile_map(map, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        integer, intent(out) :: status

        real(dp) :: derivative, inverse_s, minimum_derivative
        real(dp) :: psihat_axis, psihat_edge, segment_psihat
        integer :: i, count

        status = EQDSK_FLUX_MAP_INVALID_CERTIFICATE
        if (.not. map%initialized) then
            status = EQDSK_FLUX_MAP_UNINITIALIZED
            return
        end if
        if (.not. map%endpoints_certified .or. &
                .not. map%strict_monotonicity_certified) return
        if (map%certificate_id /= EQDSK_FLUX_MAP_CERTIFICATE_ID) return
        if (.not. allocated(map%s_tor) .or. &
                .not. allocated(map%scaled_psi)) return
        count = size(map%s_tor)
        if (count < 2 .or. size(map%scaled_psi) /= count) return
        if (.not. all(ieee_is_finite([map%s_tor, map%scaled_psi, &
                map%field_scale, map%psi_sep, map%psihat_axis, &
                map%psihat_edge, map%minimum_dpsihat_dstor]))) return
        if (map%field_scale <= 0.0_dp .or. map%psi_sep <= 0.0_dp) return
        if (map%s_tor(1) /= 0.0_dp .or. map%s_tor(count) /= 1.0_dp) return
        do i = 1, count-1
            if (map%s_tor(i+1) <= map%s_tor(i)) return
            call evaluate_neort_eqdsk_flux_profile_segment(map%s_tor(i), &
                0.0_dp, map%s_tor(i), map%s_tor(i+1), &
                map%scaled_psi(i), map%scaled_psi(i+1), map%field_scale, &
                map%psi_sep, segment_psihat, derivative, inverse_s)
            if (.not. all(ieee_is_finite([segment_psihat, derivative, &
                    inverse_s]))) return
            if (derivative <= 0.0_dp) return
            if (i == 1) then
                minimum_derivative = derivative
                psihat_axis = segment_psihat
            else
                minimum_derivative = min(minimum_derivative, derivative)
            end if
        end do
        if (map%minimum_dpsihat_dstor <= 0.0_dp) return
        if (map%minimum_dpsihat_dstor /= minimum_derivative) return
        call evaluate_neort_eqdsk_flux_profile_segment(map%s_tor(count), &
            1.0_dp, map%s_tor(count-1), map%s_tor(count), &
            map%scaled_psi(count-1), map%scaled_psi(count), &
            map%field_scale, map%psi_sep, psihat_edge, derivative, inverse_s)
        if (.not. all(ieee_is_finite([psihat_axis, psihat_edge, &
                derivative, inverse_s]))) return
        if (map%psihat_axis /= psihat_axis) return
        if (map%psihat_edge /= psihat_edge) return
        if (map%psihat_axis /= 0.0_dp .or. map%psihat_edge /= 1.0_dp) then
            if (abs(map%psihat_axis) > 256.0_dp*epsilon(1.0_dp) .or. &
                    abs(map%psihat_edge-1.0_dp) > &
                    256.0_dp*epsilon(1.0_dp)) return
        end if
        status = EQDSK_FLUX_MAP_SUCCESS
    end subroutine validate_eqdsk_flux_profile_map

    subroutine map_eqdsk_s_tor_to_psihat(map, s_tor, psihat, &
            dpsihat_dstor, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: s_tor
        real(dp), intent(out) :: psihat, dpsihat_dstor
        integer, intent(out) :: status

        real(dp) :: unused_inverse
        integer :: segment

        psihat = 0.0_dp
        dpsihat_dstor = 0.0_dp
        call validate_eqdsk_flux_profile_map(map, status)
        if (status /= EQDSK_FLUX_MAP_SUCCESS) return
        if (.not. ieee_is_finite(s_tor)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        segment = locate_s_segment(map, s_tor)
        if (segment < 1) then
            status = EQDSK_FLUX_MAP_OUT_OF_RANGE
            return
        end if
        call evaluate_neort_eqdsk_flux_profile_segment(s_tor, 0.0_dp, &
            map%s_tor(segment), map%s_tor(segment+1), &
            map%scaled_psi(segment), map%scaled_psi(segment+1), &
            map%field_scale, map%psi_sep, psihat, dpsihat_dstor, &
            unused_inverse)
        if (.not. all(ieee_is_finite([psihat, dpsihat_dstor]))) then
            psihat = 0.0_dp
            dpsihat_dstor = 0.0_dp
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        if (dpsihat_dstor <= 0.0_dp) then
            status = EQDSK_FLUX_MAP_INVALID_CERTIFICATE
            return
        end if
        status = EQDSK_FLUX_MAP_SUCCESS
    end subroutine map_eqdsk_s_tor_to_psihat

    subroutine map_eqdsk_psihat_to_s_tor(map, psihat, s_tor, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: psihat
        real(dp), intent(out) :: s_tor
        integer, intent(out) :: status

        real(dp) :: endpoint_psihat, unused_derivative
        integer :: segment, count

        s_tor = 0.0_dp
        call validate_eqdsk_flux_profile_map(map, status)
        if (status /= EQDSK_FLUX_MAP_SUCCESS) return
        if (.not. ieee_is_finite(psihat)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        if (psihat < 0.0_dp .or. psihat > 1.0_dp) then
            status = EQDSK_FLUX_MAP_OUT_OF_RANGE
            return
        end if
        count = size(map%s_tor)
        segment = 0
        do segment = 1, count-1
            call evaluate_neort_eqdsk_flux_profile_segment( &
                map%s_tor(segment+1), psihat, map%s_tor(segment), &
                map%s_tor(segment+1), map%scaled_psi(segment), &
                map%scaled_psi(segment+1), map%field_scale, map%psi_sep, &
                endpoint_psihat, unused_derivative, s_tor)
            if (.not. all(ieee_is_finite([endpoint_psihat, s_tor]))) then
                s_tor = 0.0_dp
                status = EQDSK_FLUX_MAP_NONFINITE
                return
            end if
            if (psihat <= endpoint_psihat .or. segment == count-1) exit
        end do
        if (.not. ieee_is_finite(s_tor)) then
            s_tor = 0.0_dp
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        status = EQDSK_FLUX_MAP_SUCCESS
    end subroutine map_eqdsk_psihat_to_s_tor

    subroutine map_eqdsk_psihat_interval_to_s_tor(map, psihat, s_tor, &
            segments_covered, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        type(gc_outward_interval_t), intent(in) :: psihat
        type(gc_outward_interval_t), intent(out) :: s_tor
        integer, intent(out) :: segments_covered
        integer, intent(out) :: status

        type(gc_outward_interval_t) :: segment_psihat, derivative
        type(gc_outward_interval_t) :: inverse_s, overlap, unused_inverse
        integer :: i, count
        logical :: first_piece

        s_tor = gc_outward_interval(0.0_dp, 0.0_dp)
        segments_covered = 0
        call validate_eqdsk_flux_profile_map(map, status)
        if (status /= EQDSK_FLUX_MAP_SUCCESS) return
        if (.not. gc_outward_interval_is_valid(psihat)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        if (psihat%lo < 0.0_dp .or. psihat%hi > 1.0_dp) then
            status = EQDSK_FLUX_MAP_OUT_OF_RANGE
            return
        end if

        count = size(map%s_tor)
        first_piece = .true.
        do i = 1, count-1
            call evaluate_neort_eqdsk_flux_profile_segment_interval( &
                gc_outward_interval(map%s_tor(i), map%s_tor(i+1)), &
                gc_outward_interval(0.0_dp, 0.0_dp), &
                gc_outward_interval(map%s_tor(i), map%s_tor(i)), &
                gc_outward_interval(map%s_tor(i+1), map%s_tor(i+1)), &
                gc_outward_interval(map%scaled_psi(i), &
                    map%scaled_psi(i)), &
                gc_outward_interval(map%scaled_psi(i+1), &
                    map%scaled_psi(i+1)), &
                gc_outward_interval(map%field_scale, map%field_scale), &
                gc_outward_interval(map%psi_sep, map%psi_sep), &
                segment_psihat, derivative, unused_inverse)
            if (.not. gc_outward_interval_is_valid(segment_psihat)) then
                status = EQDSK_FLUX_MAP_NONFINITE
                return
            end if
            if (.not. gc_outward_interval_is_valid(derivative)) then
                status = EQDSK_FLUX_MAP_NONFINITE
                return
            end if
            if (derivative%lo <= 0.0_dp) then
                status = EQDSK_FLUX_MAP_INVALID_CERTIFICATE
                return
            end if
            overlap%lo = max(psihat%lo, segment_psihat%lo)
            overlap%hi = min(psihat%hi, segment_psihat%hi)
            if (overlap%hi < overlap%lo) cycle
            call evaluate_neort_eqdsk_flux_profile_segment_interval( &
                gc_outward_interval(map%s_tor(i), map%s_tor(i)), overlap, &
                gc_outward_interval(map%s_tor(i), map%s_tor(i)), &
                gc_outward_interval(map%s_tor(i+1), map%s_tor(i+1)), &
                gc_outward_interval(map%scaled_psi(i), &
                    map%scaled_psi(i)), &
                gc_outward_interval(map%scaled_psi(i+1), &
                    map%scaled_psi(i+1)), &
                gc_outward_interval(map%field_scale, map%field_scale), &
                gc_outward_interval(map%psi_sep, map%psi_sep), &
                segment_psihat, derivative, inverse_s)
            if (.not. gc_outward_interval_is_valid(inverse_s)) then
                status = EQDSK_FLUX_MAP_NONFINITE
                return
            end if
            inverse_s%lo = max(inverse_s%lo, map%s_tor(i))
            inverse_s%hi = min(inverse_s%hi, map%s_tor(i+1))
            if (inverse_s%hi < inverse_s%lo) cycle
            if (first_piece) then
                s_tor = inverse_s
                first_piece = .false.
            else
                s_tor%lo = min(s_tor%lo, inverse_s%lo)
                s_tor%hi = max(s_tor%hi, inverse_s%hi)
            end if
            segments_covered = segments_covered+1
        end do
        if (first_piece) then
            status = EQDSK_FLUX_MAP_OUT_OF_RANGE
            return
        end if
        if (s_tor%lo < 0.0_dp .or. s_tor%hi > 1.0_dp) then
            status = EQDSK_FLUX_MAP_INVALID_CERTIFICATE
            return
        end if
        status = EQDSK_FLUX_MAP_SUCCESS
    end subroutine map_eqdsk_psihat_interval_to_s_tor

    subroutine map_eqdsk_scaled_psi_to_s_tor(map, scaled_psi, s_tor, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: scaled_psi
        real(dp), intent(out) :: s_tor
        integer, intent(out) :: status

        real(dp) :: psihat

        s_tor = 0.0_dp
        call validate_eqdsk_flux_profile_map(map, status)
        if (status /= EQDSK_FLUX_MAP_SUCCESS) return
        if (.not. ieee_is_finite(scaled_psi)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        call evaluate_neort_eqdsk_scaled_flux_normalization(scaled_psi, &
            map%field_scale, map%psi_sep, psihat)
        if (.not. ieee_is_finite(psihat)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        call map_eqdsk_psihat_to_s_tor(map, psihat, s_tor, status)
    end subroutine map_eqdsk_scaled_psi_to_s_tor

    subroutine map_eqdsk_rho_tor_to_psihat(map, rho_tor, s_tor, psihat, &
            dpsihat_drho_tor, status)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: rho_tor
        real(dp), intent(out) :: s_tor, psihat, dpsihat_drho_tor
        integer, intent(out) :: status

        real(dp) :: dstor_drho_tor, dpsihat_dstor

        s_tor = 0.0_dp
        psihat = 0.0_dp
        dpsihat_drho_tor = 0.0_dp
        if (.not. ieee_is_finite(rho_tor)) then
            status = EQDSK_FLUX_MAP_NONFINITE
            return
        end if
        if (rho_tor < 0.0_dp .or. rho_tor > 1.0_dp) then
            status = EQDSK_FLUX_MAP_OUT_OF_RANGE
            return
        end if
        call evaluate_neort_eqdsk_rho_tor_map(rho_tor, s_tor, &
            dstor_drho_tor)
        call map_eqdsk_s_tor_to_psihat(map, s_tor, psihat, &
            dpsihat_dstor, status)
        if (status /= EQDSK_FLUX_MAP_SUCCESS) return
        call evaluate_neort_eqdsk_flux_profile_rho_chain(dstor_drho_tor, &
            dpsihat_dstor, dpsihat_drho_tor)
        if (.not. ieee_is_finite(dpsihat_drho_tor)) then
            dpsihat_drho_tor = 0.0_dp
            status = EQDSK_FLUX_MAP_NONFINITE
        end if
    end subroutine map_eqdsk_rho_tor_to_psihat

    integer function locate_s_segment(map, value) result(segment)
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: value

        integer :: low, high, middle, count

        segment = 0
        count = size(map%s_tor)
        if (value < map%s_tor(1) .or. value > map%s_tor(count)) return
        if (value == map%s_tor(count)) then
            segment = count-1
            return
        end if
        low = 1
        high = count
        do while (high-low > 1)
            middle = (low+high)/2
            if (map%s_tor(middle) <= value) then
                low = middle
            else
                high = middle
            end if
        end do
        segment = low
    end function locate_s_segment

end module neort_gc_eqdsk_flux_profile_map
