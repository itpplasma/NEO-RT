module neort_gc_eqdsk_allowed_region_interval
    !! Certified allowed-energy and canonical-coordinate jets on an Eq. 13
    !! cut graph.  Runtime code only validates, intersects profile cells, and
    !! hulls outward intervals.  Every physical expression is emitted by
    !! Fortsym.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_allowed_energy_interval_symbolic, only: &
        evaluate_neort_eqdsk_allowed_energy_interval
    use neort_eqdsk_canonical_cut_interval_symbolic, only: &
        evaluate_neort_eqdsk_canonical_cut_interval
    use neort_eqdsk_physical_flux_map_interval_symbolic, only: &
        evaluate_neort_eqdsk_physical_flux_map_interval
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t
    use neort_profile_potential_map_interval_symbolic, only: &
        evaluate_neort_profile_potential_map_interval
    implicit none
    private

    integer, parameter, public :: EQDSK_ALLOWED_INTERVAL_SUCCESS = 0
    integer, parameter, public :: EQDSK_ALLOWED_INTERVAL_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_ALLOWED_INTERVAL_PROFILE_GAP = 2
    integer, parameter, public :: EQDSK_ALLOWED_INTERVAL_UNCERTIFIED_CUT = 3
    integer, parameter, public :: EQDSK_ALLOWED_INTERVAL_NONFINITE = 4

    type, public :: eqdsk_allowed_interval_result_t
        type(gc_outward_interval_t) :: field_norm_squared
        type(gc_outward_interval_t) :: psi_physical
        type(gc_outward_interval_t) :: dpsi_physical_dR
        type(gc_outward_interval_t) :: bmod
        type(gc_outward_interval_t) :: dbmod_dR
        type(gc_outward_interval_t) :: d2bmod_dR2
        type(gc_outward_interval_t) :: omega_c
        type(gc_outward_interval_t) :: domega_c_dR
        type(gc_outward_interval_t) :: d2omega_c_dR2
        type(gc_outward_interval_t) :: energy_margin
        type(gc_outward_interval_t) :: denergy_margin_dR
        type(gc_outward_interval_t) :: d2energy_margin_dR2
        type(gc_outward_interval_t) :: v_parallel_squared
        type(gc_outward_interval_t) :: dv_parallel_squared_dR
        type(gc_outward_interval_t) :: d2v_parallel_squared_dR2
        type(gc_outward_interval_t) :: bphi_covariant
        type(gc_outward_interval_t) :: dbphi_covariant_dR
        type(gc_outward_interval_t) :: v_parallel
        type(gc_outward_interval_t) :: dv_parallel_dR
        type(gc_outward_interval_t) :: psi_star
        type(gc_outward_interval_t) :: dpsi_star_dR
        logical :: canonical_chart_certified = .false.
        integer :: cut_enclosures_covered = 0
        integer :: profile_segments_covered = 0
    end type eqdsk_allowed_interval_result_t

    public :: evaluate_eqdsk_allowed_region_interval
    public :: evaluate_eqdsk_profile_potential_interval

contains

    subroutine evaluate_eqdsk_profile_potential_interval(psi, psi_nodes, &
            phi_nodes, omega_nodes, c_light, potential, dpotential_dpsi, &
            d2potential_dpsi2, segments_covered, status)
        type(gc_outward_interval_t), intent(in) :: psi
        real(dp), intent(in) :: psi_nodes(:), phi_nodes(:), omega_nodes(:)
        real(dp), intent(in) :: c_light
        type(gc_outward_interval_t), intent(out) :: potential
        type(gc_outward_interval_t), intent(out) :: dpotential_dpsi
        type(gc_outward_interval_t), intent(out) :: d2potential_dpsi2
        integer, intent(out) :: segments_covered, status

        type(gc_outward_interval_t) :: query, local_potential
        type(gc_outward_interval_t) :: local_first, local_second
        integer :: i, node_count
        logical :: have_segment

        potential = gc_outward_interval(0.0_dp, 0.0_dp)
        dpotential_dpsi = potential
        d2potential_dpsi2 = potential
        segments_covered = 0
        status = EQDSK_ALLOWED_INTERVAL_INVALID_INPUT
        if (.not. gc_outward_interval_is_valid(psi)) return
        if (.not. ieee_is_finite(c_light) .or. c_light <= 0.0_dp) return
        node_count = size(psi_nodes)
        if (node_count < 2) return
        if (size(phi_nodes) /= node_count) return
        if (size(omega_nodes) /= node_count) return
        if (.not. all(ieee_is_finite(psi_nodes))) return
        if (.not. all(ieee_is_finite(phi_nodes))) return
        if (.not. all(ieee_is_finite(omega_nodes))) return
        do i = 1, node_count-1
            if (psi_nodes(i+1) <= psi_nodes(i)) return
        end do
        status = EQDSK_ALLOWED_INTERVAL_PROFILE_GAP
        if (psi%lo < psi_nodes(1)) return
        if (psi%hi > psi_nodes(node_count)) return

        have_segment = .false.
        do i = 1, node_count-1
            query = gc_outward_interval(max(psi%lo, psi_nodes(i)), &
                min(psi%hi, psi_nodes(i+1)))
            if (query%hi < query%lo) cycle
            call evaluate_neort_profile_potential_map_interval(query, &
                point(psi_nodes(i)), point(psi_nodes(i+1)), &
                point(phi_nodes(i)), point(omega_nodes(i)), &
                point(omega_nodes(i+1)), point(c_light), local_potential, &
                local_first, local_second)
            if (.not. valid_intervals([local_potential, local_first, &
                    local_second])) then
                status = EQDSK_ALLOWED_INTERVAL_NONFINITE
                return
            end if
            if (.not. have_segment) then
                potential = local_potential
                dpotential_dpsi = local_first
                d2potential_dpsi2 = local_second
                have_segment = .true.
            else
                potential = interval_hull(potential, local_potential)
                dpotential_dpsi = interval_hull(dpotential_dpsi, local_first)
                d2potential_dpsi2 = interval_hull(d2potential_dpsi2, &
                    local_second)
            end if
            segments_covered = segments_covered+1
        end do
        if (.not. have_segment) return
        status = EQDSK_ALLOWED_INTERVAL_SUCCESS
    end subroutine evaluate_eqdsk_profile_potential_interval

    subroutine evaluate_eqdsk_allowed_region_interval(radius, field_scale, &
            separatrix_psi, cut_enclosures, psi_nodes, phi_nodes, omega_nodes, &
            h0, j_k, mass, charge, c_light, parallel_sign, result, status)
        type(gc_outward_interval_t), intent(in) :: radius
        real(dp), intent(in) :: field_scale, separatrix_psi
        type(eqdsk_cut_interval_result_t), intent(in) :: cut_enclosures(:)
        real(dp), intent(in) :: psi_nodes(:), phi_nodes(:), omega_nodes(:)
        real(dp), intent(in) :: h0, j_k, mass, charge, c_light
        integer, intent(in) :: parallel_sign
        type(eqdsk_allowed_interval_result_t), intent(out) :: result
        integer, intent(out) :: status

        type(eqdsk_allowed_interval_result_t) :: candidate
        type(gc_outward_interval_t) :: potential, dpotential_dpsi
        type(gc_outward_interval_t) :: d2potential_dpsi2
        type(gc_outward_interval_t) :: mapped_psi, mapped_dpsi_dR
        integer :: i, local_status, segments
        logical :: have_candidate

        result = eqdsk_allowed_interval_result_t()
        status = EQDSK_ALLOWED_INTERVAL_INVALID_INPUT
        if (.not. gc_outward_interval_is_valid(radius)) return
        if (radius%lo <= 0.0_dp) return
        if (.not. all(ieee_is_finite([field_scale, separatrix_psi, h0, &
                j_k, mass, charge, c_light]))) return
        if (field_scale <= 0.0_dp .or. separatrix_psi <= 0.0_dp) return
        if (mass <= 0.0_dp .or. c_light <= 0.0_dp) return
        if (abs(charge) <= tiny(charge) .or. j_k < 0.0_dp) return
        if (abs(parallel_sign) /= 1) return
        if (size(cut_enclosures) < 1) return

        have_candidate = .false.
        result%canonical_chart_certified = .true.
        do i = 1, size(cut_enclosures)
            if (.not. cut_enclosures(i)%denominator_positive_certified .or. &
                    .not. cut_enclosures(i)%r_chart_certified) then
                status = EQDSK_ALLOWED_INTERVAL_UNCERTIFIED_CUT
                return
            end if
            call evaluate_neort_eqdsk_physical_flux_map_interval( &
                cut_enclosures(i)%psi, cut_enclosures(i)%psi_R, &
                cut_enclosures(i)%psi_Z, cut_enclosures(i)%dZ_dR, &
                point(field_scale), mapped_psi, mapped_dpsi_dR)
            if (.not. valid_intervals([mapped_psi, mapped_dpsi_dR])) then
                status = EQDSK_ALLOWED_INTERVAL_NONFINITE
                return
            end if
            call evaluate_eqdsk_profile_potential_interval(mapped_psi, &
                psi_nodes, phi_nodes, omega_nodes, c_light, potential, &
                dpotential_dpsi, d2potential_dpsi2, segments, local_status)
            if (local_status /= EQDSK_ALLOWED_INTERVAL_SUCCESS) then
                status = local_status
                return
            end if
            candidate = eqdsk_allowed_interval_result_t()
            call evaluate_energy_candidate(radius, field_scale, &
                separatrix_psi, cut_enclosures(i), potential, &
                dpotential_dpsi, d2potential_dpsi2, h0, j_k, mass, charge, &
                c_light, parallel_sign, candidate, local_status)
            if (local_status /= EQDSK_ALLOWED_INTERVAL_SUCCESS) then
                status = local_status
                return
            end if
            candidate%profile_segments_covered = segments
            if (.not. have_candidate) then
                result = candidate
                have_candidate = .true.
            else
                call merge_allowed_results(result, candidate)
            end if
            result%cut_enclosures_covered = i
        end do
        if (.not. have_candidate) return
        if (.not. result%canonical_chart_certified) then
            result%v_parallel = point(0.0_dp)
            result%dv_parallel_dR = point(0.0_dp)
            result%psi_star = point(0.0_dp)
            result%dpsi_star_dR = point(0.0_dp)
        end if
        status = EQDSK_ALLOWED_INTERVAL_SUCCESS
    end subroutine evaluate_eqdsk_allowed_region_interval

    subroutine evaluate_energy_candidate(radius, field_scale, separatrix_psi, &
            cut, potential, dpotential_dpsi, d2potential_dpsi2, h0, j_k, &
            mass, charge, c_light, parallel_sign, result, status)
        type(gc_outward_interval_t), intent(in) :: radius
        real(dp), intent(in) :: field_scale, separatrix_psi
        type(eqdsk_cut_interval_result_t), intent(in) :: cut
        type(gc_outward_interval_t), intent(in) :: potential
        type(gc_outward_interval_t), intent(in) :: dpotential_dpsi
        type(gc_outward_interval_t), intent(in) :: d2potential_dpsi2
        real(dp), intent(in) :: h0, j_k, mass, charge, c_light
        integer, intent(in) :: parallel_sign
        type(eqdsk_allowed_interval_result_t), intent(out) :: result
        integer, intent(out) :: status

        result = eqdsk_allowed_interval_result_t()
        call evaluate_neort_eqdsk_allowed_energy_interval(radius, &
            point(field_scale), cut%psi, cut%psi_R, cut%psi_Z, cut%psi_RR, &
            cut%psi_RZ, cut%psi_ZZ, cut%psi_RRR, cut%psi_RRZ, &
            cut%psi_RZZ, cut%psi_ZZZ, cut%F, cut%dF_dpsi_hat, &
            cut%d2F_dpsi_hat2, point(separatrix_psi), cut%dZ_dR, &
            cut%d2Z_dR2, point(h0), point(j_k), point(mass), point(charge), &
            point(c_light), potential, dpotential_dpsi, d2potential_dpsi2, &
            result%field_norm_squared, result%psi_physical, &
            result%dpsi_physical_dR, result%bmod, result%dbmod_dR, &
            result%d2bmod_dR2, result%omega_c, result%domega_c_dR, &
            result%d2omega_c_dR2, result%energy_margin, &
            result%denergy_margin_dR, result%d2energy_margin_dR2, &
            result%v_parallel_squared, result%dv_parallel_squared_dR, &
            result%d2v_parallel_squared_dR2, result%bphi_covariant, &
            result%dbphi_covariant_dR)
        status = EQDSK_ALLOWED_INTERVAL_NONFINITE
        if (.not. valid_energy_result(result)) return
        result%canonical_chart_certified = &
            result%v_parallel_squared%lo > 0.0_dp
        if (result%canonical_chart_certified) then
            call evaluate_neort_eqdsk_canonical_cut_interval( &
                result%v_parallel_squared, result%dv_parallel_squared_dR, &
                point(mass), point(charge), point(c_light), &
                point(real(parallel_sign, dp)), result%psi_physical, &
                result%dpsi_physical_dR, result%bphi_covariant, &
                result%dbphi_covariant_dR, result%v_parallel, &
                result%dv_parallel_dR, result%psi_star, result%dpsi_star_dR)
            if (.not. valid_intervals([result%v_parallel, &
                    result%dv_parallel_dR, result%psi_star, &
                    result%dpsi_star_dR])) return
        end if
        status = EQDSK_ALLOWED_INTERVAL_SUCCESS
    end subroutine evaluate_energy_candidate

    subroutine merge_allowed_results(result, candidate)
        type(eqdsk_allowed_interval_result_t), intent(inout) :: result
        type(eqdsk_allowed_interval_result_t), intent(in) :: candidate

        result%field_norm_squared = interval_hull( &
            result%field_norm_squared, candidate%field_norm_squared)
        result%psi_physical = interval_hull(result%psi_physical, &
            candidate%psi_physical)
        result%dpsi_physical_dR = interval_hull(result%dpsi_physical_dR, &
            candidate%dpsi_physical_dR)
        result%bmod = interval_hull(result%bmod, candidate%bmod)
        result%dbmod_dR = interval_hull(result%dbmod_dR, candidate%dbmod_dR)
        result%d2bmod_dR2 = interval_hull(result%d2bmod_dR2, &
            candidate%d2bmod_dR2)
        result%omega_c = interval_hull(result%omega_c, candidate%omega_c)
        result%domega_c_dR = interval_hull(result%domega_c_dR, &
            candidate%domega_c_dR)
        result%d2omega_c_dR2 = interval_hull(result%d2omega_c_dR2, &
            candidate%d2omega_c_dR2)
        result%energy_margin = interval_hull(result%energy_margin, &
            candidate%energy_margin)
        result%denergy_margin_dR = interval_hull( &
            result%denergy_margin_dR, candidate%denergy_margin_dR)
        result%d2energy_margin_dR2 = interval_hull( &
            result%d2energy_margin_dR2, candidate%d2energy_margin_dR2)
        result%v_parallel_squared = interval_hull( &
            result%v_parallel_squared, candidate%v_parallel_squared)
        result%dv_parallel_squared_dR = interval_hull( &
            result%dv_parallel_squared_dR, &
            candidate%dv_parallel_squared_dR)
        result%d2v_parallel_squared_dR2 = interval_hull( &
            result%d2v_parallel_squared_dR2, &
            candidate%d2v_parallel_squared_dR2)
        result%bphi_covariant = interval_hull(result%bphi_covariant, &
            candidate%bphi_covariant)
        result%dbphi_covariant_dR = interval_hull( &
            result%dbphi_covariant_dR, candidate%dbphi_covariant_dR)
        result%profile_segments_covered = &
            result%profile_segments_covered+candidate%profile_segments_covered
        if (result%canonical_chart_certified .and. &
                candidate%canonical_chart_certified) then
            result%v_parallel = interval_hull(result%v_parallel, &
                candidate%v_parallel)
            result%dv_parallel_dR = interval_hull(result%dv_parallel_dR, &
                candidate%dv_parallel_dR)
            result%psi_star = interval_hull(result%psi_star, &
                candidate%psi_star)
            result%dpsi_star_dR = interval_hull(result%dpsi_star_dR, &
                candidate%dpsi_star_dR)
        else
            result%canonical_chart_certified = .false.
        end if
    end subroutine merge_allowed_results

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    pure function interval_hull(left, right) result(hull)
        type(gc_outward_interval_t), intent(in) :: left, right
        type(gc_outward_interval_t) :: hull

        hull = gc_outward_interval(min(left%lo, right%lo), &
            max(left%hi, right%hi))
    end function interval_hull

    pure logical function valid_intervals(values)
        type(gc_outward_interval_t), intent(in) :: values(:)

        valid_intervals = all(gc_outward_interval_is_valid(values))
    end function valid_intervals

    pure logical function valid_energy_result(result)
        type(eqdsk_allowed_interval_result_t), intent(in) :: result

        valid_energy_result = valid_intervals([result%field_norm_squared, &
            result%psi_physical, result%dpsi_physical_dR, result%bmod, &
            result%dbmod_dR, result%d2bmod_dR2, result%omega_c, &
            result%domega_c_dR, result%d2omega_c_dR2, &
            result%energy_margin, result%denergy_margin_dR, &
            result%d2energy_margin_dR2, result%v_parallel_squared, &
            result%dv_parallel_squared_dR, &
            result%d2v_parallel_squared_dR2, result%bphi_covariant, &
            result%dbphi_covariant_dR])
    end function valid_energy_result

end module neort_gc_eqdsk_allowed_region_interval
