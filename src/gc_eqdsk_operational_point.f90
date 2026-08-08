module neort_gc_eqdsk_operational_point
    !! Point-level orchestration for the direct EQDSK operational chart.
    !!
    !! This module owns only input validation, kernel sequencing, and status
    !! handling.  The cut geometry, allowed-energy jet, and canonical jet are
    !! evaluated exclusively by Fortsym-generated kernels.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_allowed_energy_symbolic, only: &
        evaluate_neort_eqdsk_allowed_energy
    use neort_eqdsk_canonical_cut_symbolic, only: &
        evaluate_neort_eqdsk_canonical_cut
    use neort_eqdsk_cut_r_flux_curvature_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_curvature
    use neort_gc_eqdsk_cut_jet, only: eqdsk_cut_jet_t
    implicit none
    private

    integer, parameter, public :: EQDSK_OPERATIONAL_POINT_SUCCESS = 0
    integer, parameter, public :: EQDSK_OPERATIONAL_POINT_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET = 2
    integer, parameter, public :: &
        EQDSK_OPERATIONAL_POINT_INVALID_CUT_DENOMINATOR = 3
    integer, parameter, public :: &
        EQDSK_OPERATIONAL_POINT_INVALID_CANONICAL_SQRT = 4
    integer, parameter, public :: EQDSK_OPERATIONAL_POINT_NONFINITE = 5

    type, public :: eqdsk_operational_point_t
        !! Values and physical-R jets returned by the three generated kernels.
        integer :: status = EQDSK_OPERATIONAL_POINT_INVALID_INPUT
        real(dp) :: dZ_dR = 0.0_dp
        real(dp) :: d2Z_dR2 = 0.0_dp
        real(dp) :: dpsihat_dR = 0.0_dp
        real(dp) :: d2psihat_dR2 = 0.0_dp
        real(dp) :: field_norm_squared = 0.0_dp
        real(dp) :: psi_physical = 0.0_dp
        real(dp) :: dpsi_physical_dR = 0.0_dp
        real(dp) :: d2psi_physical_dR2 = 0.0_dp
        real(dp) :: bmod = 0.0_dp
        real(dp) :: dbmod_dR = 0.0_dp
        real(dp) :: d2bmod_dR2 = 0.0_dp
        real(dp) :: omega_c = 0.0_dp
        real(dp) :: domega_c_dR = 0.0_dp
        real(dp) :: d2omega_c_dR2 = 0.0_dp
        real(dp) :: energy_margin = 0.0_dp
        real(dp) :: denergy_margin_dR = 0.0_dp
        real(dp) :: d2energy_margin_dR2 = 0.0_dp
        real(dp) :: v_parallel_squared = 0.0_dp
        real(dp) :: dv_parallel_squared_dR = 0.0_dp
        real(dp) :: d2v_parallel_squared_dR2 = 0.0_dp
        real(dp) :: bphi_covariant = 0.0_dp
        real(dp) :: dbphi_covariant_dR = 0.0_dp
        real(dp) :: d2bphi_covariant_dR2 = 0.0_dp
        real(dp) :: v_parallel = 0.0_dp
        real(dp) :: dv_parallel_dR = 0.0_dp
        real(dp) :: d2v_parallel_dR2 = 0.0_dp
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: dpsi_star_dR = 0.0_dp
        real(dp) :: d2psi_star_dR2 = 0.0_dp
    end type eqdsk_operational_point_t

    public :: evaluate_eqdsk_operational_point

contains

    subroutine evaluate_eqdsk_operational_point(cut_jet, radius, field_scale, &
            raw_psi_sep, h, j_k, mass, charge, c_light, electrostatic_potential, &
            dphi_dpsi, d2phi_dpsi2, sigma, result, status)
        type(eqdsk_cut_jet_t), intent(in) :: cut_jet
        real(dp), intent(in) :: radius, field_scale, raw_psi_sep, h, j_k
        real(dp), intent(in) :: mass, charge, c_light, electrostatic_potential
        real(dp), intent(in) :: dphi_dpsi, d2phi_dpsi2
        integer, intent(in) :: sigma
        type(eqdsk_operational_point_t), intent(out) :: result
        integer, intent(out) :: status

        real(dp) :: dZ_dR, d2Z_dR2, dpsihat_dR, d2psihat_dR2
        real(dp) :: field_norm_squared, psi_physical, dpsi_physical_dR
        real(dp) :: d2psi_physical_dR2, bmod, dbmod_dR, d2bmod_dR2
        real(dp) :: omega_c, domega_c_dR, d2omega_c_dR2
        real(dp) :: energy_margin, denergy_margin_dR, d2energy_margin_dR2
        real(dp) :: v_parallel_squared, dv_parallel_squared_dR
        real(dp) :: d2v_parallel_squared_dR2, bphi_covariant
        real(dp) :: dbphi_covariant_dR, d2bphi_covariant_dR2
        real(dp) :: v_parallel, dv_parallel_dR, d2v_parallel_dR2
        real(dp) :: psi_star, dpsi_star_dR, d2psi_star_dR2
        real(dp) :: scalar_inputs(12), cut_inputs(3)

        result = eqdsk_operational_point_t()
        status = EQDSK_OPERATIONAL_POINT_INVALID_INPUT
        scalar_inputs = [radius, field_scale, raw_psi_sep, h, j_k, mass, &
            charge, c_light, electrostatic_potential, dphi_dpsi, d2phi_dpsi2, &
            real(sigma, dp)]
        if (.not. all(ieee_is_finite(scalar_inputs))) return
        if (radius <= 0.0_dp .or. field_scale <= 0.0_dp .or. &
                raw_psi_sep <= 0.0_dp .or. mass <= 0.0_dp .or. &
                c_light <= 0.0_dp .or. j_k < 0.0_dp .or. &
                abs(charge) <= tiny(1.0_dp)) return
        if (abs(sigma) /= 1) return

        cut_inputs = [cut_jet%d_cut_numerator_d_R, &
            cut_jet%d_cut_numerator_d_Z, cut_jet%cut_numerator]
        if (.not. all(ieee_is_finite(cut_jet%psi_jet))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET)
            return
        end if
        if (.not. all(ieee_is_finite(cut_jet%f_jet))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET)
            return
        end if
        if (.not. all(ieee_is_finite([cut_inputs, &
                cut_jet%d2_cut_numerator_d_R2, &
                cut_jet%d2_cut_numerator_d_RdZ, &
                cut_jet%d2_cut_numerator_d_Z2, cut_jet%cut_value, &
                cut_jet%d_cut_d_R, cut_jet%d_cut_d_arc_phi, &
                cut_jet%d_cut_d_Z, cut_jet%cut_rate, &
                cut_jet%absolute_cut_rate, cut_jet%orientation_scalar]))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_INVALID_CUT_JET)
            return
        end if
        if (abs(cut_jet%d_cut_numerator_d_Z) <= &
                64.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
                abs(cut_jet%d_cut_numerator_d_R), &
                abs(cut_jet%d_cut_numerator_d_Z))) then
            call fail(result, status, &
                EQDSK_OPERATIONAL_POINT_INVALID_CUT_DENOMINATOR)
            return
        end if

        call evaluate_neort_eqdsk_cut_r_flux_curvature( &
            cut_jet%d_cut_numerator_d_R, cut_jet%d_cut_numerator_d_Z, &
            cut_jet%d2_cut_numerator_d_R2, cut_jet%d2_cut_numerator_d_RdZ, &
            cut_jet%d2_cut_numerator_d_Z2, cut_jet%psi_jet(2), &
            cut_jet%psi_jet(3), cut_jet%psi_jet(4), cut_jet%psi_jet(5), &
            cut_jet%psi_jet(6), raw_psi_sep, dZ_dR, d2Z_dR2, dpsihat_dR, &
            d2psihat_dR2)
        if (.not. all(ieee_is_finite([dZ_dR, d2Z_dR2, dpsihat_dR, &
                d2psihat_dR2]))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_NONFINITE)
            return
        end if

        call evaluate_neort_eqdsk_allowed_energy(radius, field_scale, &
            cut_jet%psi_jet(1), cut_jet%psi_jet(2), cut_jet%psi_jet(3), &
            cut_jet%psi_jet(4), cut_jet%psi_jet(5), cut_jet%psi_jet(6), &
            cut_jet%psi_jet(7), cut_jet%psi_jet(8), cut_jet%psi_jet(9), &
            cut_jet%psi_jet(10), cut_jet%f_jet(1), cut_jet%f_jet(2), &
            cut_jet%f_jet(3), raw_psi_sep, dZ_dR, d2Z_dR2, h, j_k, mass, &
            charge, c_light, electrostatic_potential, dphi_dpsi, &
            d2phi_dpsi2, field_norm_squared, psi_physical, dpsi_physical_dR, &
            bmod, dbmod_dR, d2bmod_dR2, omega_c, domega_c_dR, &
            d2omega_c_dR2, energy_margin, denergy_margin_dR, &
            d2energy_margin_dR2, v_parallel_squared, &
            dv_parallel_squared_dR, d2v_parallel_squared_dR2, &
            bphi_covariant, dbphi_covariant_dR, d2psi_physical_dR2, &
            d2bphi_covariant_dR2)
        if (.not. all(ieee_is_finite([field_norm_squared, psi_physical, &
                dpsi_physical_dR, d2psi_physical_dR2, bmod, dbmod_dR, &
                d2bmod_dR2, omega_c, domega_c_dR, d2omega_c_dR2, &
                energy_margin, denergy_margin_dR, d2energy_margin_dR2, &
                v_parallel_squared, dv_parallel_squared_dR, &
                d2v_parallel_squared_dR2, bphi_covariant, &
                dbphi_covariant_dR, d2bphi_covariant_dR2]))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_NONFINITE)
            return
        end if
        if (v_parallel_squared <= 0.0_dp) then
            call fail(result, status, &
                EQDSK_OPERATIONAL_POINT_INVALID_CANONICAL_SQRT)
            return
        end if

        call evaluate_neort_eqdsk_canonical_cut(v_parallel_squared, &
            dv_parallel_squared_dR, mass, charge, c_light, real(sigma, dp), &
            psi_physical, dpsi_physical_dR, bphi_covariant, &
            dbphi_covariant_dR, d2v_parallel_squared_dR2, &
            d2psi_physical_dR2, d2bphi_covariant_dR2, v_parallel, &
            dv_parallel_dR, psi_star, dpsi_star_dR, d2v_parallel_dR2, &
            d2psi_star_dR2)
        if (.not. all(ieee_is_finite([v_parallel, dv_parallel_dR, &
                d2v_parallel_dR2, psi_star, dpsi_star_dR, &
                d2psi_star_dR2]))) then
            call fail(result, status, EQDSK_OPERATIONAL_POINT_NONFINITE)
            return
        end if

        result%dZ_dR = dZ_dR
        result%d2Z_dR2 = d2Z_dR2
        result%dpsihat_dR = dpsihat_dR
        result%d2psihat_dR2 = d2psihat_dR2
        result%field_norm_squared = field_norm_squared
        result%psi_physical = psi_physical
        result%dpsi_physical_dR = dpsi_physical_dR
        result%d2psi_physical_dR2 = d2psi_physical_dR2
        result%bmod = bmod
        result%dbmod_dR = dbmod_dR
        result%d2bmod_dR2 = d2bmod_dR2
        result%omega_c = omega_c
        result%domega_c_dR = domega_c_dR
        result%d2omega_c_dR2 = d2omega_c_dR2
        result%energy_margin = energy_margin
        result%denergy_margin_dR = denergy_margin_dR
        result%d2energy_margin_dR2 = d2energy_margin_dR2
        result%v_parallel_squared = v_parallel_squared
        result%dv_parallel_squared_dR = dv_parallel_squared_dR
        result%d2v_parallel_squared_dR2 = d2v_parallel_squared_dR2
        result%bphi_covariant = bphi_covariant
        result%dbphi_covariant_dR = dbphi_covariant_dR
        result%d2bphi_covariant_dR2 = d2bphi_covariant_dR2
        result%v_parallel = v_parallel
        result%dv_parallel_dR = dv_parallel_dR
        result%d2v_parallel_dR2 = d2v_parallel_dR2
        result%psi_star = psi_star
        result%dpsi_star_dR = dpsi_star_dR
        result%d2psi_star_dR2 = d2psi_star_dR2
        status = EQDSK_OPERATIONAL_POINT_SUCCESS
        result%status = status
    end subroutine evaluate_eqdsk_operational_point

    subroutine fail(result, status, code)
        type(eqdsk_operational_point_t), intent(inout) :: result
        integer, intent(out) :: status
        integer, intent(in) :: code

        result = eqdsk_operational_point_t()
        result%status = code
        status = code
    end subroutine fail

end module neort_gc_eqdsk_operational_point
