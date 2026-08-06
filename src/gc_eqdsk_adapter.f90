module neort_gc_eqdsk_adapter
    !! Adapter from the initialized direct-GEQDSK backend to the coordinate-
    !! general guiding-center field interface.  Sampling never changes the
    !! legacy do_magfie current-surface cache.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: sample_eqdsk_field
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, GC_MODEL_INVALID_STATE, gc_field_t

    implicit none
    private

    type, extends(gc_field_t), public :: eqdsk_gc_field_t
        real(dp) :: field_scale = 1.0_dp
    contains
        procedure :: evaluate => evaluate_eqdsk_gc_field
        procedure :: radial_profiles => evaluate_eqdsk_profiles
    end type eqdsk_gc_field_t

contains

    subroutine evaluate_eqdsk_gc_field(self, position, sample, status)
        class(eqdsk_gc_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp) :: q_value, dqds_value, psi_pol, grad_psi_pol(3), psi_tor_edge

        sample = gc_field_sample_t()
        status = GC_MODEL_INVALID_STATE
        if (position(1) < 0.0_dp .or. position(1) > 1.0_dp) return
        if (self%field_scale <= 0.0_dp) return

        call sample_eqdsk_field(position, self%field_scale, sample%bmod, &
            sample%sqrtg, sample%grad_log_b, sample%hcov, sample%hcon, &
            sample%curl_h, q_value, dqds_value, psi_pol, grad_psi_pol, &
            psi_tor_edge)
        associate (unused_q => q_value, unused_dqds => dqds_value, &
                unused_tor_flux => psi_tor_edge)
        end associate
        sample%psi = psi_pol
        sample%grad_psi = grad_psi_pol
        if (.not. all(ieee_is_finite([sample%bmod, sample%sqrtg, &
            sample%grad_log_b, sample%hcov, sample%hcon, &
            sample%curl_h, sample%psi, sample%grad_psi]))) return
        if (sample%bmod <= 0.0_dp .or. &
            abs(sample%sqrtg) <= tiny(sample%sqrtg)) return
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_eqdsk_gc_field

    subroutine evaluate_eqdsk_profiles(self, radial, q_value, dqds_value, &
            psi_pol, dpsi_pol_ds, psi_tor_edge, status)
        class(eqdsk_gc_field_t), intent(in) :: self
        real(dp), intent(in) :: radial
        real(dp), intent(out) :: q_value, dqds_value, psi_pol
        real(dp), intent(out) :: dpsi_pol_ds, psi_tor_edge
        integer, intent(out) :: status

        real(dp) :: bmod, sqrtg, grad_log_b(3), hcov(3), hcon(3), curl_h(3)
        real(dp) :: grad_psi_pol(3)

        q_value = 0.0_dp
        dqds_value = 0.0_dp
        psi_pol = 0.0_dp
        dpsi_pol_ds = 0.0_dp
        psi_tor_edge = 0.0_dp
        status = GC_MODEL_INVALID_STATE
        if (radial < 0.0_dp .or. radial > 1.0_dp) return
        if (self%field_scale <= 0.0_dp) return
        call sample_eqdsk_field([radial, 0.0_dp, 0.0_dp], self%field_scale, &
            bmod, sqrtg, grad_log_b, hcov, hcon, curl_h, q_value, &
            dqds_value, psi_pol, grad_psi_pol, psi_tor_edge)
        dpsi_pol_ds = grad_psi_pol(1)
        if (.not. all(ieee_is_finite([bmod, sqrtg, q_value, dqds_value, &
            psi_pol, dpsi_pol_ds, psi_tor_edge]))) return
        if (bmod <= 0.0_dp .or. abs(sqrtg) <= tiny(sqrtg) &
            .or. abs(q_value) <= tiny(q_value)) return
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_eqdsk_profiles

end module neort_gc_eqdsk_adapter
