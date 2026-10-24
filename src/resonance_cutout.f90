! Optional regularization of the collisionless resonance weight near folds.
!
! The collisionless quadrature in compute_transport_integral replaces
! delta(R(eta)) by 1/|R_eta| at each root of R = mth*Omega_theta +
! mph*Omega_phi. Near a fold of R (two roots merging where R_eta = 0) and next
! to the trapped-passing boundary this weight diverges. The delta-function
! limit is only valid when R stays linear across the collisional width of the
! resonance. With pitch-angle scattering D_eta, the characteristic widths are
!   delta_lin = (D_eta / |R_eta|)**(1/3),
! and a quadratic fold, R = R_eta_eta*(eta - eta_f)**2/2 + ..., has the width
!   delta_fold = (D_eta / |R_eta_eta|)**(1/4).
! Their order-unity normalization is fixed by the linear-layer-fit criterion,
! rather than a solved collisional response. A root is inside when its distance,
! d = |R_eta / R_eta_eta|, is smaller than delta_fold. Both widths agree on
! the layer boundary d = delta_fold, so the criterion is the point where the
! linear (delta-function) layer stops fitting between the two merging roots.
!
! resonance_cutout selects what replaces 1/|R_eta| for roots inside the layer
! (cutout_weight returns the replacement and whether it applies):
!   0  off (default, bit-identical to the plain collisionless quadrature),
!   1  zero contribution,
!   2  boundary-value fill, 1/|R_eta| frozen at its value on the layer boundary,
!      1/(|R_eta_eta| * w * delta_fold).
! resonance_cutout_width = w >= 1 widens the layer, e.g. as a bound on the
! additional finite-orbit-width broadening that this local model omits.
! See doc/resonance_cutout.md.
module neort_resonance_cutout
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_collisional_layer, only: comp_ellint_k, comp_ellint_e

    implicit none
    private

    integer, parameter, public :: CUTOUT_OFF = 0, CUTOUT_ZERO = 1, CUTOUT_FILL = 2

    integer, public :: resonance_cutout = CUTOUT_OFF
    real(dp), public :: resonance_cutout_width = 1.0_dp

    public :: validate_resonance_cutout, reset_resonance_cutout
    public :: bounce_pitch_diffusion, fold_layer_width, cutout_weight, cutout_curvature

    ! Keeps K(m) finite next to the trapped-passing boundary (m -> 1).
    real(dp), parameter :: M_MAX = 1.0_dp - 1.0e-12_dp

contains

    subroutine reset_resonance_cutout()
        resonance_cutout = CUTOUT_OFF
        resonance_cutout_width = 1.0_dp
    end subroutine reset_resonance_cutout

    subroutine validate_resonance_cutout()
        if (resonance_cutout < CUTOUT_OFF .or. resonance_cutout > CUTOUT_FILL) &
            error stop 'resonance_cutout must be 0 (off), 1 (zero) or 2 (fill)'
        if (.not. ieee_is_finite(resonance_cutout_width)) &
            error stop 'resonance_cutout_width must be finite'
        if (resonance_cutout_width < 1.0_dp) &
            error stop 'resonance_cutout_width must be at least 1'
    end subroutine validate_resonance_cutout

    pure function bounce_pitch_diffusion(nu_d, eta, etatp, etadt) result(D_eta)
        ! Orbit-averaged pitch diffusion D_eta = 2 nu_d eta <(1 - eta B)/B>
        ! in the large-aspect-ratio field B = B0 (1 - eps cos theta), where
        ! 1/etadt = Bmin = B0 (1 - eps) and 1/etatp = Bmax = B0 (1 + eps).
        ! With q = 2 eps eta B0 = eta (1/etatp - 1/etadt) and
        ! k2 = q / (1 - eta/etadt) (k2 = 1 on the trapped-passing boundary):
        !   trapped (k2 > 1), m = 1/k2:  <(1-eta B)/B> = (q/B0) (E(m)/K(m) - (1-m))
        !   passing (k2 < 1):            <(1-eta B)/B> = (q/(k2 B0)) E(k2)/K(k2).
        ! The trapped branch equals the D_eta_ba of layer_delta_eta. Both tend
        ! to zero logarithmically at the boundary.
        real(dp), intent(in) :: nu_d, eta, etatp, etadt
        real(dp) :: D_eta
        real(dp) :: q, k2, m, inv_B0, average, gap, ratio

        call cutout_pitch_parameters(eta, etatp, etadt, inv_B0, q, m, gap)
        if (eta >= etatp) then
            m = min(max(m, 0.0_dp), M_MAX)
            ratio = comp_ellint_e(m) / comp_ellint_k(m)
            call cutout_trapped_average(q, inv_B0, ratio, m, average)
        else
            call cutout_passing_modulus(q, gap, k2)
            k2 = min(max(k2, 0.0_dp), M_MAX)
            ratio = comp_ellint_e(k2) / comp_ellint_k(k2)
            call cutout_passing_average(gap, inv_B0, ratio, average)
        end if
        call cutout_pitch_diffusion(nu_d, eta, max(average, 0.0_dp), D_eta)
    end function bounce_pitch_diffusion

    pure function fold_layer_width(D_eta, abs_R_eta_eta) result(width)
        ! delta_fold = (D_eta / |R_eta_eta|)**(1/4); zero (no layer) unless
        ! both inputs are positive and finite.
        real(dp), intent(in) :: D_eta, abs_R_eta_eta
        real(dp) :: width

        width = 0.0_dp
        if (.not. (ieee_is_finite(D_eta) .and. ieee_is_finite(abs_R_eta_eta))) return
        if (D_eta > 0.0_dp .and. abs_R_eta_eta > 0.0_dp) &
            call cutout_fold_width(D_eta, abs_R_eta_eta, width)
    end function fold_layer_width

    pure subroutine cutout_weight(mode, width_factor, abs_R_eta, abs_R_eta_eta, &
                                  D_eta, weight, inside_layer)
        ! Evaluate the replacement weight before division by R_eta, so that
        ! the zero and fill variants remain finite at an exact fold.
        integer, intent(in) :: mode
        real(dp), intent(in) :: width_factor, abs_R_eta, abs_R_eta_eta, D_eta
        real(dp), intent(out) :: weight
        logical, intent(out) :: inside_layer
        real(dp) :: layer

        weight = 0.0_dp
        inside_layer = .false.
        if (mode /= CUTOUT_OFF) then
            layer = width_factor * fold_layer_width(D_eta, abs_R_eta_eta)
            if (layer > 0.0_dp) then
                if (abs_R_eta < abs_R_eta_eta * layer) then
                    select case (mode)
                    case (CUTOUT_ZERO)
                        inside_layer = .true.
                        weight = 0.0_dp
                        return
                    case (CUTOUT_FILL)
                        inside_layer = .true.
                        call cutout_boundary_fill(abs_R_eta_eta, layer, weight)
                        return
                    end select
                end if
            end if
        end if
    end subroutine cutout_weight

    include "generated/resonance_cutout_leaf.inc"

end module neort_resonance_cutout
