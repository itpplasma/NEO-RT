module neort_gc_axisymmetric_physical_flow_classifier
    !! Runtime boundary for the generated physical fixed-point Jacobian.
    !!
    !! The caller owns the stationary-point certificate.  This module owns
    !! only physical-domain validation, the singularity gate, invocation of
    !! the Fortsym kernel, and the sign-preserving O/X classification.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_axisymmetric_physical_flow_domain_symbolic, only: &
        evaluate_neort_gc_axisymmetric_physical_flow_domain
    use neort_gc_axisymmetric_physical_flow_jacobian_symbolic, only: &
        evaluate_neort_gc_axisymmetric_physical_flow_jacobian
    implicit none
    private

    integer, parameter, public :: GC_PHYSICAL_FLOW_SUCCESS = 0
    integer, parameter, public :: GC_PHYSICAL_FLOW_INVALID_INPUT = 1
    integer, parameter, public :: GC_PHYSICAL_FLOW_MISSING_STATIONARY_CERTIFICATE = 2
    integer, parameter, public :: GC_PHYSICAL_FLOW_OUT_OF_DOMAIN = 3
    integer, parameter, public :: GC_PHYSICAL_FLOW_TURNING = 4
    integer, parameter, public :: GC_PHYSICAL_FLOW_SINGULAR_DENOMINATOR = 5
    integer, parameter, public :: GC_PHYSICAL_FLOW_NONFINITE = 6
    integer, parameter, public :: GC_PHYSICAL_FLOW_DEGENERATE = 7

    integer, parameter, public :: GC_PHYSICAL_FLOW_UNRESOLVED = 0
    integer, parameter, public :: GC_PHYSICAL_FLOW_O = 1
    integer, parameter, public :: GC_PHYSICAL_FLOW_X = 2

    type, public :: gc_axisymmetric_physical_flow_input_t
        real(dp) :: mass = 0.0_dp
        real(dp) :: charge = 0.0_dp
        real(dp) :: c_light = 0.0_dp
        real(dp) :: radius = 0.0_dp
        real(dp) :: h0 = 0.0_dp
        real(dp) :: j_k = 0.0_dp
        real(dp) :: sigma = 0.0_dp
        real(dp) :: psi = 0.0_dp
        real(dp) :: psi_r = 0.0_dp
        real(dp) :: psi_z = 0.0_dp
        real(dp) :: psi_rr = 0.0_dp
        real(dp) :: psi_rz = 0.0_dp
        real(dp) :: psi_zz = 0.0_dp
        real(dp) :: psi_rrr = 0.0_dp
        real(dp) :: psi_rrz = 0.0_dp
        real(dp) :: psi_rzz = 0.0_dp
        real(dp) :: psi_zzz = 0.0_dp
        real(dp) :: psi_sep = 0.0_dp
        real(dp) :: f = 0.0_dp
        real(dp) :: df_dpsihat = 0.0_dp
        real(dp) :: d2f_dpsihat2 = 0.0_dp
        real(dp) :: phi = 0.0_dp
        real(dp) :: dphi_dpsi = 0.0_dp
        real(dp) :: d2phi_dpsi2 = 0.0_dp
    end type gc_axisymmetric_physical_flow_input_t

    type, public :: gc_axisymmetric_stationary_point_certificate_t
        logical :: certified = .false.
        integer :: certificate_id = 0
        real(dp) :: residual_r = 0.0_dp
        real(dp) :: residual_z = 0.0_dp
        real(dp) :: tolerance = 0.0_dp
    end type gc_axisymmetric_stationary_point_certificate_t

    type, public :: gc_axisymmetric_physical_flow_result_t
        integer :: status = GC_PHYSICAL_FLOW_INVALID_INPUT
        integer :: classification = GC_PHYSICAL_FLOW_UNRESOLVED
        integer :: stationary_certificate_id = 0
        logical :: stationary_certificate_accepted = .false.
        real(dp) :: bmod = 0.0_dp
        real(dp) :: bparallel_star = 0.0_dp
        real(dp) :: p_parallel = 0.0_dp
        real(dp) :: energy_margin = 0.0_dp
        real(dp) :: j11 = 0.0_dp
        real(dp) :: j12 = 0.0_dp
        real(dp) :: j21 = 0.0_dp
        real(dp) :: j22 = 0.0_dp
        real(dp) :: trace = 0.0_dp
        real(dp) :: determinant = 0.0_dp
        real(dp) :: discriminant = 0.0_dp
    end type gc_axisymmetric_physical_flow_result_t

    public :: classify_gc_axisymmetric_physical_fixed_point

contains

    subroutine classify_gc_axisymmetric_physical_fixed_point(input, &
            stationary, result, status)
        type(gc_axisymmetric_physical_flow_input_t), intent(in) :: input
        type(gc_axisymmetric_stationary_point_certificate_t), intent(in) :: &
            stationary
        type(gc_axisymmetric_physical_flow_result_t), intent(out) :: result
        integer, intent(out) :: status

        real(dp) :: bmod, bparallel_star, p_parallel, energy_margin
        real(dp) :: p_parallel_squared
        real(dp) :: j11, j12, j21, j22, trace, determinant, discriminant

        result = gc_axisymmetric_physical_flow_result_t()
        status = GC_PHYSICAL_FLOW_INVALID_INPUT

        if (.not. finite_input(input)) return
        if (.not. valid_stationary_certificate(stationary)) then
            status = GC_PHYSICAL_FLOW_MISSING_STATIONARY_CERTIFICATE
            result%status = status
            return
        end if
        result%stationary_certificate_accepted = .true.
        result%stationary_certificate_id = stationary%certificate_id

        if (input%mass <= 0.0_dp .or. input%c_light <= 0.0_dp .or. &
                input%radius <= 0.0_dp .or. input%psi_sep <= 0.0_dp .or. &
                input%j_k < 0.0_dp .or. &
                singular_denominator(input%charge, 1.0_dp) .or. &
                abs(input%sigma) /= 1.0_dp) then
            status = GC_PHYSICAL_FLOW_OUT_OF_DOMAIN
            result%status = status
            return
        end if

        call evaluate_neort_gc_axisymmetric_physical_flow_domain(input%mass, &
            input%charge, input%c_light, input%radius, input%h0, input%j_k, &
            input%psi, input%psi_r, input%psi_z, input%psi_sep, input%f, &
            input%phi, bmod, energy_margin, p_parallel_squared)
        if (.not. all(ieee_is_finite([bmod, energy_margin, &
                p_parallel_squared]))) then
            status = GC_PHYSICAL_FLOW_NONFINITE
            result%status = status
            return
        end if
        result%bmod = bmod
        result%energy_margin = energy_margin
        if (bmod <= 0.0_dp) then
            status = GC_PHYSICAL_FLOW_SINGULAR_DENOMINATOR
            result%status = status
            return
        end if
        if (energy_margin <= 0.0_dp .or. p_parallel_squared <= 0.0_dp) then
            status = GC_PHYSICAL_FLOW_TURNING
            result%status = status
            return
        end if

        call evaluate_neort_gc_axisymmetric_physical_flow_jacobian( &
            input%mass, input%charge, input%c_light, input%radius, input%h0, &
            input%j_k, input%sigma, input%psi, input%psi_r, input%psi_z, &
            input%psi_rr, input%psi_rz, input%psi_zz, input%psi_rrr, &
            input%psi_rrz, input%psi_rzz, input%psi_zzz, input%psi_sep, &
            input%f, input%df_dpsihat, input%d2f_dpsihat2, input%phi, &
            input%dphi_dpsi, input%d2phi_dpsi2, j11, j12, j21, j22, trace, &
            determinant, discriminant, bmod, p_parallel, bparallel_star)

        if (.not. all(ieee_is_finite([j11, j12, j21, j22, trace, &
                determinant, discriminant, bmod, p_parallel, &
                bparallel_star]))) then
            status = GC_PHYSICAL_FLOW_NONFINITE
            result%status = status
            return
        end if

        if (singular_denominator(p_parallel, sqrt(p_parallel_squared))) then
            status = GC_PHYSICAL_FLOW_TURNING
            result%status = status
            return
        end if
        if (singular_denominator(bparallel_star, bmod)) then
            status = GC_PHYSICAL_FLOW_SINGULAR_DENOMINATOR
            result%status = status
            return
        end if
        result%bmod = bmod
        result%bparallel_star = bparallel_star
        result%p_parallel = p_parallel
        result%energy_margin = energy_margin
        result%j11 = j11
        result%j12 = j12
        result%j21 = j21
        result%j22 = j22
        result%trace = trace
        result%determinant = determinant
        result%discriminant = discriminant

        if (discriminant < 0.0_dp) then
            result%classification = GC_PHYSICAL_FLOW_O
            status = GC_PHYSICAL_FLOW_SUCCESS
        else if (discriminant > 0.0_dp) then
            result%classification = GC_PHYSICAL_FLOW_X
            status = GC_PHYSICAL_FLOW_SUCCESS
        else
            result%classification = GC_PHYSICAL_FLOW_UNRESOLVED
            status = GC_PHYSICAL_FLOW_DEGENERATE
        end if
        result%status = status
    end subroutine classify_gc_axisymmetric_physical_fixed_point

    logical function finite_input(input)
        type(gc_axisymmetric_physical_flow_input_t), intent(in) :: input

        finite_input = all(ieee_is_finite([input%mass, input%charge, &
            input%c_light, input%radius, input%h0, input%j_k, input%sigma, &
            input%psi, input%psi_r, input%psi_z, input%psi_rr, input%psi_rz, &
            input%psi_zz, input%psi_rrr, input%psi_rrz, input%psi_rzz, &
            input%psi_zzz, input%psi_sep, input%f, input%df_dpsihat, &
            input%d2f_dpsihat2, input%phi, input%dphi_dpsi, &
            input%d2phi_dpsi2]))
    end function finite_input

    logical function valid_stationary_certificate(certificate)
        type(gc_axisymmetric_stationary_point_certificate_t), intent(in) :: &
            certificate

        valid_stationary_certificate = .false.
        if (.not. certificate%certified) return
        if (certificate%certificate_id <= 0) return
        if (.not. all(ieee_is_finite([certificate%residual_r, &
                certificate%residual_z, certificate%tolerance]))) return
        if (certificate%tolerance <= 0.0_dp) return
        if (abs(certificate%residual_r) > certificate%tolerance) return
        if (abs(certificate%residual_z) > certificate%tolerance) return
        valid_stationary_certificate = .true.
    end function valid_stationary_certificate

    logical function singular_denominator(value, scale)
        real(dp), intent(in) :: value, scale
        real(dp) :: comparison_scale

        singular_denominator = .true.
        if (.not. ieee_is_finite(value)) return
        if (.not. ieee_is_finite(scale)) return
        comparison_scale = max(1.0_dp, abs(scale))
        singular_denominator = abs(value) <= 64.0_dp*epsilon(1.0_dp)* &
            comparison_scale
    end function singular_denominator

end module neort_gc_axisymmetric_physical_flow_classifier
