module neort_conserving_basis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_fp_tensor, only: neort_fp_q1_shape
    implicit none
    private

    type, public :: conserving_basis_t
        integer :: speed_steps = 0, pitch_steps = 0, npitch = 0
        integer :: ngraph = 0, ndof = 0, flow = 0, energy = 0
        real(dp) :: umax = 0.0_dp, etatp = 0.0_dp, etadt = 0.0_dp
        real(dp) :: vth_m_s = 0.0_dp
    end type conserving_basis_t
    public :: initialize_conserving_basis, evaluate_conserving_basis

contains

    subroutine initialize_conserving_basis(nu, np, umax, etatp, etadt, vth, &
                                           basis, info)
        integer, intent(in) :: nu, np
        real(dp), intent(in) :: umax, etatp, etadt, vth
        type(conserving_basis_t), intent(out) :: basis
        integer, intent(out) :: info

        info = -1
        if (min(nu, np) < 2) return
        if (.not. all(ieee_is_finite([umax, etatp, etadt, vth]))) return
        if (min(umax, etatp, vth) <= 0.0_dp .or. etadt <= etatp) return
        basis%speed_steps = nu
        basis%pitch_steps = np
        basis%npitch = 1 + 3*np
        basis%ngraph = 1 + nu*basis%npitch
        basis%ndof = basis%ngraph + 2
        basis%flow = basis%ngraph + 1
        basis%energy = basis%ngraph + 2
        basis%umax = umax
        basis%etatp = etatp
        basis%etadt = etadt
        basis%vth_m_s = vth
        info = 0
    end subroutine initialize_conserving_basis

    ! Gyrotropic basis about local z. Graph tails continue their outer trace
    ! constantly in speed, so no discontinuity is inserted at the speed wall.
    pure subroutine evaluate_conserving_basis(basis, bmod, velocity, values, &
                                              gradient, branch, dvalues, &
                                              htheta_cm, dlogb_dtheta)
        type(conserving_basis_t), intent(in) :: basis
        real(dp), intent(in) :: bmod, velocity(3)
        real(dp), intent(out) :: values(:), gradient(:, :)
        integer, intent(out), optional :: branch
        real(dp), intent(out), optional :: dvalues(:)
        real(dp), intent(in), optional :: htheta_cm, dlogb_dtheta
        real(dp) :: v, u, xi, eta, coordinate, span_eta, span_u, eu, ee
        real(dp) :: shape(4), du(4), de(4), grad(3), analytic_gradient(3)
        real(dp) :: flow, energy, flow_gradient, flow_dt, ht, db
        integer :: arm, cu, cp, offsets(3), indices(4), k

        values = 0.0_dp
        gradient = 0.0_dp
        if (present(dvalues)) dvalues = 0.0_dp
        v = sqrt(sum(velocity**2))
        if (v == 0.0_dp) then
            values(1) = 1.0_dp
            values(basis%energy) = -1.5_dp
            gradient(3, basis%flow) = 1.0_dp/basis%vth_m_s
            if (present(branch)) branch = 1
            return
        end if
        u = v/basis%vth_m_s
        xi = velocity(3)/v
        eta = max(0.0_dp, min(basis%etadt, (1.0_dp - xi**2)/bmod))
        if (eta > basis%etatp) then
            arm = 1
            span_eta = (basis%etadt - basis%etatp)/basis%pitch_steps
            coordinate = (eta - basis%etatp)/span_eta
        else
            arm = 2
            if (xi < 0.0_dp) arm = 3
            span_eta = -basis%etatp/basis%pitch_steps
            coordinate = (eta - basis%etatp)/span_eta
        end if
        if (present(branch)) branch = arm
        cp = max(1, min(basis%pitch_steps, int(coordinate) + 1))
        ee = 2.0_dp*(coordinate - real(cp - 1, dp)) - 1.0_dp
        span_u = basis%umax/basis%speed_steps
        cu = min(basis%speed_steps, int(u/span_u) + 1)
        eu = 2.0_dp*(min(u, basis%umax)/span_u - real(cu - 1, dp)) - 1
        call neort_fp_q1_shape(eu, ee, span_u, span_eta, shape, du, de)
        if (u >= basis%umax) du = 0.0_dp
        offsets = [1, 1 + basis%pitch_steps, 1 + 2*basis%pitch_steps]
        indices = [node(cu, cp), node(cu + 1, cp), &
                   node(cu + 1, cp + 1), node(cu, cp + 1)]
        ht = 0.0_dp
        db = 0.0_dp
        if (present(htheta_cm)) ht = htheta_cm
        if (present(dlogb_dtheta)) db = dlogb_dtheta
        do k = 1, 4
            call conserving_velocity_leaf(velocity(1), velocity(2), velocity(3), &
                 basis%vth_m_s, bmod, du(k), de(k), ht, db, grad(1), grad(2), grad(3), &
                             flow, energy, analytic_gradient(1), analytic_gradient(2), &
                                          analytic_gradient(3), flow_gradient, flow_dt)
            values(indices(k)) = values(indices(k)) + shape(k)
            gradient(:, indices(k)) = gradient(:, indices(k)) + grad
        end do
        values(basis%flow) = flow
        values(basis%energy) = energy
        gradient(3, basis%flow) = flow_gradient
        gradient(:, basis%energy) = analytic_gradient
        if (present(dvalues)) dvalues(basis%flow) = flow_dt

    contains

        pure integer function node(iu, ie) result(index)
            integer, intent(in) :: iu, ie
            integer :: pitch
            index = 1
            if (iu == 1) return
            pitch = 1
            if (ie > 1) pitch = offsets(arm) + ie - 1
            index = 1 + (iu - 2)*basis%npitch + pitch
        end function node
    end subroutine evaluate_conserving_basis
    include 'generated/conserving_basis_leaf.inc'
end module neort_conserving_basis
