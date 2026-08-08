program test_gc_cut_topology_primitives
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cut_topology_flow_jacobian_symbolic, only: &
        evaluate_neort_gc_cut_topology_flow_jacobian
    use neort_gc_cut_topology_partner_symbolic, only: &
        evaluate_neort_gc_cut_topology_partner
    use neort_gc_cut_topology_section_flow_symbolic, only: &
        evaluate_neort_gc_cut_topology_section_flow
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: tolerance = 1.0e-13_dp
    real(dp) :: flow, trace, determinant, discriminant
    real(dp) :: partner_residual, partner_derivative
    real(dp) :: radius, height, partner_radius, partner_height
    real(dp) :: p_star, p_star_x, dp_star_d_r, dp_star_d_z, z_r
    logical :: classification_unresolved

    call evaluate_neort_gc_cut_topology_section_flow(2.0_dp, 5.0_dp, 1.5_dp, flow)
    if (abs(flow - 2.0_dp) > tolerance) then
        error stop "normal section flow has the wrong graph-slope contraction"
    end if

    call evaluate_neort_gc_cut_topology_flow_jacobian(0.0_dp, -1.0_dp, 1.0_dp, &
        0.0_dp, trace, determinant, discriminant)
    if (abs(discriminant + 4.0_dp) > tolerance) then
        error stop "elliptic O-flow discriminant is wrong"
    end if
    if (discriminant >= 0.0_dp) error stop "O-flow was not classified elliptic"

    call evaluate_neort_gc_cut_topology_flow_jacobian(1.0_dp, 0.0_dp, 0.0_dp, &
        -1.0_dp, trace, determinant, discriminant)
    if (abs(discriminant - 4.0_dp) > tolerance) then
        error stop "hyperbolic X-flow discriminant is wrong"
    end if
    if (discriminant <= 0.0_dp) error stop "X-flow was not classified hyperbolic"

    call evaluate_neort_gc_cut_topology_flow_jacobian(0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, trace, determinant, discriminant)
    classification_unresolved = abs(discriminant) <= tolerance
    if (.not. classification_unresolved) then
        error stop "degenerate flow was not left unresolved"
    end if

    radius = 1.0_dp
    height = 2.0_dp
    partner_radius = 2.0_dp
    partner_height = -1.0_dp
    p_star = radius**2 + height
    p_star_x = partner_radius**2 + partner_height
    dp_star_d_r = 2.0_dp*radius
    dp_star_d_z = 1.0_dp
    z_r = 0.5_dp
    call evaluate_neort_gc_cut_topology_partner(p_star, p_star_x, dp_star_d_r, &
        dp_star_d_z, z_r, partner_residual, partner_derivative)
    if (abs(p_star - p_star_x) > tolerance) then
        error stop "manufactured partner points do not share p_star"
    end if
    if (abs(partner_residual) > tolerance) then
        error stop "same-invariant partner residual is wrong"
    end if
    if (abs(partner_derivative - (2.0_dp*radius + z_r)) > tolerance) then
        error stop "graph-cut partner derivative is wrong"
    end if

    call pass_test
end program test_gc_cut_topology_primitives
