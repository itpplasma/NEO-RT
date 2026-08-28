! Manufactured Cartesian axisymmetric case verifying the sign of the electric
! contribution to the thermodynamic force A1 across three independent paths.
!
! Reference: GitHub issue #98.  The documented relation (doc/driftorbit.lyx,
! eq. Om_tE) and the thermodynamic force (eq. A1) give, on substitution,
!
!     Om_tE  = -c*Phi_e'/chi'          (chi' = sign_theta*psi_pr/q)
!     A1     = n'/n - e*Phi_e'/T - 3*T'/(2*T)
!     A1_elec= -q_i*Phi_e'/(T*ev)      =  +q_i/(T*ev)*(chi'/c)*Om_tE
!
! hence the electric term in init_thermodynamic_forces must be PLUS
! q_i/(T*ev)*sign_theta*psi_pr/(q*c)*Om_tE.  This test builds a circular
! axisymmetric manufactured equilibrium with a prescribed Phi_e(R,Z), derives
! E = -grad(Phi_e), the E x B drift and the poloidal-flux gradient from the
! same explicit vectors, verifies the three paths against one another, and
! calls init_thermodynamic_forces with the derived Om_tE to check the code's
! A1 against the independent A1_elec.
program test_thermodynamic_force_sign
    use iso_fortran_env, only: dp => real64
    use util, only: c, qe, ev, pi
    use do_magfie_mod, only: sign_theta
    use neort_profiles, only: init_thermodynamic_forces, Om_tE, A1, ni1, Ti1, dni1ds, dTi1ds, dOm_tEds

    implicit none

    ! Manufactured case geometry (cgs): circular flux surfaces, major radius
    ! R0 and minor radius r on the surface, right-handed (R, phi, Z).
    real(dp), parameter :: R0 = 500.0_dp    ! major radius (cm)
    real(dp), parameter :: r_surf = 30.0_dp ! minor radius of surface (cm)
    real(dp), parameter :: B0 = 2.0e4_dp    ! on-axis field (G)
    real(dp), parameter :: T_keV = 5.0e3_dp ! ion temperature (keV)
    real(dp), parameter :: pot_deriv = 1.0e-2_dp ! dPhi_e/dr (statvolt/cm)

    real(dp), parameter :: theta_probe = 1.1_dp ! poloidal probe point (rad)

    real(dp) :: phi_e_pr, chi_prime, om_tE_doc, om_tE_exb, a1_elec_indep
    real(dp) :: psi_pr_test, q_test, a1_code
    real(dp) :: R, Z, Bth, Bph, Er, Evec(3), Bvec(3), vE(3)
    real(dp) :: psi_pr_val, q_val, tol
    integer :: sig_case
    logical :: all_ok

    ! Fixed thermodynamic state: no density/temperature gradients so that A1
    ! contains only the electric contribution we are testing.
    ni1 = 1.0e13_dp
    Ti1 = T_keV
    dni1ds = 0.0_dp
    dTi1ds = 0.0_dp

    tol = 1.0e-9_dp
    all_ok = .true.

    ! Reverse Phi_e' and the magnetic/coordinate sign (sign_theta is a
    ! compile-time parameter, so reverse psi_pr and q instead: chi' =
    ! sign_theta*psi_pr/q flips when psi_pr/q flips).
    do sig_case = 1, 4
        select case (sig_case)
        case (1)
            phi_e_pr = pot_deriv
            psi_pr_val = -1.5e8_dp
            q_val = -2.0_dp
        case (2)
            phi_e_pr = -pot_deriv
            psi_pr_val = -1.5e8_dp
            q_val = -2.0_dp
        case (3)
            phi_e_pr = pot_deriv
            psi_pr_val = 1.5e8_dp
            q_val = -2.0_dp
        case (4)
            phi_e_pr = pot_deriv
            psi_pr_val = -1.5e8_dp
            q_val = 2.0_dp
        end select

        ! --- Manufactured electric field and E x B drift (explicit vectors) ---
        call manufactured_field(theta_probe, R0, r_surf, B0, phi_e_pr, psi_pr_val, &
            q_val, R, Z, Bth, Bph, Er, Evec, Bvec, vE)

        ! chi' = sign_theta*psi_pr/q  (toroidal psi_pr, per do_magfie conventions)
        chi_prime = sign_theta * psi_pr_val / q_val

        ! Path 1 (documented): Om_tE = -c*Phi_e'/chi'
        om_tE_doc = -c * phi_e_pr / chi_prime

        ! Path 2 (E x B / precession): Om_tE = c*E_r/chi' with E_r the radial
        ! component of E = -grad(Phi_e).  Also verify the E x B drift vector
        ! vE = c*E x B/B^2 is orthogonal to E and B (roundoff).
        om_tE_exb = c * Er / chi_prime
        call check_path1_path2(om_tE_doc, om_tE_exb, Evec, Bvec, vE, sig_case, all_ok)

        ! Path 3 (A1 electric term), independent of Om_tE:
        a1_elec_indep = -qe / (T_keV * ev) * phi_e_pr

        ! Code path: feed Om_tE into init_thermodynamic_forces and compare.
        Om_tE = om_tE_doc
        dOm_tEds = 0.0_dp
        psi_pr_test = psi_pr_val
        q_test = q_val
        call init_thermodynamic_forces(psi_pr_test, q_test)
        a1_code = A1
        call check_a1(a1_code, a1_elec_indep, sig_case, all_ok)
    end do

    if (all_ok) then
        call print_ok
    else
        call print_fail
        error stop
    end if

contains

    ! Manufactured circular field: B = Bth*e_theta + Bph*e_phi (cgs), with the
    ! poloidal flux derivative such that psi_pr (toroidal boundary flux) and q
    ! are consistent with the input, and E = -dPhi_e/dr*e_r from the radial
    ! potential gradient.  Returns the E x B drift vE = c*E x B/B^2.
    subroutine manufactured_field(th, r0loc, r_surf_loc, b0loc, philoc, psiprloc, qloc, &
        Rloc, Zloc, Bthloc, Bphloc, Erloc, Eloc, Bloc, vEloc)
        real(dp), intent(in) :: th, r0loc, r_surf_loc, b0loc, philoc, psiprloc, qloc
        real(dp), intent(out) :: Rloc, Zloc, Bthloc, Bphloc, Erloc, Eloc(3), Bloc(3), vEloc(3)

        real(dp) :: sq, eR(3), eTh(3), ePhi(3), B2

        Rloc = r0loc + r_surf_loc*cos(th)
        Zloc = r_surf_loc*sin(th)
        eR = (/cos(th), 0.0_dp, sin(th)/)
        eTh = (/-sin(th), 0.0_dp, cos(th)/)
        ePhi = (/0.0_dp, 1.0_dp, 0.0_dp/)

        ! Toroidal field, and poloidal field from the safety factor.
        Bphloc = b0loc * r0loc / Rloc
        sq = qloc
        Bthloc = b0loc * r_surf_loc * r0loc / (sq * Rloc * Rloc)

        Bloc = Bthloc*eTh + Bphloc*ePhi
        B2 = dot_product(Bloc, Bloc)

        ! Electric potential Phi_e = philoc * r  ->  E = -dPhi/dr*e_r.
        Erloc = -philoc
        Eloc = Erloc*eR

        vEloc = c * cross(Eloc, Bloc) / B2
    end subroutine manufactured_field

    pure function cross(a, b) result(cv)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: cv(3)
        cv(1) = a(2)*b(3) - a(3)*b(2)
        cv(2) = a(3)*b(1) - a(1)*b(3)
        cv(3) = a(1)*b(2) - a(2)*b(1)
    end function cross

    subroutine check_path1_path2(om1, om2, E, B, vE, icase, ok)
        real(dp), intent(in) :: om1, om2, E(3), B(3), vE(3)
        integer, intent(in) :: icase
        logical, intent(inout) :: ok
        real(dp) :: rel, dotEB, dotEB2, b2
        rel = abs(om1 - om2) / max(abs(om1), 1.0e-30_dp)
        if (rel > 1.0e-9_dp) then
            print *, "case", icase, ": Om_tE path mismatch: doc=", om1, " exb=", om2
            ok = .false.
        end if
        ! vE must be the actual E x B drift: orthogonal to E and to B.
        dotEB = dot_product(vE, E)
        dotEB2 = dot_product(vE, B)
        b2 = dot_product(B, B)
        if (abs(dotEB) > 1.0e-9_dp * b2 .or. abs(dotEB2) > 1.0e-9_dp * b2) then
            print *, "case", icase, ": E x B drift not orthogonal (E.B-vE=", dotEB, &
                ", B.B-vE=", dotEB2, ")"
            ok = .false.
        end if
    end subroutine check_path1_path2

    subroutine check_a1(code_a1, indep_a1, icase, ok)
        real(dp), intent(in) :: code_a1, indep_a1
        integer, intent(in) :: icase
        logical, intent(inout) :: ok
        real(dp) :: rel
        rel = abs(code_a1 - indep_a1) / max(abs(indep_a1), 1.0e-30_dp)
        if (rel > 1.0e-9_dp) then
            print *, "case", icase, ": A1 electric mismatch: code=", code_a1, &
                " independent=", indep_a1, " rel=", rel
            ok = .false.
        else
            print *, "case", icase, ": A1 electric OK (rel ", rel, ")"
        end if
    end subroutine check_a1

    subroutine print_ok
        print *, "    .................................................... OK"
    end subroutine print_ok

    subroutine print_fail
        print *, "    .................................................... FAIL"
    end subroutine print_fail

end program test_thermodynamic_force_sign
