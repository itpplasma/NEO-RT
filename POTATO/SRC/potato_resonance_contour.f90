program potato_resonance_contour
    use parmot_mod, only : rmu, ro0
    use orbit_dim_mod, only : neqm, next, numbasef
    use global_invariants, only : dtau, cE_ref, Phi_eff
    use phielec_of_psi_mod, only : polyphi, polydens, polytemp
    use poicut_mod, only : npc, rpc_arr
    use potato_input_mod, only : read_potato_input, E_alpha, A_alpha, Z_alpha, &
        rho_pol_max, scalfac_energy, scalfac_efield, &
        Rmax_orbit, ntimstep, npoicut, profile_file, &
        edge_extension, contour_rho_min, contour_rho_max, contour_nrho, &
        contour_jperp_min, contour_jperp_max, contour_njperp, &
        contour_enkin, contour_sigma
    use field_eq_mod, only : allow_sol, psi_axis, psi_sep
    use field_sub, only : psif
    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    double precision :: v0

    external :: find_bounce, velo

    call read_potato_input("potato.in")
    call initialize_contour(v0)
    call write_contour(v0)

contains

    subroutine initialize_contour(v0)
        double precision, intent(out) :: v0
        double precision :: bmod_ref

        allow_sol = edge_extension
        E_alpha = E_alpha/scalfac_energy
        rmu = 1.d30
        v0 = sqrt(2.d0*E_alpha*ev/(p_mass*A_alpha))
        bmod_ref = 1.d0
        ro0 = v0*p_mass*A_alpha*c/(e_charge*Z_alpha*bmod_ref)
        cE_ref = E_alpha*ev
        Phi_eff = c*E_alpha*ev/(e_charge*Z_alpha*v0)
        dtau = Rmax_orbit/dble(ntimstep)
        numbasef = 0
        next = 0

        call load_profiles
        call find_poicut(rho_pol_max, npoicut)
    end subroutine initialize_contour

    subroutine write_contour(v0)
        double precision, intent(in) :: v0
        integer :: i, j, u
        double precision :: rho, jperp

        call validate_contour_input
        open(newunit=u, file="potato_resonance_contour.dat", &
            status="replace", action="write")
        write(u, '(A)') "# id rho_pol R_gc Z_gc p xi sigma H0 J_perp " &
            // "psi_star psi_axis psi_edge phi_elec v0_cm_s taub delphi " &
            // "omega_b omega_phi ierr"
        do i = 1, contour_nrho
            rho = grid_value(contour_rho_min, contour_rho_max, contour_nrho, i)
            do j = 1, contour_njperp
                jperp = grid_value(contour_jperp_min, contour_jperp_max, &
                    contour_njperp, j)
                call write_contour_point(u, (i - 1)*contour_njperp + j, &
                    rho, jperp, v0)
            end do
            write(u, *)
        end do
        close(u)
    end subroutine write_contour

    subroutine write_contour_point(u, id, rho, jperp, v0)
        integer, intent(in) :: u, id
        double precision, intent(in) :: rho, jperp, v0
        integer :: ierr
        double precision :: psi, phi_elec, dphi_dpsi, h0, rst, psi_star
        double precision :: dpsi_star_dr, taub, delphi, omega_b, omega_phi
        double precision :: z(neqm)

        psi = psi_axis + rho*rho*(psi_sep - psi_axis)
        call phielec_of_psi(psi, phi_elec, dphi_dpsi)
        h0 = contour_enkin + phi_elec
        rst = R_from_psi_lfs(psi)
        z = 0d0
        psi_star = 0d0
        dpsi_star_dr = 0d0
        call starter_doublecount(h0, jperp, dble(contour_sigma), rst, &
            psi_star, dpsi_star_dr, z, ierr)
        call trace_contour_point(z, taub, delphi, omega_b, omega_phi, v0, ierr)
        write(u, '(I8,17ES18.9,I8)') id, rho, z(1), z(3), z(4), z(5), &
            dble(contour_sigma), h0, jperp, psi_star, psi_axis, psi_sep, &
            phi_elec, v0, taub, delphi, omega_b, omega_phi, ierr
    end subroutine write_contour_point

    subroutine trace_contour_point(z, taub, delphi, omega_b, omega_phi, v0, ierr)
        double precision, intent(inout) :: z(neqm)
        double precision, intent(out) :: taub, delphi, omega_b, omega_phi
        double precision, intent(in) :: v0
        integer, intent(inout) :: ierr
        double precision :: extraset(1)

        taub = 0d0
        delphi = 0d0
        omega_b = 0d0
        omega_phi = 0d0
        if (ierr /= 0) return
        extraset = 0d0
        call find_bounce(next, velo, dtau, z, taub, delphi, extraset, ierr)
        if (ierr /= 0) return
        if (taub <= 0d0) then
            ierr = 1
            return
        endif
        omega_b = 2d0*acos(-1d0)*v0/taub
        omega_phi = delphi*v0/taub
    end subroutine trace_contour_point

    subroutine validate_contour_input
        if (contour_nrho < 1) error stop "contour_nrho must be positive"
        if (contour_njperp < 1) error stop "contour_njperp must be positive"
        if (contour_enkin <= 0d0) error stop "contour_enkin must be positive"
        if (abs(contour_sigma) /= 1) error stop "contour_sigma must be -1 or 1"
    end subroutine validate_contour_input

    double precision function grid_value(lower, upper, count, index)
        double precision, intent(in) :: lower, upper
        integer, intent(in) :: count, index

        if (count == 1) then
            grid_value = lower
        else
            grid_value = lower + (upper - lower)*dble(index - 1)/dble(count - 1)
        endif
    end function grid_value

    subroutine load_profiles
        integer :: iunit

        open(newunit=iunit, file=trim(profile_file), status="old", action="read")
        read(iunit, *)
        read(iunit, *)
        read(iunit, *) polydens
        read(iunit, *) polytemp
        read(iunit, *) polytemp
        read(iunit, *) polyphi
        close(iunit)
        polytemp = polytemp/scalfac_energy
        polyphi = polyphi/scalfac_efield
        polytemp = polytemp/E_alpha
        polyphi = polyphi*Z_alpha*e_charge/(E_alpha*ev)
    end subroutine load_profiles

    double precision function R_from_psi_lfs(psi_target)
        double precision, intent(in) :: psi_target
        integer :: k, best
        double precision :: Z, dZ_dR, x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        double precision :: err, best_err

        best = max(0, npc - 1)
        best_err = huge(1.d0)
        do k = 0, npc
            if (rpc_arr(k) < rpc_arr(npc/2)) cycle
            call get_poicut(rpc_arr(k), Z, dZ_dR)
            x = [rpc_arr(k), 0.d0, Z]
            call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            err = abs(psif - psi_target)
            if (err < best_err) then
                best = k
                best_err = err
            end if
        end do
        R_from_psi_lfs = rpc_arr(best)
    end function R_from_psi_lfs

end program potato_resonance_contour
