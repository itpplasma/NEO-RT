program potato_resonance_contour
    use parmot_mod, only : rmu, ro0
    use orbit_dim_mod, only : neqm, next, numbasef
    use global_invariants, only : dtau, cE_ref, Phi_eff
    use phielec_of_psi_mod, only : polyphi, polydens, polytemp
    use poicut_mod, only : npc, rpc_arr
    use potato_input_mod, only : read_potato_input, E_alpha, A_alpha, Z_alpha, &
        rho_pol_max, scalfac_energy, scalfac_efield, &
        Rmax_orbit, ntimstep, npoicut, profile_file, &
        edge_extension
    use field_eq_mod, only : allow_sol, psi_axis, psi_sep
    use field_sub, only : psif
    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    integer, parameter :: nrho = 70, neta = 90
    double precision, parameter :: ux = 1.5d0

    integer :: i, j, u, ierr
    double precision :: rho, eta, psi, dens, temp, ddens, dtemp, phi_elec, dPhi_dpsi
    double precision :: enkin, toten, perpinv, Rst, psiast, dpsiast_dRst
    double precision :: taub, delphi, extraset(1), z(neqm), v0, bmod_ref

    external :: find_bounce, velo

    call read_potato_input("potato.in")
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

    open(newunit=u, file="potato_resonance_contour.dat", status="replace", action="write")
    write(u, '(A)') &
        "# rho_pol eta ux delphi taub ierr psiast Rst toten perpinv R_gc Z_gc"
    do i = 1, nrho
        rho = 0.05d0 + 0.93d0*dble(i - 1)/dble(nrho - 1)
        psi = psi_axis + rho*rho*(psi_sep - psi_axis)
        call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
        call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
        enkin = ux*ux*temp
        toten = enkin + phi_elec
        Rst = R_from_psi_lfs(psi)
        do j = 1, neta
            eta = 1.d-6 + (8.d-5 - 1.d-6)*dble(j - 1)/dble(neta - 1)
            perpinv = eta*enkin
            call starter_doublecount(toten, perpinv, 1.d0, Rst, psiast, dpsiast_dRst, z, ierr)
            taub = 0.d0
            delphi = 0.d0
            if (ierr == 0) then
                extraset = 0.d0
                call find_bounce(next, velo, dtau, z, taub, delphi, extraset, ierr)
            end if
            write(u, '(12ES18.9)') rho, eta, ux, delphi, taub, dble(ierr), &
                psiast, Rst, toten, perpinv, z(1), z(3)
        end do
        write(u, *)
    end do
    close(u)

contains

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
