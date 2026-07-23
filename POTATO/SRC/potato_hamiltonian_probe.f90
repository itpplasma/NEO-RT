program potato_hamiltonian_probe
    use parmot_mod, only: rmu, ro0
    use orbit_dim_mod, only: neqm, next, numbasef
    use global_invariants, only: dtau, cE_ref, Phi_eff
    use phielec_of_psi_mod, only: polyphi, polydens, polytemp
    use poicut_mod, only: npc, rpc_arr
    use potato_input_mod, only: read_potato_input, E_alpha, A_alpha, Z_alpha, &
        rho_pol_max, scalfac_energy, scalfac_efield, Rmax_orbit, ntimstep, &
        npoicut, profile_file, edge_extension, probe_rho_pol, probe_ux, &
        probe_eta, probe_m, probe_n
    use field_eq_mod, only: allow_sol, psi_axis, psi_sep
    use field_sub, only: psif
    use resint_mod, only: twopim2, rm3, taub_new, delphi_new, &
        toten_orb, perpinv_orb
    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    double precision, parameter :: pi = 3.14159265358979d0

    character(len=128) :: argument
    integer :: ierr, ierr_h
    double precision :: phi_shift, psi, dens, temp, ddens, dtemp
    double precision :: phi_elec, dPhi_dpsi, enkin, toten, perpinv
    double precision :: Rst, psiast, dpsiast_dRst, v0, bmod_ref
    double precision :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    double precision :: taub, delphi, absHn2, absHn2_production
    double precision :: extraset0(1), extraset3(3)
    double precision :: z(neqm), z_work(neqm), zext(neqm + 3), vzext(neqm + 3)
    complex(8) :: bmod_n, h_integral

    external :: find_bounce, velo, velo_res, pertham

    phi_shift = 0.d0
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) read(argument, *) phi_shift

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

    call load_profiles
    call find_poicut(rho_pol_max, npoicut)

    psi = psi_axis + probe_rho_pol**2*(psi_sep - psi_axis)
    call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
    call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
    enkin = probe_ux**2*temp
    toten = enkin + phi_elec
    perpinv = probe_eta*enkin
    Rst = R_from_psi_lfs(psi)
    call starter_doublecount(toten, perpinv, 1.d0, Rst, &
        psiast, dpsiast_dRst, z, ierr)
    if (ierr /= 0) error stop "starter_doublecount failed"
    z(2) = phi_shift

    call get_bmod_and_Phi(z(1:3), bmod, phi_elec)
    call magfie(z(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    call bmod_pert(z(1), z(3), bmod_n)

    toten_orb = z(4)**2 + phi_elec
    perpinv_orb = z(4)**2*(1.d0 - z(5)**2)/bmod
    twopim2 = 2.d0*pi*dble(probe_m)
    rm3 = dble(probe_n)

    next = 0
    extraset0 = 0.d0
    z_work = z
    call find_bounce(next, velo, dtau, z_work, taub, delphi, extraset0, ierr)
    if (ierr /= 0) error stop "unperturbed find_bounce failed"
    taub_new = taub
    delphi_new = delphi

    next = 3
    zext = 0.d0
    zext(1:neqm) = z
    call velo_res(0.d0, zext, vzext)
    extraset3 = 0.d0
    z_work = z
    call find_bounce(next, velo_res, dtau, z_work, taub, delphi, &
        extraset3, ierr_h)
    if (ierr_h /= 0) error stop "Hamiltonian find_bounce failed"
    h_integral = cmplx(extraset3(2)/taub, extraset3(3)/taub, kind=8)
    absHn2 = abs(h_integral)**2

    z_work = z
    twopim2 = 2.d0*pi*dble(probe_m)
    rm3 = dble(probe_n)
    call pertham(z_work, absHn2_production)

    write(*, '(A,2I8,12ES25.16)') "POTATO_HAMILTONIAN_STATE ", &
        probe_m, probe_n, phi_shift, z, psi, psiast, bmod, phi_elec, &
        taub_new, delphi_new
    write(*, '(A,8ES25.16)') "POTATO_HAMILTONIAN_FIELD ", &
        sqrtg, hctrvr, real(bmod_n), aimag(bmod_n), psif, psi_sep
    write(*, '(A,6ES25.16)') "POTATO_HAMILTONIAN_NATIVE ", &
        vzext(7), vzext(8), real(h_integral), aimag(h_integral), &
        absHn2, absHn2_production

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
        double precision :: height, dZ_dR, x(3), bmag, jacobian
        double precision :: bder_local(3), hcov_local(3)
        double precision :: hcon_local(3), hcurl_local(3)
        double precision :: error, best_error

        best = max(0, npc - 1)
        best_error = huge(1.d0)
        do k = 0, npc
            if (rpc_arr(k) < rpc_arr(npc/2)) cycle
            call get_poicut(rpc_arr(k), height, dZ_dR)
            x = [rpc_arr(k), 0.d0, height]
            call magfie(x, bmag, jacobian, bder_local, hcov_local, &
                hcon_local, hcurl_local)
            error = abs(psif - psi_target)
            if (error < best_error) then
                best = k
                best_error = error
            end if
        end do
        R_from_psi_lfs = rpc_arr(best)
    end function R_from_psi_lfs

end program potato_hamiltonian_probe
