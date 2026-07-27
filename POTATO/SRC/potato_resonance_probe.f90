program potato_resonance_probe
    use parmot_mod, only : rmu, ro0
    use orbit_dim_mod, only : numbasef, &
        orbit_clip_resonance_classes => clip_resonance_classes
    use global_invariants, only : dtau, toten, perpinv, cE_ref, Phi_eff
    use bounds_fixpoints_mod, only : region_set_t
    use form_classes_doublecount_mod, only : nclasses, ifuntype, sigma_class, &
        R_class_beg, R_class_end
    use get_matrix_mod, only : iclass, delphi_max
    use resonance_mode_bounds_mod, only : resonant_delphi_bound
    use sample_matrix_mod, only : n1, npoi, xarr
    use phielec_of_psi_mod, only : polyphi, polydens, polytemp
    use potato_input_mod, only : read_potato_input, E_alpha, A_alpha, Z_alpha, &
        rho_pol, rho_pol_max, scalfac_energy, scalfac_efield, Rmax_orbit, &
        ntimstep, npoicut, profile_file, edge_extension, probe_rho_pol, &
        probe_ux, probe_eta, probe_toten, probe_perpinv, probe_m, probe_n, &
        input_clip_resonance_classes => clip_resonance_classes
    use field_eq_mod, only : allow_sol, psi_axis, psi_sep
    implicit none

    double precision, parameter :: c = 2.9979d10
    double precision, parameter :: e_charge = 4.8032d-10
    double precision, parameter :: p_mass = 1.6726d-24
    double precision, parameter :: ev = 1.6022d-12
    double precision, parameter :: pi = 3.14159265358979d0

    type(region_set_t) :: regions
    integer :: ierr, unit_out
    integer, dimension(1) :: probe_modes_m, probe_modes_n
    double precision :: psi, dens, temp, ddens, dtemp, phi_elec, dPhi_dpsi
    double precision :: enkin, v0, bmod_ref

    call read_potato_input("potato.in")
    orbit_clip_resonance_classes = input_clip_resonance_classes
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
    call rhopol_boundary(rho_pol)
    probe_modes_m = probe_m
    probe_modes_n = probe_n
    delphi_max = resonant_delphi_bound(probe_modes_m, probe_modes_n)

    psi = psi_axis + probe_rho_pol**2*(psi_sep - psi_axis)
    call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
    call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
    if (probe_toten >= 0.d0 .or. probe_perpinv >= 0.d0) then
        if (probe_toten < 0.d0 .or. probe_perpinv < 0.d0) then
            error stop 'probe_toten and probe_perpinv must be supplied together'
        endif
        toten = probe_toten
        perpinv = probe_perpinv
        enkin = toten - phi_elec
    else
        enkin = probe_ux**2*temp
        toten = enkin + phi_elec
        perpinv = probe_eta*enkin
    endif

    call find_bounds_fixpoints(regions, ierr)

    open(newunit=unit_out, file="potato_resonance_probe.dat", &
        status="replace", action="write")
    write(unit_out, '(A)') &
        "# target rho ux eta toten perpinv enkin nclasses ierr"

    if (ierr /= 0) then
        write(unit_out, '(7ES18.9,I8,I8)') probe_rho_pol, probe_ux, probe_eta, &
            toten, perpinv, enkin, 0.d0, 0, ierr
        close(unit_out)
        stop
    endif

    call write_regions(unit_out, regions)
    call form_classes_doublecount(regions, .false., ierr)
    write(unit_out, '(7ES18.9,I8,I8)') probe_rho_pol, probe_ux, probe_eta, &
        toten, perpinv, enkin, dble(nclasses), nclasses, ierr
    write(unit_out, '(A,I0,A,I0)') "# residual: delphi + 2*pi*m/n, m=", &
        probe_m, " n=", probe_n
    write(unit_out, '(A)') &
        "# iclass ifuntype sigma R_beg R_end x residual dresdx psiast taub delphi"

    do iclass = 1, nclasses
        call sample_class_doublecount(1, ierr)
        if (ierr /= 0) then
            write(unit_out, '(A,I0,A,I0)') "# class_error iclass=", iclass, &
                " ierr=", ierr
            call write_partial_class_grid(unit_out)
            cycle
        endif
        call write_class_probe(unit_out)
    enddo

    close(unit_out)

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

    subroutine write_class_probe(iunit)
        integer, intent(in) :: iunit
        integer :: i
        double precision :: residual, dresdx, psiast, taub, delphi
        double precision :: prev_x, prev_residual
        logical :: have_prev

        have_prev = .false.
        do i = 1, npoi
            call probe_residual(xarr(i), residual, dresdx, psiast, taub, delphi)
            write(iunit, '(I8,I8,9ES18.9)') iclass, ifuntype(iclass), &
                sigma_class(iclass), R_class_beg(iclass), R_class_end(iclass), &
                xarr(i), residual, dresdx, psiast, taub, delphi
            if (have_prev) then
                if (prev_residual*residual < 0.d0) then
                    write(iunit, '(A,I0,4ES18.9)') "# bracket ", iclass, &
                        prev_x, xarr(i), prev_residual, residual
                endif
            endif
            prev_x = xarr(i)
            prev_residual = residual
            have_prev = .true.
        enddo
    end subroutine write_class_probe

    subroutine write_partial_class_grid(iunit)
        use sample_matrix_mod, only : amat_arr

        integer, intent(in) :: iunit
        integer :: i

        write(iunit, '(A)') &
            "# partial iclass i x psiast taub delphi"
        do i = 1, npoi
            write(iunit, '(A,2I8,4ES24.15)') "# partial", iclass, i, xarr(i), &
                amat_arr(1:3,1,i)
        enddo
    end subroutine write_partial_class_grid

    subroutine write_regions(iunit, regions_in)
        integer, intent(in) :: iunit
        type(region_set_t), intent(in) :: regions_in
        integer :: isig, ireg, ix

        write(iunit, '(A,I0)') "# regions ", regions_in%nregions
        do ireg = 1, regions_in%nregions
            do isig = 1, 2
                write(iunit, '(A,2I6,L2,2I6,4ES18.9)') "# region ", isig, ireg, &
                    regions_in%all_regions(isig, ireg)%within_rhopol, &
                    regions_in%all_regions(isig, ireg)%n_o, &
                    regions_in%all_regions(isig, ireg)%n_x, &
                    regions_in%all_regions(isig, ireg)%R_b, &
                    regions_in%all_regions(isig, ireg)%R_e, &
                    regions_in%all_regions(isig, ireg)%psiast_b, &
                    regions_in%all_regions(isig, ireg)%psiast_e
                if (regions_in%all_regions(isig, ireg)%n_sep > 0) then
                    do ix = 1, regions_in%all_regions(isig, ireg)%n_sep
                        write(iunit, '(A,3I6,L2,ES18.9)') "# sep ", isig, ireg, &
                            ix, regions_in%all_regions(isig, ireg)%xpoint(ix), &
                            regions_in%all_regions(isig, ireg)%R_sep(ix)
                    enddo
                endif
            enddo
        enddo
    end subroutine write_regions

    subroutine probe_residual(x, residual, dresdx, psiast, taub, delphi)
        double precision, intent(in) :: x
        double precision, intent(out) :: residual, dresdx, psiast, taub, delphi
        double precision, dimension(n1) :: vec, dvec

        call interpolate_class_doublecount(x, vec, dvec)
        psiast = vec(1)
        taub = vec(2)
        delphi = vec(3)
        residual = delphi + 2.d0*pi*dble(probe_m)/dble(probe_n)
        dresdx = dvec(3)
    end subroutine probe_residual

end program potato_resonance_probe
