!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! modules
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module resint_mod
    use resonance_status_mod, only : resonance_status_success
    !
    type respoints_fix_jperp_mode_class
    integer :: nrespoi
    integer :: scan_status=resonance_status_success
    double precision :: toten_res,perpinv_res
    double precision, dimension(:),   allocatable :: w_res
    double precision, dimension(:,:), allocatable :: z_res
    double precision, dimension(:),   allocatable :: taub
end type
!
type respoint_single
double precision :: toten_res,perpinv_res,w_res
double precision, dimension(5) :: z_res
double precision :: taub
end type
!
type respoints_fix_jperp
integer :: nclasses
type(respoints_fix_jperp_mode_class), dimension(:,:), allocatable :: respoints_jp
end type
!
integer          :: nmodes,nperp_max=0
double precision :: twopim2,rm3,taub_new,delphi_new
!$omp threadprivate(twopim2,rm3,taub_new,delphi_new)
! Energy and perpendicular invariant of the resonant orbit, set by pertham and
! read by velo_res for the perturbed-Hamiltonian Fourier integral.  Held separate
! from the global class invariants toten,perpinv so the per-mode resonance loop
! never writes those -- they stay the shared, read-only class values, which lets
! the loop be parallelized without making toten,perpinv threadprivate.
double precision :: toten_orb,perpinv_orb
!$omp threadprivate(toten_orb,perpinv_orb)
integer :: perturbation_status=0
!$omp threadprivate(perturbation_status)
! By-products of the resonance condition get_rescond ($\psi^\ast$, bounce time,
! toroidal shift), recovered for the orbit weight.  get_rescond is an internal
! procedure of integrate_class_resonances that is also passed as the dummy root
! function to find_all_roots; with gfortran's trampoline for that dummy-procedure
! call, a host-local PRIVATE variable written through it is not seen as private,
! so these must be threadprivate module state instead of privatized host locals.
double precision :: psiast_res,taub_res,delphi_res
!$omp threadprivate(psiast_res,taub_res,delphi_res)
integer, dimension(:), allocatable :: marr,narr
double precision, dimension(:), allocatable :: delint_mode
!
type(respoints_fix_jperp_mode_class), dimension(:,:), allocatable :: respoints_jp
type(respoints_fix_jperp),            dimension(:),   allocatable :: respoints_all, &
    respoints_all_tmp
type(respoint_single),                dimension(:),   allocatable :: respoint
integer :: resline_unit=31415,resline_diag_unit=31416
logical :: resline_unit_is_private=.false.,resline_diag_unit_is_private=.false.
! Per-energy conservation ledger.  A searched-zero harmonic is counted
! separately from a failed search; neither is reconstructed from nrespoi alone.
integer :: ledger_class_evaluations=0,ledger_harmonic_scans=0
integer :: ledger_searched_zero=0,ledger_root_count=0,ledger_failure_count=0
integer :: ledger_jperp_samples=0
! Per-energy-slice resonance bookkeeping.  The energy loop may run slices in
! parallel; each slice owns its J_perp grid, extracted resonant points, and
! temporary resonance-line unit.  nmodes,marr,narr are read-only after setup.
!$omp threadprivate(nperp_max,delint_mode,respoints_jp,respoints_all, &
!$omp               respoints_all_tmp,respoint,resline_unit,resline_diag_unit, &
!$omp               resline_unit_is_private,resline_diag_unit_is_private)
!$omp threadprivate(ledger_class_evaluations,ledger_harmonic_scans, &
!$omp               ledger_searched_zero,ledger_root_count,ledger_failure_count, &
!$omp               ledger_jperp_samples)

end module resint_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! rotines
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine velo_res(dtau,z,vz)
    !
    ! Computes the extended RHS for equations of motion and Fourier amplitude
    ! of the perturbed Hamiltonian $\hat H_\bm$, Eq.(102) (former Eq.(93))
    !
    use orbit_dim_mod,     only : neqm,next,numbasef
    use resint_mod,        only : twopim2,rm3,taub_new,delphi_new,toten_orb,perpinv_orb, &
                                  perturbation_status
    use bmod_pert_mod,     only : bmod_pert_success
    use potato_symbolic_kernel_mod, only : potato_hm_kernel
    !
    implicit none
    !
    double precision :: dtau,bmod,phi_elec,mode_phase
    double precision :: hm_real,hm_imag,hm_squared,hamiltonian_coefficient
    integer :: perturbation_ierr
    double precision, dimension(neqm+next) :: z,vz
    complex(8) :: bmod_n
    external :: bmod_pert_status
    !
    call velo(dtau,z(1:neqm),vz(1:neqm))
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
    call bmod_pert_status(z(1),z(3),bmod_n,perturbation_ierr)
    if(perturbation_ierr.ne.bmod_pert_success) then
        perturbation_status=perturbation_ierr
        vz(7)=0.d0
        vz(8)=0.d0
        return
    endif
    !
    mode_phase=rm3*z(2)-(twopim2+delphi_new*rm3)*z(6)/taub_new
    call potato_hm_kernel(toten_orb,phi_elec,bmod,perpinv_orb, &
        real(bmod_n),aimag(bmod_n),mode_phase,hm_real,hm_imag,hm_squared, &
        hamiltonian_coefficient)
    !
    vz(6)=1.d0
    vz(7)=hm_real
    vz(8)=hm_imag
    !
end subroutine velo_res
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine pertham(z,absHn2,ierr_status)
    !
    ! Computes modulus squared of the Fourier amplitude of the normalized perturbed
    ! Hamiltoninan, $|\hat H_\bm|^2$, with $\hat H_\bm$ defined by Eq.(102) (former Eq.(93)).
    !
    use orbit_dim_mod,     only : neqm,next     ,write_orb,iunit1
    use global_invariants, only : dtau
    use resint_mod,        only : taub_new,delphi_new,toten_orb,perpinv_orb, &
                                  perturbation_status
    use bmod_pert_mod,      only : bmod_pert_success
    use resonance_status_mod, only : resonance_status_success, resonance_is_finite
    !
    !
    implicit none
    !
    complex(8), parameter :: imun=(0.d0,1.d0)
    !
    integer :: ierr
    integer, intent(out) :: ierr_status
    double precision :: absHn2,bmod,phi_elec,taub,delphi
    double precision, dimension(neqm) :: z
    double precision, dimension(:), allocatable :: extraset
    !
    external :: velo,velo_res
    !
    ierr_status=resonance_status_success
    absHn2=0.d0
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
    !
    toten_orb=z(4)**2+phi_elec
    perpinv_orb=z(4)**2*(1.d0-z(5)**2)/bmod
    !
    next=0
    allocate(extraset(next))
    !
    call find_bounce(next,velo,dtau,z,taub,delphi,extraset,ierr)
    if(ierr.ne.0) then
        ! Preserve POTATO's established zero-amplitude handling for an
        ! unintegrable Fourier sample.  Topology and root failures remain
        ! strict; this is only the perturbation weight at one sample.
        deallocate(extraset)
        return
    endif
    if(.not.resonance_is_finite(taub) .or. taub.le.0.d0 .or. &
       .not.resonance_is_finite(delphi)) then
        deallocate(extraset)
        return
    endif
    !
    taub_new=taub
    delphi_new=delphi
    deallocate(extraset)
    !
    next=3
    allocate(extraset(next))
    extraset=0.d0
    perturbation_status=bmod_pert_success
    !
    call find_bounce(next,velo_res,dtau,z,taub,delphi,extraset,ierr)
    if(perturbation_status.ne.bmod_pert_success) then
        deallocate(extraset)
        return
    endif
    if(ierr.ne.0) then
        deallocate(extraset)
        return
    endif
    if(.not.resonance_is_finite(taub) .or. taub.le.0.d0 .or. &
       .not.resonance_is_finite(delphi)) then
        deallocate(extraset)
        return
    endif
    !
    absHn2=(extraset(2)/taub)**2+(extraset(3)/taub)**2
    if(.not.resonance_is_finite(absHn2) .or. absHn2.lt.0.d0) then
        deallocate(extraset)
        return
    endif
    !
    deallocate(extraset)
    !
end subroutine pertham
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine integrate_class_resonances(ierr_out,append)
    !
    ! Computes sum over resonances $x=x^{res}_{(\bm,k)}$ in Eq.(104) for a given class $k$
    !
    use find_all_roots_mod, only : customgrid,ncustom,niter,relerr_allroots, &
        xcustom,nroots,roots
    use get_matrix_mod,     only : iclass
    use form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end,sigma_class
    use resint_mod,         only : nmodes,marr,narr,twopim2,rm3,delint_mode,respoints_jp, &
        psiast_res,taub_res,delphi_res,resline_unit, &
        resline_unit_is_private,resline_diag_unit,ledger_harmonic_scans, &
        ledger_searched_zero,ledger_root_count,ledger_failure_count
    use orbit_dim_mod,      only : neqm
    use global_invariants,  only : dtau,toten,perpinv,cE_ref,Phi_eff
    use sample_matrix_mod,  only : npoi,xarr
    use logging_mod,        only : tee_message
    use interp_cache_mod,   only : interp_cache_reset
    use field_sub,          only : psif
    use field_eq_mod,       only : psi_axis,psi_sep
    use resonance_mode_bounds_mod, only : canonical_flux_outside_lcfs
    use resonance_status_mod, only : resonance_status_success, &
        resonance_status_root_failure, &
        resonance_status_starter_failure,resonance_status_from_root_error, &
        resonance_status_nonfinite_weight,classify_resonance_root, &
        resonance_is_finite
    use potato_symbolic_kernel_mod, only : potato_resonance_torque_kernel
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    !
    implicit none
    !
    double precision, parameter :: pi=3.14159265358979d0, twopi=2.d0*pi
    !
    integer          :: mode,iroot,ierr,ierr_pertham,status_weight
    integer, intent(out) :: ierr_out
    logical, intent(in) :: append
    integer :: old_nrespoi
    double precision :: xbeg,xend
    double precision :: rescond,dresconddx,dpsiastdx
    double precision :: one_res,sigma,delta_R,Rst,xi,dxi_dx,dpsiast_dRst,absHn2
    double precision :: root_factor
    double precision :: bmod_st,phi_elec_st
    double precision :: fmaxw,A1ast,A2ast,thermodynamic_force
    double precision :: delta_root_weight
    double precision :: dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi
    double precision, allocatable :: old_w_res(:),old_z_res(:,:),old_taub(:)
    double precision, allocatable :: combined_w_res(:),combined_z_res(:,:), &
                                     combined_taub(:)
    double precision, dimension(neqm) :: z
    logical :: wall_zero
    character(len=256) :: msg
    !
    ierr_out=resonance_status_success
    sigma=sigma_class(iclass)
    delta_R=R_class_end(iclass)-R_class_beg(iclass)
    ! sample_class_doublecount may trim either class endpoint for a finite
    ! harmonic guard, an orbit-return boundary, or wall closure.  Rebuilding
    ! the nominal class chart here would reintroduce the excluded singular or
    ! invalid endpoint into the resonance root search.  The completed class
    ! grid is the authoritative physical domain for interpolation and roots.
    if(npoi.le.1) then
        write(msg,'(A,I0,A,I0,A,I0)') &
            'resonance root search has no class interval class=',iclass, &
            ' iftype=',ifuntype(iclass),' npoi=',npoi
        !$omp critical (reslines_log)
        call tee_message(trim(msg))
        !$omp end critical (reslines_log)
        ierr_out=resonance_status_root_failure
        return
    endif
    xbeg=xarr(1)
    xend=xarr(npoi)
    customgrid=.true.
    ncustom=npoi
    allocate(xcustom(ncustom))
    xcustom=xarr
    !
    ! Mode loop.  The outer energy loop owns the concurrency; this per-class
    ! scan is deliberately serial so a failed harmonic can abort the class
    ! without an OpenMP worksharing loop hiding the status.
    ! Each iteration runs its own find_all_roots and per-root starter+pertham state.
    ! toten,perpinv and the root-search grid are copied into each worker as read-only
    ! class invariants. The threadprivate form-class bounds it reads
    ! (ifuntype,R_class_beg,R_class_end,sigma_class) are copyin'd.  The get_rescond by-products psiast_res,taub_res,delphi_res are
    ! threadprivate module state (resint_mod), NOT host-local privates: get_rescond
    ! is passed as the dummy root function to find_all_roots, and gfortran's
    ! trampoline for that call does not honor a host-local private.  reslines
    ! write(31415) and tee_message go in a critical section.  Each thread resets its
    ! interp cache at entry so it drops the previous class's grid before reusing it.
    call interp_cache_reset
    ledger_harmonic_scans=ledger_harmonic_scans+nmodes
    do mode=1,nmodes
        if(ierr_out.ne.resonance_status_success) exit
        twopim2=twopi*dble(marr(mode))
        rm3=dble(narr(mode))
        if(.not.append) delint_mode(mode)=0.d0
        old_nrespoi=0
        if(append) then
            old_nrespoi=respoints_jp(mode,iclass)%nrespoi
            if(old_nrespoi.gt.0) then
                allocate(old_w_res(old_nrespoi),old_z_res(5,old_nrespoi), &
                         old_taub(old_nrespoi))
                old_w_res=respoints_jp(mode,iclass)%w_res
                old_z_res=respoints_jp(mode,iclass)%z_res
                old_taub=respoints_jp(mode,iclass)%taub
            endif
            if(allocated(respoints_jp(mode,iclass)%w_res)) then
                deallocate(respoints_jp(mode,iclass)%w_res, &
                           respoints_jp(mode,iclass)%z_res, &
                           respoints_jp(mode,iclass)%taub)
            endif
        endif
        respoints_jp(mode,iclass)%nrespoi=0
        respoints_jp(mode,iclass)%scan_status=resonance_status_success
        !
        call find_all_roots_bracketed(get_rescond,xbeg,xend,ierr)
        !
        if(ierr.ne.0) then
            write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES14.6,A,ES14.6)') &
                'resonance root search failed ierr=',ierr, &
                ' class=',iclass,' iftype=',ifuntype(iclass), &
                ' mode=',marr(mode),'/',narr(mode), &
                ' xbeg=',xbeg,' xend=',xend
            !$omp critical (reslines_log)
            call tee_message( &
                'integrate_class_resonances: error in find_all_roots')
            call tee_message(trim(msg))
            !$omp end critical (reslines_log)
            ierr_out=resonance_status_from_root_error(ierr)
            ledger_failure_count=ledger_failure_count+1
            respoints_jp(mode,iclass)%scan_status=ierr_out
            respoints_jp(mode,iclass)%toten_res=toten
            respoints_jp(mode,iclass)%perpinv_res=perpinv
            if(allocated(old_w_res)) deallocate(old_w_res)
            if(allocated(old_z_res)) deallocate(old_z_res)
            if(allocated(old_taub)) deallocate(old_taub)
            exit
        endif
        !
        respoints_jp(mode,iclass)%nrespoi=nroots
        respoints_jp(mode,iclass)%toten_res=toten
        respoints_jp(mode,iclass)%perpinv_res=perpinv
        if(nroots.eq.0) then
            ledger_searched_zero=ledger_searched_zero+1
            if(old_nrespoi.gt.0) then
                call move_alloc(old_w_res,respoints_jp(mode,iclass)%w_res)
                call move_alloc(old_z_res,respoints_jp(mode,iclass)%z_res)
                call move_alloc(old_taub,respoints_jp(mode,iclass)%taub)
                respoints_jp(mode,iclass)%nrespoi=old_nrespoi
            endif
            cycle
        endif
        ledger_root_count=ledger_root_count+nroots
        allocate(respoints_jp(mode,iclass)%w_res(nroots), &
            respoints_jp(mode,iclass)%z_res(5,nroots), &
            respoints_jp(mode,iclass)%taub(nroots))
        !
        do iroot=1,nroots
            !
            call get_rescond(roots(iroot),rescond,dresconddx)
            call xi_func(ifuntype(iclass),roots(iroot),xi,dxi_dx)
            !
            Rst=R_class_beg(iclass)+delta_R*xi
            !
            call starter_doublecount(toten,perpinv,sigma,Rst,   &
                psiast_res,dpsiast_dRst,z,ierr)
            !
            if(ierr.ne.0) then
                !$omp critical (reslines_log)
                call tee_message( &
                    'integrate_class_resonances: error in starter_doublecount')
                !$omp end critical (reslines_log)
                ierr_out=resonance_status_starter_failure
                ledger_failure_count=ledger_failure_count+1
                respoints_jp(mode,iclass)%scan_status=ierr_out
                if(allocated(respoints_jp(mode,iclass)%w_res)) then
                    deallocate(respoints_jp(mode,iclass)%w_res, &
                               respoints_jp(mode,iclass)%z_res, &
                               respoints_jp(mode,iclass)%taub)
                endif
                exit
            endif
            !
            if(.true.) then
                if(resline_unit_is_private) then
                    write(resline_unit,*) toten,perpinv,psiast_res,marr(mode),narr(mode)
                    call get_bmod_and_Phi(z(1:3),bmod_st,phi_elec_st)
                    write(resline_diag_unit,*) toten,perpinv,psiast_res,marr(mode),narr(mode), &
                        Rst,z(3),psif,psiast_res-psif,roots(iroot),iclass,taub_res,delphi_res, &
                        z(4),z(5),bmod_st,phi_elec_st
                else
                    !$omp critical (reslines_log)
                    write(resline_unit,*) toten,perpinv,psiast_res,marr(mode),narr(mode)
                    call get_bmod_and_Phi(z(1:3),bmod_st,phi_elec_st)
                    write(resline_diag_unit,*) toten,perpinv,psiast_res,marr(mode),narr(mode), &
                        Rst,z(3),psif,psiast_res-psif,roots(iroot),iclass,taub_res,delphi_res, &
                        z(4),z(5),bmod_st,phi_elec_st
                    !$omp end critical (reslines_log)
                endif
            endif
            !
            respoints_jp(mode,iclass)%z_res(:,iroot)=z
            !
            dpsiastdx=dpsiast_dRst*delta_R*dxi_dx !$\difp{\psi^\ast}{x}$
            !
            wall_zero=canonical_flux_outside_lcfs(psiast_res,psi_axis,psi_sep)
            if(wall_zero) then
                ! The wall/SOL exclusion is a certified zero only after the
                ! resonance root factor has been checked below.  It must not
                ! mask a tangent root or an invalid derivative.
                absHn2=0.d0
            else
                call pertham(z,absHn2,ierr_pertham)
                if(ierr_pertham.ne.0) then
                    !$omp critical (reslines_log)
                    call tee_message( &
                        'integrate_class_resonances: find_bounce/pertham failed')
                    !$omp end critical (reslines_log)
                    ierr_out=ierr_pertham
                    ledger_failure_count=ledger_failure_count+1
                    respoints_jp(mode,iclass)%scan_status=ierr_out
                    deallocate(respoints_jp(mode,iclass)%w_res, &
                               respoints_jp(mode,iclass)%z_res, &
                               respoints_jp(mode,iclass)%taub)
                    exit
                endif
            endif
            call classify_resonance_root(dresconddx,dpsiastdx,absHn2, &
                wall_zero,root_factor,status_weight)
            if(.not.resonance_is_finite(rescond) .or. &
               .not.resonance_is_finite(taub_res)) then
                status_weight=resonance_status_nonfinite_weight
            endif
            if(status_weight.ne.0) then
                !$omp critical (reslines_log)
                call tee_message( &
                    'integrate_class_resonances: unresolved resonance root weight')
                !$omp end critical (reslines_log)
                print *,'integrate_class_resonances: resonance status = ',status_weight
                ierr_out=status_weight
                ledger_failure_count=ledger_failure_count+1
                respoints_jp(mode,iclass)%scan_status=ierr_out
                deallocate(respoints_jp(mode,iclass)%w_res, &
                           respoints_jp(mode,iclass)%z_res, &
                           respoints_jp(mode,iclass)%taub)
                exit
            endif
            if(wall_zero) then
                ! Exact wall-zero: the finite delta/root factor was certified,
                ! so the excluded residence/perturbation contributes exactly 0.
                one_res=0.d0
            else
                call equilmaxw(psiast_res,fmaxw)
                call denstemp_of_psi(psiast_res,dens,temp,ddens,dtemp)
                call phielec_of_psi(psiast_res,phi_elec,dPhi_dpsi)
                !
                ! toten,perpinv are never written in this loop (pertham writes toten_orb,
                ! perpinv_orb instead), so they still hold the class invariants set on entry --
                ! no save/restore needed.
                !
                ! Non-local thermodynamic forces Eq.(95) (former Eq.(87)):
                A2ast=dtemp/temp
                A1ast=ddens/dens+dPhi_dpsi/temp-1.5d0*A2ast
                !
                thermodynamic_force=A1ast+A2ast*(toten-phi_elec)/temp
                ! Eq. (17) delta reduction, the leading -pi^(3/2)/4, and the
                ! two bounce-time factors are generated together.  The mode
                ! number enters as abs(n), so a conjugate representation does
                ! not reverse the torque.
                call potato_resonance_torque_kernel(dpsiastdx,dresconddx,rm3,absHn2, &
                    fmaxw,Phi_eff,taub_res,thermodynamic_force,cE_ref, &
                    delta_root_weight,one_res)
                if(.not.ieee_is_finite(one_res)) then
                    !$omp critical (reslines_log)
                    call tee_message( &
                        'integrate_class_resonances: nonfinite resonance weight')
                    !$omp end critical (reslines_log)
                    ierr_out=resonance_status_nonfinite_weight
                    ledger_failure_count=ledger_failure_count+1
                    respoints_jp(mode,iclass)%scan_status=ierr_out
                    deallocate(respoints_jp(mode,iclass)%w_res, &
                               respoints_jp(mode,iclass)%z_res, &
                               respoints_jp(mode,iclass)%taub)
                    exit
                endif
            endif
            !
            respoints_jp(mode,iclass)%w_res(iroot)=one_res
            respoints_jp(mode,iclass)%taub(iroot)=taub_res
            delint_mode(mode)=delint_mode(mode)+one_res
        enddo
        if(ierr_out.eq.resonance_status_success) then
            if(old_nrespoi.gt.0) then
                allocate(combined_w_res(old_nrespoi+nroots), &
                         combined_z_res(5,old_nrespoi+nroots), &
                         combined_taub(old_nrespoi+nroots))
                combined_w_res(1:old_nrespoi)=old_w_res
                combined_z_res(:,1:old_nrespoi)=old_z_res
                combined_taub(1:old_nrespoi)=old_taub
                combined_w_res(old_nrespoi+1:)=respoints_jp(mode,iclass)%w_res
                combined_z_res(:,old_nrespoi+1:)=respoints_jp(mode,iclass)%z_res
                combined_taub(old_nrespoi+1:)=respoints_jp(mode,iclass)%taub
                deallocate(respoints_jp(mode,iclass)%w_res, &
                           respoints_jp(mode,iclass)%z_res, &
                           respoints_jp(mode,iclass)%taub)
                call move_alloc(combined_w_res,respoints_jp(mode,iclass)%w_res)
                call move_alloc(combined_z_res,respoints_jp(mode,iclass)%z_res)
                call move_alloc(combined_taub,respoints_jp(mode,iclass)%taub)
                respoints_jp(mode,iclass)%nrespoi=old_nrespoi+nroots
            endif
        endif
        if(allocated(old_w_res)) deallocate(old_w_res)
        if(allocated(old_z_res)) deallocate(old_z_res)
        if(allocated(old_taub)) deallocate(old_taub)
    enddo
    !
    customgrid=.false.
    deallocate(xcustom)
    !
    !------------
contains
    !------------
    !
    subroutine get_rescond(x,rescond,dresconddx)
        !
        ! Computes resonance condition $F(x)=\Delta\varphi_b+2\pi m_2/m_3$
        ! and its derivative $F^\prime(x)$ for $F(x)=0$ root finding.
        ! Computes as by-products normalized toroidal momentum $\psi^\ast$,
        ! bounce time $\tau_b$ and toroidal displacement $\Delta\varphi_b$.
        !
        use sample_matrix_mod, only : n1
        !
        implicit none
        !
        double precision :: x,rescond,dresconddx
        double precision, dimension(n1) :: vec,dvec
        !
        call interpolate_class_doublecount(x,vec,dvec)
        !
        psiast_res=vec(1) ! $\psi^\ast$
        taub_res=vec(2) ! $\tau_b$
        delphi_res=vec(3) ! $\Delta\varphi_b$
        rescond=delphi_res+twopim2/rm3
        dresconddx=dvec(3)
        !
    end subroutine get_rescond
    !
    !------------
    !
end subroutine integrate_class_resonances
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_matrix_res
    !
    use sample_matrix_out_mod,        only : n1,n2,x,amat,icount, &
        topology_signature,topology_error,topology_context_h, &
        topology_signature_of_classes,topology_probe_only
    use sample_matrix_mod, only : sample_npoi => npoi, sample_x => x, &
        sample_xbeg => xbeg, sample_xend => xend, matrix_eval_error, &
        sample_boundary_error => matrix_boundary_error
    use global_invariants,            only : toten,perpinv
    use form_classes_doublecount_mod, only : nclasses,ifuntype,sigma_class, &
        R_class_beg,R_class_end
    use get_matrix_mod,               only : iclass,delphi_max, &
        class_return_failure_reason
    use sample_class_status_mod,       only : sample_class_success, &
        sample_class_no_resonance
    use resonance_status_mod,          only : resonance_status_success
    use resint_mod,                   only : nmodes,nperp_max,respoints_jp, &
        respoints_all,respoints_all_tmp,ledger_class_evaluations
    use cc_mod,                       only : wrbounds,dowrite
    use orbit_dim_mod,                only : write_orb,orbit_failure_stage,orbit_wall_loss
    use logging_mod,                  only : tee_message
    use bounds_fixpoints_mod,         only : region_set_t
    !
    logical :: classes_talk
    !
    integer :: ierr,mode,ierr_resonance,ierr_segment
    logical :: has_next
    character(len=256) :: msg
    type(region_set_t) :: regions
    external :: sample_class_next_segment
    !
    wrbounds=.false.
    dowrite=.false.
    write_orb=.false.
    classes_talk=.false.
    topology_signature=0
    topology_error=0
    topology_context_h=toten
    ! Topology probes run before sample_matrix_out allocates its matrix
    ! workspace; they need only the class signature.  The production callback
    ! owns an allocated amat and initializes it below.
    if(.not.topology_probe_only) amat=0.d0
    !
    perpinv=x
    !
    call find_bounds_fixpoints(regions,ierr)
    !
    if(ierr.ne.0) then
        write(msg, '(A,I0)') &
            'get_matrix_res: find_bounds_fixpoints ierr = ', ierr
        call tee_message(trim(msg))
        topology_error=ierr
        return
    endif
    !
    call form_classes_doublecount(regions,classes_talk,ierr)
    !
    if(ierr.ne.0) then
        write(msg, '(A,I0)') &
            'get_matrix_res: form_classes ierr = ', ierr
        call tee_message(trim(msg))
        topology_error=ierr
        return
    endif
    if(nclasses.le.0) then
        ! The requested J_perp slice can have no orbit intersecting the
        ! clipped rho_pol domain.  That is a certified zero contribution,
        ! not a callback failure or an unresolved topology transition.
        topology_signature=topology_signature_of_classes(0,ifuntype,sigma_class)
        return
    endif
    topology_signature=topology_signature_of_classes(nclasses,ifuntype,sigma_class)
    if(topology_probe_only) return
    ledger_class_evaluations=ledger_class_evaluations+nclasses
    !
    allocate(respoints_jp(nmodes,nclasses))
    amat=0.d0
    do iclass=1,nclasses
        do mode=1,nmodes
            respoints_jp(mode,iclass)%nrespoi=0
            respoints_jp(mode,iclass)%scan_status=resonance_status_success
            respoints_jp(mode,iclass)%toten_res=toten
            respoints_jp(mode,iclass)%perpinv_res=perpinv
        enddo
    enddo
    !
    do iclass=1,nclasses
        !
        call sample_class_doublecount(1,ierr)
        !
        if(ierr.eq.sample_class_success) then
            !
            call integrate_class_resonances(ierr_resonance,.false.)
            if(ierr_resonance.ne.0) then
                write(msg, '(A,I0)') &
                    'get_matrix_res: unresolved resonance status ',ierr_resonance
                call tee_message(trim(msg))
                topology_error=ierr_resonance
                deallocate(respoints_jp)
                return
            endif
            do
                call sample_class_next_segment(ierr_segment,has_next)
                if(ierr_segment.ne.sample_class_success) then
                    write(msg, '(A,I0,A,I0)') &
                        'get_matrix_res: next class segment failed status=', &
                        ierr_segment,' class=',iclass
                    call tee_message(trim(msg))
                    topology_error=ierr_segment
                    deallocate(respoints_jp)
                    return
                endif
                if(.not.has_next) exit
                call integrate_class_resonances(ierr_resonance,.true.)
                if(ierr_resonance.ne.0) then
                    write(msg, '(A,I0)') &
                        'get_matrix_res: unresolved resonance status ',ierr_resonance
                    call tee_message(trim(msg))
                    topology_error=ierr_resonance
                    deallocate(respoints_jp)
                    return
                endif
            enddo
            !
            do mode=1,nmodes
                ! amat(:,1) - sum over classes and resonances within each class:
                if(respoints_jp(mode,iclass)%nrespoi.gt.0) then
                    amat(mode,1)=amat(mode,1)+sum(respoints_jp(mode,iclass)%w_res)
                endif
            enddo
        else if(ierr.eq.sample_class_no_resonance) then
            ! The harmonic guard proved that this class has no resonance;
            ! its initialized zero records remain part of the class topology.
            cycle
        else
            write(msg, '(A,I0,A,I0,A,I0,A,F6.1,A,2ES14.6,A,2ES14.6,A,ES14.6,A,I0,A,I0,A,I0,A,L1)') &
                'get_matrix_res: sample_class error=',ierr, &
                ' class=',iclass,' iftype=',ifuntype(iclass), &
                ' sigma=',sigma_class(iclass), &
                ' Rbeg,Rend=',R_class_beg(iclass),R_class_end(iclass), &
                ' xbeg,xend=',sample_xbeg,sample_xend, &
                ' sample_x=',sample_x,' npoi=',sample_npoi, &
                ' matrix_error=',matrix_eval_error, &
                ' orbit_stage=',orbit_failure_stage, &
                ' wall_loss=',orbit_wall_loss
            call tee_message(trim(msg))
            write(msg, '(A,I0)') &
                'get_matrix_res: sample_class boundary_error=', sample_boundary_error
            call tee_message(trim(msg))
            write(msg, '(A,ES14.6)') &
                'get_matrix_res: worker delphi_max=', delphi_max
            call tee_message(trim(msg))
            write(msg, '(A,I0)') &
                'get_matrix_res: return_failure_reason=', class_return_failure_reason
            call tee_message(trim(msg))
            topology_error=ierr
            deallocate(respoints_jp)
            return
        endif
    enddo

    !
    icount=icount+1
    if(icount.gt.nperp_max) then
        if(nperp_max.le.0) then
            nperp_max=1
            if(allocated(respoints_all)) deallocate(respoints_all)
            allocate(respoints_all(nperp_max))
        else
            allocate(respoints_all_tmp(nperp_max))
            respoints_all_tmp=respoints_all
            deallocate(respoints_all)
            allocate(respoints_all(2*nperp_max))
            respoints_all(1:nperp_max)=respoints_all_tmp
            deallocate(respoints_all_tmp)
            nperp_max=2*nperp_max
        endif
    endif
    respoints_all(icount)%nclasses=nclasses
    allocate(respoints_all(icount)%respoints_jp(nmodes,nclasses))
    respoints_all(icount)%respoints_jp=respoints_jp
    deallocate(respoints_jp)
    !print *,icount  !uncomment to see some life
    !
end subroutine get_matrix_res
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_matrix_res_topology_boundaries(nmax,ncandidates,candidates,ierr)
    use global_invariants, only : toten
    use sample_matrix_out_mod, only : topology_context_h
    integer, intent(in) :: nmax
    integer, intent(out) :: ncandidates,ierr
    double precision, intent(out) :: candidates(nmax)
    external :: find_jperp_topology_boundaries

    topology_context_h=toten
    call find_jperp_topology_boundaries(candidates,nmax,ncandidates,ierr)
end subroutine get_matrix_res_topology_boundaries
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine get_matrix_res_contribution_envelope(gap_lo,gap_hi,envelope,ierr)
    use sample_matrix_out_mod, only : sample_matrix_contribution_unresolved
    use topology_gap_quadrature_mod, only : estimate_topology_gap_from_samples
    implicit none

    double precision, intent(in) :: gap_lo,gap_hi
    double precision, intent(out) :: envelope
    integer, intent(out) :: ierr
! The sampler consumes an average envelope, which is multiplied by the gap
! width in its existing callback ABI.  The provider uses the already sampled
! one-sided branch data and the generated inverse-square-root limiting form;
! it does not reintegrate the full resonance solver at new J_perp points.
    call estimate_topology_gap_from_samples(gap_lo,gap_hi,envelope,ierr)
    if(ierr.ne.0 .and. ierr.ne.sample_matrix_contribution_unresolved) then
        ierr=sample_matrix_contribution_unresolved
    endif
end subroutine get_matrix_res_contribution_envelope
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine resonant_torque
    !
    !
    use field_sub, only : psif
    use field_eq_mod, only : nrad,nzet,rad,zet,psi_sep
    use poicut_mod,        only : rmagaxis,zmagaxis,psimagaxis,psi_bou,rhopol_bou
    use global_invariants, only : dtau,toten,perpinv,sigma
    use poicut_mod,        only : Rbou_lfs,Zbou_lfs
    use get_matrix_mod,    only : iclass,delphi_max,relerror,relmargin
    use form_classes_doublecount_mod, only : nclasses
    use orbit_dim_mod,     only : numbasef,write_orb,next,orbit_wall_loss, &
        orbit_failure_stage, &
        clip_resonance_classes
    use cc_mod,             only : dowrite,wrbounds
    use find_all_roots_mod, only : customgrid,ncustom,niter,nsearch_min, &
        relerr_allroots,root_eval_valid,root_eval_error
    use potato_boundary_scan_mod, only : fixedpoint_scan_sigma, &
        fixedpoint_scan_branch,fixedpoint_scan_left,fixedpoint_scan_right
    use resonance_mode_bounds_mod, only : set_resonance_harmonic_guard, &
        write_resonance_harmonic_set, harmonic_guard_success
    use potato_boundary_scan_mod, only : jperpmax_success,jperpmax_status, &
        jperpmax_certified,jperpmax_witness,jperpmax_upper_bound
    use resint_mod,        only : nmodes,marr,narr,delint_mode,respoints_jp,respoints_all,nperp_max, &
        respoints_all_tmp,respoint,resline_unit,resline_diag_unit, &
        resline_unit_is_private,resline_diag_unit_is_private, &
        ledger_class_evaluations,ledger_harmonic_scans,ledger_searched_zero, &
        ledger_root_count,ledger_failure_count,ledger_jperp_samples
    use sample_matrix_out_mod, only : nlagr,n1,n2,npoi,itermax,x,amat,icount,xbeg,xend,eps, &
        ind_hist,xarr,amat_arr,topology_arr,topology_error,topology_signature, &
        topology_partition_tol, &
        topology_context_h,topology_probe_only,sample_matrix_preserve_history, &
        topology_candidate_count,topology_transition_count,topology_gap_bound, &
        topology_gap_measure,topology_gap_geometric_bound, &
        topology_contribution_error_bound,topology_contribution_error_certified
    use potato_input_mod,  only : nbox, nenerg_input => nenerg, &
        thermen_max_input => thermen_max, &
        adaptive_jperp, npoi_init, nlagr_sampling, &
        eps_sampling, itermax_sampling, &
        topology_partition_tol_input => topology_partition_tol, &
        topology_partition_tol_refined, topology_refinement_lane, &
        require_topology_contribution_bound
    use logging_mod,       only : tee_message
    !$ use omp_lib, only : omp_set_max_active_levels, omp_set_num_threads

    implicit none
    !
    double precision, parameter :: pi=3.14159265358979d0
    integer :: nr,nz,ir,iz,i,k,iperp,nperp,ierr,nprof,ienerg,nenerg
    integer :: omp_threads_env_len
    integer :: nrespoints,unit1901,unit1902,ienerg_begin
    double precision :: rbeg,hr,zbeg,hz,weight,psi,psipow
    double precision :: bmod,phi_elec,phi_elec_min,phi_elec_max
    double precision :: toten_min,toten_max,thermen_max,toten_range
    double precision :: omdens,trapez_fac,perpinv_max
    double precision :: torque_int,torque_int_loc
    double precision :: xjperp,xenerg,totxint,step_energ
    double precision :: time_beg,time_end
    double precision :: dens, temp, ddens, dtemp
    character(len=256) :: msg
    double precision, dimension(:), allocatable :: torque_int_modes
    double precision, dimension(:), allocatable :: torque_int_modes_loc
    double precision, dimension(:), allocatable :: sbox
    double precision, dimension(:), allocatable :: taubox
    double precision, dimension(:), allocatable :: torquebox
    double precision, dimension(:), allocatable :: torquebox_loc
    double precision, dimension(:), allocatable :: torque_int_energy
    double precision, dimension(:,:), allocatable :: torque_int_modes_energy
    double precision, dimension(:,:), allocatable :: torquebox_energy
    integer, dimension(:), allocatable :: ledger_classes_energy,ledger_harmonics_energy, &
        ledger_searched_zero_energy,ledger_roots_energy,ledger_failures_energy, &
        ledger_jperp_samples_energy,ledger_jperpmax_status_energy
    double precision, dimension(:), allocatable :: ledger_gap_energy,ledger_gap_bound_energy, &
        ledger_contribution_bound_energy,ledger_jperp_witness_energy, &
        ledger_jperp_upper_energy
    logical, dimension(:), allocatable :: ledger_contribution_certified_energy
    ! Per-resonance-point box times, filled by the parallel time_in_box pass and
    ! consumed by the serial accumulate/write pass below.
    double precision, dimension(:,:), allocatable :: taubox_all
    !
    external :: get_matrix_res,get_matrix_res_topology_boundaries, &
        get_matrix_res_contribution_envelope, sample_matrix_out_partitioned_certified
    !
    allocate(sbox(nbox), taubox(nbox), torquebox(nbox))
    !
    numbasef=0 !no extra integrals sampled, pure orbit integration
    call linspace(1d0/nbox, 1d0, nbox, sbox)
    !
    ! Bound only the optional torque root search with the exact finite supplied
    ! harmonic set.  The signed m/n remains in g=Delta_phi+2*pi*m/n; the
    ! generated guard uses |m/n| only for the symmetric search extent.
    call set_resonance_harmonic_guard(marr,narr,delphi_max,ierr)
    if(ierr.ne.harmonic_guard_success) then
        call tee_message('resonant_torque: invalid harmonic guard/set')
        error stop 2
    endif
    write(msg, '(A,ES14.6)') &
        'torque root search bounded to |delphi_b| <= ', delphi_max
    call tee_message(trim(msg))
    !
    ! Find minimum and maximum values of electrostatic potential in the computation domain:
    !
    call find_Phiminmax(phi_elec_min,phi_elec_max)
    !
    thermen_max=thermen_max_input
    !
    call denstemp_of_psi(psimagaxis, dens, temp, ddens, dtemp)
    !
    thermen_max=thermen_max*temp !maximum kinetic energy in units of reference energy
    !
    ! Energy integration limits:
    toten_min=phi_elec_min
    toten_max=thermen_max+phi_elec_max
    toten_range=toten_max-toten_min
    !
    write(msg, '(A,ES14.6)') &
        'maximum kinetic energy = ', thermen_max
    call tee_message(trim(msg))
    write(msg, '(A,ES14.6)') &
        'minimum potential energy = ', phi_elec_min
    call tee_message(trim(msg))
    write(msg, '(A,ES14.6)') &
        'maximum potential energy = ', phi_elec_max
    call tee_message(trim(msg))
    write(msg, '(A,ES14.6)') &
        'minimum total energy = ', toten_min
    call tee_message(trim(msg))
    write(msg, '(A,ES14.6)') &
        'maximum total energy = ', toten_max
    call tee_message(trim(msg))
    !
    nenerg=nenerg_input
    !
    torque_int=0.d0
    allocate(torque_int_modes(nmodes))
    torque_int_modes=0.d0
    allocate(torque_int_energy(nenerg),torque_int_modes_energy(nmodes,nenerg), &
        torquebox_energy(nbox,nenerg))
    allocate(ledger_classes_energy(nenerg),ledger_harmonics_energy(nenerg), &
        ledger_searched_zero_energy(nenerg),ledger_roots_energy(nenerg), &
        ledger_failures_energy(nenerg),ledger_jperp_samples_energy(nenerg), &
        ledger_jperpmax_status_energy(nenerg),ledger_gap_energy(nenerg), &
        ledger_gap_bound_energy(nenerg),ledger_contribution_bound_energy(nenerg), &
        ledger_jperp_witness_energy(nenerg),ledger_jperp_upper_energy(nenerg), &
        ledger_contribution_certified_energy(nenerg))
    torque_int_energy=0.d0
    torque_int_modes_energy=0.d0
    torquebox_energy=0.d0
    ledger_classes_energy=0
    ledger_harmonics_energy=0
    ledger_searched_zero_energy=0
    ledger_roots_energy=0
    ledger_failures_energy=0
    ledger_jperp_samples_energy=0
    ledger_jperpmax_status_energy=0
    ledger_gap_energy=0.d0
    ledger_gap_bound_energy=0.d0
    ledger_contribution_bound_energy=huge(1.d0)
    ledger_jperp_witness_energy=0.d0
    ledger_jperp_upper_energy=huge(1.d0)
    ledger_contribution_certified_energy=.false.
    !
    torquebox=0.d0
    !
    step_energ=toten_range/dble(nenerg) !integration step over total energy
    ienerg_begin=2
    !step_energ=0.22521463755624047d0
    !
    omp_threads_env_len=0
    !$ call get_environment_variable('OMP_NUM_THREADS',length=omp_threads_env_len)
    !$ if(omp_threads_env_len.eq.0) call omp_set_num_threads(16)
    !$ call omp_set_max_active_levels(1)
    !$omp parallel do default(shared) schedule(dynamic) &
    !$omp   private(xenerg,perpinv_max,trapez_fac,nrespoints,i,k,iperp,ierr,nperp) &
    !$omp   private(time_beg,time_end,xjperp,torque_int_loc,msg,taubox_all) &
    !$omp   private(unit1901,unit1902,torque_int_modes_loc,torquebox_loc) &
    !$omp   copyin(dtau,toten,perpinv,sigma,relerror,relmargin,delphi_max,iclass, &
        !$omp         write_orb,next,orbit_wall_loss,orbit_failure_stage, &
        !$omp         dowrite,wrbounds, &
    !$omp         customgrid,ncustom,niter,nsearch_min,relerr_allroots, &
    !$omp         root_eval_valid,root_eval_error,fixedpoint_scan_sigma, &
    !$omp         fixedpoint_scan_branch,fixedpoint_scan_left,fixedpoint_scan_right, &
    !$omp         topology_partition_tol,topology_error,topology_signature, &
    !$omp         topology_context_h,topology_probe_only, &
    !$omp         sample_matrix_preserve_history,topology_candidate_count, &
    !$omp         topology_transition_count,topology_gap_measure, &
    !$omp         topology_gap_geometric_bound,topology_gap_bound, &
    !$omp         topology_contribution_error_bound, &
    !$omp         topology_contribution_error_certified)
    do ienerg=ienerg_begin,nenerg
        xenerg=(dble(ienerg)-0.5d0)/dble(nenerg)
        toten=toten_min+toten_range*xenerg
        ! size of result matrix in get_matrix_res:
        n1=nmodes
        n2=1
        ledger_class_evaluations=0
        ledger_harmonic_scans=0
        ledger_searched_zero=0
        ledger_root_count=0
        ledger_failure_count=0
        ledger_jperp_samples=0
        ! Set adaptive sampling parameters from input:
        nlagr=nlagr_sampling
        eps=eps_sampling
        itermax=itermax_sampling
        topology_partition_tol=topology_partition_tol_input
        if(allocated(delint_mode)) deallocate(delint_mode)
        allocate(delint_mode(nmodes))
        allocate(torque_int_modes_loc(nmodes),torquebox_loc(nbox))
        torque_int_modes_loc=0.d0
        torquebox_loc=0.d0
        resline_unit_is_private=.true.
        resline_diag_unit_is_private=.true.
        call open_energy_outputs(ienerg,adaptive_jperp,unit1901,unit1902)
        !toten =   -3.1211921097605737d0
        write(msg, '(A,ES22.14)') 'toten = ', toten
        call tee_message(trim(msg))
        !
        call cpu_time(time_beg)
        !
        ! find maximum possible value of J_perp for given total energy:
        !
        call find_Jperpmax(perpinv_max)
        if(jperpmax_status.ne.jperpmax_success .or. &
           .not.jperpmax_certified) then
            write(msg, '(A,I0,A,ES22.14,A,ES22.14)') &
                'resonant_torque: J_perp domain enclosure uncertified status = ', &
                jperpmax_status, ' witness = ', jperpmax_witness, &
                ' upper_bound = ', jperpmax_upper_bound
            call tee_message(trim(msg))
            error stop 2
        endif
        !
        if(adaptive_jperp) then
            !
            ! New, adaptive integration:
            !
            nperp_max=0 !initial size of respoints_all(nperp_max) - will be increased by sampling
            if(allocated(respoints_all)) deallocate(respoints_all)
            xbeg=0.d0 !lower integration limit over normalized J_perp
            xend=perpinv_max*0.9999d0 !upper integration limit over normalized J_perp
            npoi=npoi_init !initial equidistant grid size over J_perp
            icount=0 !initialize the counter of resonant points
            !
            ! Generate J_perp grid with function values:
            !
            call sample_matrix_out_partitioned_certified( &
                get_matrix_res,get_matrix_res_topology_boundaries,ierr, &
                get_matrix_res_contribution_envelope)
            if(ierr.ne.0) then
                write(msg, '(A,I0,A,ES22.14)') &
                    'resonant_torque: adaptive J_perp sampling failed ierr = ', &
                    ierr, ' at H = ', toten
                call tee_message(trim(msg))
                ! Do not reorder or integrate a truncated/non-certified grid.
                ! A bare STOP exits successfully with gfortran and used to make
                ! this failure look like a completed physics run.
                error stop 2
            endif
            if(require_topology_contribution_bound .and. &
               .not.topology_contribution_error_certified) then
                write(msg, '(A,ES22.14,A,ES22.14,A,ES22.14)') &
                    'resonant_torque: geometric topology gap has no certified contribution bound at H = ', &
                    toten, ' gap = ', topology_gap_measure, &
                    ' geometric_bound = ', topology_gap_geometric_bound
                call tee_message(trim(msg))
                error stop 2
            endif
            ledger_jperp_samples= npoi
            ledger_gap_energy(ienerg)=topology_gap_measure
            ledger_gap_bound_energy(ienerg)=topology_gap_geometric_bound
            ledger_contribution_bound_energy(ienerg)=topology_contribution_error_bound
            ledger_contribution_certified_energy(ienerg)= &
                topology_contribution_error_certified
            ledger_jperpmax_status_energy(ienerg)=jperpmax_status
            ledger_jperp_witness_energy(ienerg)=jperpmax_witness
            ledger_jperp_upper_energy(ienerg)=jperpmax_upper_bound
            !
            allocate(respoints_all_tmp(npoi))
            !
            ! reorder J_perp points to the increasing sequence:
            do i=1,npoi
                respoints_all_tmp(i)=respoints_all(ind_hist(i))
            enddo
            !
            deallocate(respoints_all)
            allocate(respoints_all(npoi))
            respoints_all=respoints_all_tmp
            deallocate(respoints_all_tmp)
            !
            ! count the number of resonant points and update their weights in accordance
            ! with trapezoidal integration rule over J_perp:
            nrespoints=0
            do iperp=1,npoi
                trapez_fac=0.d0
                if(iperp.gt.1) then
                    if(topology_arr(iperp).eq.topology_arr(iperp-1)) then
                        trapez_fac=trapez_fac+0.5d0*(xarr(iperp)-xarr(iperp-1))
                    endif
                endif
                if(iperp.lt.npoi) then
                    if(topology_arr(iperp).eq.topology_arr(iperp+1)) then
                        trapez_fac=trapez_fac+0.5d0*(xarr(iperp+1)-xarr(iperp))
                    endif
                endif
                do iclass=1,respoints_all(iperp)%nclasses
                    do i=1,nmodes
                        if(respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi.gt.0) then
                            respoints_all(iperp)%respoints_jp(i,iclass)%w_res = &
                                respoints_all(iperp)%respoints_jp(i,iclass)%w_res*trapez_fac
                            nrespoints=nrespoints &
                                +respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi
                        endif
                    enddo
                enddo
            enddo
            !
            allocate(respoint(nrespoints))
            !
            ! extract separate resonant points and update their weight in accorance with
            ! integration over total energy:
            nrespoints=0
            do iperp=1,npoi
                do iclass=1,respoints_all(iperp)%nclasses
                    do i=1,nmodes
                        do k=1,respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi
                            nrespoints=nrespoints+1
                            respoint(nrespoints)%toten_res                             &
                                =respoints_all(iperp)%respoints_jp(i,iclass)%toten_res
                            respoint(nrespoints)%perpinv_res                           &
                                =respoints_all(iperp)%respoints_jp(i,iclass)%perpinv_res
                            respoint(nrespoints)%w_res                                 &
                                =respoints_all(iperp)%respoints_jp(i,iclass)%w_res(k)    &
                                *step_energ
                            respoint(nrespoints)%z_res &
                                =respoints_all(iperp)%respoints_jp(i,iclass)%z_res(:,k)
                            respoint(nrespoints)%taub &
                                =respoints_all(iperp)%respoints_jp(i,iclass)%taub(k)
                        enddo
                    enddo
                enddo
            enddo
            !
            ! Computation of the integral torque:
            torque_int_loc=0.d0
            !
            ! Box counter (time_in_box) integrates each resonant orbit through the radial
            ! boxes via dvode -- the dominant per-point cost.  Points are independent and
            ! dvode_f90_m is threadprivate/thread-safe, so run them in parallel into the
            ! per-point taubox_all.  z_res(1:5) is the orbit start at the Poincare cut.
            allocate(taubox_all(nbox,nrespoints))
            !$omp parallel do default(shared) private(i) schedule(dynamic)
            do i=1,nrespoints
                if(respoint(i)%w_res.eq.0.d0) then
                    ! Zero-weight (SOL / unintegrable) point: no torque to deposit,
                    ! and its orbit is the expensive one to trace.  Skip it.
                    taubox_all(:,i)=0.d0
                else
                    call time_in_box(respoint(i)%z_res(1:5), nbox, sbox, &
                        respoint(i)%taub, taubox_all(:,i))
                endif
            enddo
            !$omp end parallel do
            !
            ! Serial accumulate/write pass: torque_int_loc is a running prefix sum written
            ! per row to fort.1901, so it stays serial and in order; torquebox sums the
            ! per-point box times.  Width of the energy interval is already in the weight.
            do i=1,nrespoints
                torque_int_loc=torque_int_loc+respoint(i)%w_res
                write(unit1901,*) respoint(i)%toten_res, &
                    respoint(i)%perpinv_res, &
                    torque_int_loc
                if(respoint(i)%w_res.ne.0.d0 .and. respoint(i)%taub.gt.0.d0) then
                    torquebox_loc=torquebox_loc &
                        +respoint(i)%w_res*taubox_all(:,i)/respoint(i)%taub
                endif
            enddo
            deallocate(taubox_all)
            write(unit1901,*) ' '
            !
            deallocate(respoint)
            write(msg, '(A,I0)') &
                'number of J_perp points = ', npoi
            call tee_message(trim(msg))
            !
            ! Sum up contributions of energy levels:
            torque_int_energy(ienerg)=torque_int_loc
            torquebox_energy(:,ienerg)=torquebox_loc
            write(msg, '(A,ES22.14)') &
                'method 1, torque_int_loc = ', torque_int_loc
            call tee_message(trim(msg))
            !
            ! Alternative computation via sampling matrix. Here distribution over resonances is also computed:
            do i=1,npoi
                !        write(10000,*) xarr(i),amat_arr(:,1,i)
                if(i.eq.1) then
                    torque_int_loc=0.d0
                elseif(topology_arr(i).eq.topology_arr(i-1)) then
                    torque_int_loc=torque_int_loc &
                        +0.5d0*(sum(amat_arr(:,1,i))+sum(amat_arr(:,1,i-1))) &
                        *(xarr(i)-xarr(i-1))*step_energ
                    torque_int_modes_loc=torque_int_modes_loc &
                        +0.5d0*(amat_arr(:,1,i)+amat_arr(:,1,i-1)) &
                        *(xarr(i)-xarr(i-1))*step_energ
                endif
                write(unit1902,*) xarr(i),torque_int_loc,torque_int_modes_loc
            enddo
            torque_int_modes_energy(:,ienerg)=torque_int_modes_loc
            write(unit1902,*) ' '
            write(msg, '(A,ES22.14)') &
                'method 2, torque_int_loc = ', torque_int_loc
            call tee_message(trim(msg))
            !
        else
            !
            ! Old, non-adaptive integration:
            ! Not prepared for box counting, use for testing only.
            torque_int_loc=0.d0
            nperp=2500 !5000 !100    !size of the integration grid over normalized J_perp
            if(.not.allocated(amat)) allocate(amat(n1,n2))
            icount=0
            nperp_max=0
            do iperp=1,nperp
                if(iperp.eq.nperp) then
                    trapez_fac=0.5d0
                else
                    trapez_fac=1.d0
                endif
                xjperp=dble(iperp)/dble(nperp)
                trapez_fac=trapez_fac*2.d0*xjperp/dble(nperp)
                x=perpinv_max*(1.d0-xjperp**2)
                !
                call get_matrix_res
                if(topology_error.ne.0) then
                    write(msg, '(A,I0,A,ES22.14,A,ES22.14)') &
                        'resonant_torque: fixed J_perp sampling failed ierr = ', &
                        topology_error, ' at H = ', toten, ' J = ', x
                    call tee_message(trim(msg))
                    error stop 2
                endif
                !
                torque_int_loc=torque_int_loc &
                    +perpinv_max*trapez_fac*sum(amat(:,1))*step_energ
                torque_int_modes_loc=torque_int_modes_loc &
                    +perpinv_max*trapez_fac*amat(:,1)*step_energ
                ! Subintegrand of dimensional integral over energy for a total torque as function of J_perp:
                write(unit1902,*) perpinv,torque_int_loc,torque_int_modes_loc
                write(msg, '(A,I0,A,I0,A,I0,A,I0)') &
                    'perpinv:', iperp, '/', nperp, &
                    ' toten:', ienerg, '/', nenerg
                call tee_message(trim(msg))
            enddo
            torque_int_energy(ienerg)=torque_int_loc
            torque_int_modes_energy(:,ienerg)=torque_int_modes_loc
        endif
        !
        ledger_classes_energy(ienerg)=ledger_class_evaluations
        ledger_harmonics_energy(ienerg)=ledger_harmonic_scans
        ledger_searched_zero_energy(ienerg)=ledger_searched_zero
        ledger_roots_energy(ienerg)=ledger_root_count
        ledger_failures_energy(ienerg)=ledger_failure_count
        ledger_jperp_samples_energy(ienerg)=ledger_jperp_samples
        call write_energy_ledger(ienerg,toten,torque_int_energy(ienerg), &
            torque_int_loc)
        call cpu_time(time_end)
        write(msg, '(A,I0,A,I0,A,F10.2,A)') &
            ' toten:', ienerg, '/', nenerg, &
            ' cpu time = ', time_end-time_beg, ' sec'
        call tee_message(trim(msg))
        !
        close(resline_unit)
        close(resline_diag_unit)
        if(adaptive_jperp) close(unit1901)
        close(unit1902)
        deallocate(torque_int_modes_loc,torquebox_loc)
        !
    enddo
    !$omp end parallel do
    !
    torque_int=0.d0
    torque_int_modes=0.d0
    torquebox=0.d0
    do ienerg=1,nenerg
        torque_int=torque_int+torque_int_energy(ienerg)
        torque_int_modes=torque_int_modes+torque_int_modes_energy(:,ienerg)
        torquebox=torquebox+torquebox_energy(:,ienerg)
    enddo
    !
    call merge_energy_files('fort.31415.energy.', 'fort.31415', nenerg)
    call merge_energy_files('fort.31415.diagnostics.energy.', &
        'fort.31415.diagnostics', nenerg)
    call merge_energy_files('potato_topology_ledger.energy.', &
        'potato_topology_ledger.dat', nenerg)
    if(adaptive_jperp) then
        call merge_energy_files('subint_ofH0int_104_vsJperp_fromresp.dat.energy.', &
            'subint_ofH0int_104_vsJperp_fromresp.dat', nenerg)
        call merge_jperp_files('subint_ofH0int_104_vsJperp_adapt.dat.energy.', &
            'subint_ofH0int_104_vsJperp_adapt.dat', nenerg, .false.)
    else
        call merge_jperp_files('subint_ofH0int_104_vsJperp_equi.dat.energy.', &
            'subint_ofH0int_104_vsJperp_equi.dat', nenerg, .true.)
    endif
    !
    ! Integral torque:
    write(msg, '(A,ES22.14)') &
        'resonant torque  = ', torque_int
    call tee_message(trim(msg))
    open(1,file='integral_torque.dat')
    write(1,*) torque_int
    close(1)
    call write_run_manifest
    !
    !Write box-counted integral torque:
    if(adaptive_jperp) then
        open(1,file='boxcounted_torque.dat')
        write(1,*) torquebox(1), 0d0, sbox(1)
        do i=2,nbox
            write(1,*) torquebox(i), sbox(i-1), sbox(i)
        enddo
        close(1)
    endif
    !
contains
    !
    subroutine energy_path(prefix, idx, path)
        !
        implicit none
        !
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: idx
        character(len=*), intent(out) :: path
        !
        write(path,'(A,I6.6)') trim(prefix), idx
        !
    end subroutine energy_path
    !
    subroutine write_energy_ledger(idx,energy,method1,method2)
        implicit none
        integer, intent(in) :: idx
        double precision, intent(in) :: energy,method1,method2
        integer :: ledger_unit,io
        character(len=256) :: path

        call energy_path('potato_topology_ledger.energy.',idx,path)
        open(newunit=ledger_unit,file=trim(path),status='replace',action='write', &
            iostat=io)
        if(io.ne.0) return
        write(ledger_unit,'(A)') &
            '# ienerg nenerg H class_evaluations harmonic_scans searched_zero roots failures '// &
            'jperp_samples jperpmax_status jperp_witness jperp_upper_bound '// &
            'geometric_gap geometric_bound contribution_error_bound contribution_certified '// &
            'torque_method1 torque_method2 active_tol refined_tol refinement_lane'
        write(ledger_unit,*) idx,nenerg,energy,ledger_classes_energy(idx), &
            ledger_harmonics_energy(idx),ledger_searched_zero_energy(idx), &
            ledger_roots_energy(idx),ledger_failures_energy(idx), &
            ledger_jperp_samples_energy(idx),ledger_jperpmax_status_energy(idx), &
            ledger_jperp_witness_energy(idx),ledger_jperp_upper_energy(idx), &
            ledger_gap_energy(idx),ledger_gap_bound_energy(idx), &
            ledger_contribution_bound_energy(idx), &
            ledger_contribution_certified_energy(idx),method1,method2, &
            topology_partition_tol_input,topology_partition_tol_refined, &
            topology_refinement_lane
        close(ledger_unit)
    end subroutine write_energy_ledger
    !
    subroutine write_run_manifest
        implicit none
        integer :: manifest_unit,io,harmonic_ierr

        open(newunit=manifest_unit,file='potato_run_manifest.dat', &
            status='replace',action='write',iostat=io)
        if(io.ne.0) return
        write(manifest_unit,'(A)') '# POTATO topology/conservation manifest v1'
        write(manifest_unit,'(A,I0)') 'nenerg = ',nenerg
        write(manifest_unit,'(A,ES24.16)') 'topology_partition_tol = ', &
            topology_partition_tol_input
        write(manifest_unit,'(A,ES24.16)') 'topology_partition_tol_refined = ', &
            topology_partition_tol_refined
        write(manifest_unit,'(A,L1)') 'topology_refinement_lane = ', &
            topology_refinement_lane
        write(manifest_unit,'(A,L1)') 'require_topology_contribution_bound = ', &
            require_topology_contribution_bound
        write(manifest_unit,'(A,L1)') 'clip_resonance_classes = ', &
            clip_resonance_classes
        call write_resonance_harmonic_set(manifest_unit,harmonic_ierr)
        write(manifest_unit,'(A,I0)') 'harmonic_guard_status = ',harmonic_ierr
        write(manifest_unit,'(A,ES24.16)') 'max_geometric_gap = ', &
            maxval(ledger_gap_energy)
        write(manifest_unit,'(A,ES24.16)') 'max_geometric_gap_bound = ', &
            maxval(ledger_gap_bound_energy)
        write(manifest_unit,'(A,ES24.16)') 'max_contribution_error_bound = ', &
            maxval(ledger_contribution_bound_energy)
        write(manifest_unit,'(A,L1)') 'all_contribution_bounds_certified = ', &
            all(ledger_contribution_certified_energy)
        write(manifest_unit,'(A,I0)') 'total_class_evaluations = ', &
            sum(ledger_classes_energy)
        write(manifest_unit,'(A,I0)') 'total_harmonic_scans = ', &
            sum(ledger_harmonics_energy)
        write(manifest_unit,'(A,I0)') 'total_searched_zero_harmonics = ', &
            sum(ledger_searched_zero_energy)
        write(manifest_unit,'(A,I0)') 'total_roots = ',sum(ledger_roots_energy)
        write(manifest_unit,'(A,I0)') 'total_failures = ',sum(ledger_failures_energy)
        close(manifest_unit)
    end subroutine write_run_manifest
    !
    subroutine open_energy_outputs(idx, adaptive, unit_fromresp, unit_jperp)
        !
        implicit none
        !
        integer, intent(in) :: idx
        logical, intent(in) :: adaptive
        integer, intent(out) :: unit_fromresp, unit_jperp
        character(len=256) :: path
        !
        call energy_path('fort.31415.energy.', idx, path)
        open(newunit=resline_unit,file=trim(path),status='replace',action='write')
        call energy_path('fort.31415.diagnostics.energy.', idx, path)
        open(newunit=resline_diag_unit,file=trim(path),status='replace',action='write')
        write(resline_diag_unit,'(A)') &
            '# toten perpinv psiast m n Rst Zst psif psiast_minus_psif x iclass taub delphi p lambda bmod phi_elec'
        !
        unit_fromresp=-1
        if(adaptive) then
            call energy_path('subint_ofH0int_104_vsJperp_fromresp.dat.energy.', idx, path)
            open(newunit=unit_fromresp,file=trim(path),status='replace',action='write')
            call energy_path('subint_ofH0int_104_vsJperp_adapt.dat.energy.', idx, path)
        else
            call energy_path('subint_ofH0int_104_vsJperp_equi.dat.energy.', idx, path)
        endif
        open(newunit=unit_jperp,file=trim(path),status='replace',action='write')
        !
    end subroutine open_energy_outputs
    !
    subroutine merge_energy_files(prefix, final_name, count)
        !
        implicit none
        !
        character(len=*), intent(in) :: prefix, final_name
        integer, intent(in) :: count
        character(len=4096) :: line
        character(len=256) :: path
        integer :: idx,io,in_unit,out_unit
        !
        open(newunit=out_unit,file=trim(final_name),status='replace',action='write')
        do idx=ienerg_begin,count
            call energy_path(prefix,idx,path)
            open(newunit=in_unit,file=trim(path),status='old',action='read',iostat=io)
            if(io.ne.0) then
                call tee_message('missing energy output: '//trim(path))
                stop
            endif
            do
                read(in_unit,'(A)',iostat=io) line
                if(io.ne.0) exit
                write(out_unit,'(A)') trim(line)
            enddo
            close(in_unit,status='delete')
        enddo
        close(out_unit)
        !
    end subroutine merge_energy_files
    !
    subroutine merge_jperp_files(prefix, final_name, count, prefix_scalar)
        !
        implicit none
        !
        character(len=*), intent(in) :: prefix, final_name
        integer, intent(in) :: count
        logical, intent(in) :: prefix_scalar
        character(len=4096) :: line
        character(len=256) :: path
        integer :: idx,io,in_unit,out_unit
        double precision :: scalar_prefix
        double precision, dimension(:), allocatable :: mode_prefix,vals
        !
        allocate(mode_prefix(nmodes),vals(nmodes+2))
        mode_prefix=0.d0
        scalar_prefix=0.d0
        open(newunit=out_unit,file=trim(final_name),status='replace',action='write')
        do idx=ienerg_begin,count
            call energy_path(prefix,idx,path)
            open(newunit=in_unit,file=trim(path),status='old',action='read',iostat=io)
            if(io.ne.0) then
                call tee_message('missing energy output: '//trim(path))
                stop
            endif
            do
                read(in_unit,'(A)',iostat=io) line
                if(io.ne.0) exit
                if(len_trim(line).eq.0) then
                    write(out_unit,*) ' '
                else
                    read(line,*) vals
                    if(prefix_scalar) vals(2)=vals(2)+scalar_prefix
                    vals(3:nmodes+2)=vals(3:nmodes+2)+mode_prefix
                    write(out_unit,*) vals
                endif
            enddo
            close(in_unit,status='delete')
            scalar_prefix=scalar_prefix+torque_int_energy(idx)
            mode_prefix=mode_prefix+torque_int_modes_energy(:,idx)
        enddo
        close(out_unit)
        deallocate(mode_prefix,vals)
        !
    end subroutine merge_jperp_files
    !
end subroutine resonant_torque
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
