module resonance_mode_bounds_mod
    !! Exact finite-harmonic guard for the optional POTATO class clip.
    !!
    !! For each signed pair (m,n), the physical resonance is
    !! g = Delta_phi + 2*pi*m/n.  The root search keeps this signed g.  The
    !! only magnitude used by the optional endpoint clip is the generated
    !! reverse-triangle bound |Delta_phi| > |2*pi*m/n|, which is a certified
    !! no-root condition for that harmonic.  The global clip is the generated
    !! finite-set maximum over the exact supplied (m,n) list.
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use potato_symbolic_kernel_mod, only : potato_resonance_harmonic_kernel, &
        potato_resonance_extent_envelope_kernel
    implicit none

    integer, parameter :: harmonic_guard_success=0
    integer, parameter :: harmonic_guard_invalid_set=1
    integer, parameter :: harmonic_guard_invalid_value=2
    integer, parameter :: harmonic_guard_uninitialized=3

    integer, allocatable, save :: guard_m(:),guard_n(:)
    double precision, allocatable, save :: guard_target(:),guard_extent(:)
    logical, save :: guard_initialized=.false.

contains

    subroutine set_resonance_harmonic_guard(m_modes,n_modes,bound,ierr)
        integer, intent(in) :: m_modes(:),n_modes(:)
        double precision, intent(out) :: bound
        integer, intent(out) :: ierr
        integer :: i
        double precision :: target,extent,resonance_g,no_root_margin
        double precision :: next_envelope,next_envelope_new

        ierr=harmonic_guard_invalid_set
        bound=0.d0
        guard_initialized=.false.
        if(allocated(guard_m)) deallocate(guard_m,guard_n,guard_target,guard_extent)
        if(size(m_modes).le.0 .or. size(m_modes).ne.size(n_modes)) return
        if(any(n_modes.eq.0)) return

        allocate(guard_m(size(m_modes)),guard_n(size(n_modes)), &
                 guard_target(size(m_modes)),guard_extent(size(m_modes)))
        guard_m=m_modes
        guard_n=n_modes
        next_envelope=0.d0
        do i=1,size(m_modes)
            call potato_resonance_harmonic_kernel(dble(m_modes(i)), &
                dble(n_modes(i)),0.d0,target,extent,resonance_g,no_root_margin)
            if(.not.valid_harmonic_values(target,extent,resonance_g, &
                                          no_root_margin)) then
                deallocate(guard_m,guard_n,guard_target,guard_extent)
                return
            endif
            guard_target(i)=target
            guard_extent(i)=extent
            call potato_resonance_extent_envelope_kernel(next_envelope,extent, &
                                                         next_envelope_new)
            next_envelope=next_envelope_new
            if(.not.ieee_is_finite(next_envelope) .or. next_envelope.lt.0.d0) then
                deallocate(guard_m,guard_n,guard_target,guard_extent)
                ierr=harmonic_guard_invalid_value
                return
            endif
        enddo
        bound=next_envelope
        guard_initialized=.true.
        ierr=harmonic_guard_success
    end subroutine set_resonance_harmonic_guard

    function resonant_delphi_bound(m_modes,n_modes) result(bound)
        integer, intent(in) :: m_modes(:),n_modes(:)
        double precision :: bound
        integer :: i
        double precision :: target,extent,resonance_g,no_root_margin
        double precision :: next_envelope,next_envelope_new

        bound=-1.d0
        if(size(m_modes).le.0 .or. size(m_modes).ne.size(n_modes)) return
        if(any(n_modes.eq.0)) return
        next_envelope=0.d0
        do i=1,size(m_modes)
            call potato_resonance_harmonic_kernel(dble(m_modes(i)), &
                dble(n_modes(i)),0.d0,target,extent,resonance_g,no_root_margin)
            if(.not.valid_harmonic_values(target,extent,resonance_g, &
                                          no_root_margin)) return
            call potato_resonance_extent_envelope_kernel(next_envelope,extent, &
                                                         next_envelope_new)
            next_envelope=next_envelope_new
            if(.not.ieee_is_finite(next_envelope) .or. next_envelope.lt.0.d0) then
                bound=-1.d0
                return
            endif
        enddo
        bound=next_envelope
    end function resonant_delphi_bound

    subroutine resonance_harmonic_values(m_mode,n_mode,delta_phi,target,extent, &
                                         resonance_g,no_root_margin,ierr)
        integer, intent(in) :: m_mode,n_mode
        double precision, intent(in) :: delta_phi
        double precision, intent(out) :: target,extent,resonance_g,no_root_margin
        integer, intent(out) :: ierr

        ierr=harmonic_guard_invalid_set
        target=0.d0
        extent=0.d0
        resonance_g=0.d0
        no_root_margin=0.d0
        if(n_mode.eq.0 .or. .not.ieee_is_finite(delta_phi)) return
        call potato_resonance_harmonic_kernel(dble(m_mode),dble(n_mode), &
            delta_phi,target,extent,resonance_g,no_root_margin)
        if(.not.valid_harmonic_values(target,extent,resonance_g,no_root_margin)) then
            ierr=harmonic_guard_invalid_value
            return
        endif
        ierr=harmonic_guard_success
    end subroutine resonance_harmonic_values

    subroutine resonance_no_root_for_any(delta_phi,no_root,ierr)
        double precision, intent(in) :: delta_phi
        logical, intent(out) :: no_root
        integer, intent(out) :: ierr
        integer :: i,local_ierr
        double precision :: target,extent,resonance_g,no_root_margin

        no_root=.false.
        ierr=harmonic_guard_uninitialized
        if(.not.guard_initialized) return
        if(.not.ieee_is_finite(delta_phi)) then
            ierr=harmonic_guard_invalid_value
            return
        endif
        no_root=.true.
        do i=1,size(guard_m)
            call resonance_harmonic_values(guard_m(i),guard_n(i),delta_phi, &
                target,extent,resonance_g,no_root_margin,local_ierr)
            if(local_ierr.ne.harmonic_guard_success) then
                no_root=.false.
                ierr=local_ierr
                return
            endif
            ! Equality is retained: a physical root can sit exactly at the
            ! clip boundary and must not be discarded by a strict guard.
            if(no_root_margin.le.0.d0) no_root=.false.
        enddo
        ierr=harmonic_guard_success
    end subroutine resonance_no_root_for_any

    subroutine write_resonance_harmonic_set(unit,ierr)
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        integer :: i

        ierr=harmonic_guard_uninitialized
        if(.not.guard_initialized) return
        write(unit,'(A,I0)') '# harmonic_set_count = ',size(guard_m)
        write(unit,'(A)') '# harmonic_index m n target_delphi search_extent'
        do i=1,size(guard_m)
            write(unit,'(I8,2I8,2ES24.16)') i,guard_m(i),guard_n(i), &
                guard_target(i),guard_extent(i)
        enddo
        ierr=harmonic_guard_success
    end subroutine write_resonance_harmonic_set

    pure logical function valid_harmonic_values(target,extent,resonance_g, &
                                                no_root_margin)
        double precision, intent(in) :: target,extent,resonance_g,no_root_margin

        valid_harmonic_values=ieee_is_finite(target) .and. &
            ieee_is_finite(extent) .and. ieee_is_finite(resonance_g) .and. &
            ieee_is_finite(no_root_margin) .and. extent.ge.0.d0
    end function valid_harmonic_values

    pure logical function canonical_flux_outside_lcfs(psi_star, psi_axis, psi_edge)
        double precision, intent(in) :: psi_star, psi_axis, psi_edge

        ! Outside means beyond the edge in the axis-to-edge flux direction.
        canonical_flux_outside_lcfs = &
            (psi_star - psi_edge)*(psi_edge - psi_axis) > 0.d0
    end function canonical_flux_outside_lcfs

end module resonance_mode_bounds_mod
