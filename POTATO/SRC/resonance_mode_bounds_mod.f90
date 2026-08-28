module resonance_mode_bounds_mod
    implicit none

contains

    pure function resonant_delphi_bound(m_modes, n_modes) result(bound)
        integer, intent(in) :: m_modes(:), n_modes(:)
        double precision :: bound
        double precision, parameter :: pi = 3.14159265358979d0

        ! The signed toroidal mode remains in m*Omega_b+n*Omega_phi=0 and in the
        ! Fourier phase.  Only this symmetric search extent is a magnitude.  Using
        ! signed n here made the native n<0 case produce a negative extent and
        ! silently discard every otherwise valid resonance contribution.
        bound = 2.d0*pi*(maxval(abs(dble(m_modes))/abs(dble(n_modes))) &
            + 1.d0/dble(minval(abs(n_modes))))
    end function resonant_delphi_bound

    pure logical function canonical_flux_outside_lcfs(psi_star, psi_axis, psi_edge)
        double precision, intent(in) :: psi_star, psi_axis, psi_edge

        ! Outside means beyond the edge in the axis-to-edge flux direction.  The
        ! earlier psi_star < psi_edge test was valid only when psi decreased from
        ! axis to edge and rejected every ITER resonance after field_eq selected
        ! the opposite native GEQDSK gauge.
        canonical_flux_outside_lcfs = &
            (psi_star - psi_edge)*(psi_edge - psi_axis) > 0.d0
    end function canonical_flux_outside_lcfs

end module resonance_mode_bounds_mod
