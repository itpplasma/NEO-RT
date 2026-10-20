module resonance_mode_bounds_mod
    implicit none

contains

    pure function resonant_delphi_bound(m_modes, n_modes) result(bound)
        integer, intent(in) :: m_modes(:), n_modes(:)
        double precision :: bound
        double precision, parameter :: pi = 3.14159265358979d0

        ! Keep signed n in the resonance equation and Fourier phase; only this
        ! symmetric search extent uses its magnitude.
        bound = 2.d0*pi*(maxval(abs(dble(m_modes))/abs(dble(n_modes))) &
            + 1.d0/dble(minval(abs(n_modes))))
    end function resonant_delphi_bound

    pure logical function canonical_flux_outside_lcfs(psi_star, psi_axis, psi_edge)
        double precision, intent(in) :: psi_star, psi_axis, psi_edge

        ! All fluxes share the field gauge. Outside means beyond the edge in
        ! the axis-to-edge direction, independently of a common flux offset.
        canonical_flux_outside_lcfs = &
            (psi_star - psi_edge)*(psi_edge - psi_axis) > 0.d0
    end function canonical_flux_outside_lcfs

end module resonance_mode_bounds_mod
