module potato_invariants_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    private

    integer, parameter, public :: POTATO_INVARIANT_SUCCESS = 0
    integer, parameter, public :: POTATO_INVARIANT_INVALID_INPUT = 1
    integer, parameter, public :: POTATO_INVARIANT_NONZERO_POTENTIAL = 2
    integer, parameter, public :: POTATO_INVARIANT_FIELD_ERROR = 3
    real(dp), parameter :: POTENTIAL_TOLERANCE = 1.0e-12_dp

    type, public :: potato_invariant_t
        real(dp) :: p = 0.0_dp
        real(dp) :: xi = 0.0_dp
        real(dp) :: h0 = 0.0_dp
        real(dp) :: j_perp = 0.0_dp
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: psi_axis = 0.0_dp
        real(dp) :: psi_edge = 0.0_dp
        real(dp) :: phi_elec = 0.0_dp
        real(dp) :: v0 = 0.0_dp
        integer :: sigma = 0
        integer :: status = POTATO_INVARIANT_INVALID_INPUT
    end type potato_invariant_t

    public :: open_potato_invariant_handoff, potato_invariant_from_local_state
    public :: write_potato_invariant_from_local_state

contains

    subroutine open_potato_invariant_handoff(iunit)
        integer, intent(out) :: iunit

        open(newunit=iunit, file='potato_invariants.dat', status='replace', &
            action='write')
        write(iunit, '(A)') '# potato-simple-invariants-v1'
        write(iunit, '(A)') '# id R_cm Z_cm p xi sigma H0 J_perp psi_star '// &
            'psi_axis psi_edge phi_elec v0_cm_s status'
    end subroutine open_potato_invariant_handoff

    subroutine write_potato_invariant_from_local_state(iunit, identifier, &
            radius, height, p, xi, bmod, h_phi, psi, phi_elec, rho0, &
            psi_axis, psi_edge, v0, field_status)
        integer, intent(in) :: iunit, identifier, field_status
        real(dp), intent(in) :: radius, height, p, xi, bmod, h_phi, psi
        real(dp), intent(in) :: phi_elec, rho0, psi_axis, psi_edge, v0
        type(potato_invariant_t) :: invariant

        call potato_invariant_from_local_state(p, xi, bmod, h_phi, psi, &
            phi_elec, rho0, psi_axis, psi_edge, v0, invariant)
        if (field_status /= 0) invariant%status = POTATO_INVARIANT_FIELD_ERROR
        call write_potato_invariant(iunit, identifier, radius, height, invariant)
    end subroutine write_potato_invariant_from_local_state

    subroutine write_potato_invariant(iunit, identifier, radius, height, invariant)
        integer, intent(in) :: iunit, identifier
        real(dp), intent(in) :: radius, height
        type(potato_invariant_t), intent(in) :: invariant

        write(iunit, '(I0,1X,4(ES24.16,1X),I0,1X,7(ES24.16,1X),I0)') &
            identifier, radius, height, invariant%p, invariant%xi, &
            invariant%sigma, invariant%h0, invariant%j_perp, &
            invariant%psi_star, invariant%psi_axis, invariant%psi_edge, &
            invariant%phi_elec, invariant%v0, invariant%status
    end subroutine write_potato_invariant

    pure subroutine potato_invariant_from_local_state(p, xi, bmod, h_phi, &
            psi, phi_elec, rho0, psi_axis, psi_edge, v0, invariant)
        real(dp), intent(in) :: p, xi, bmod, h_phi, psi, phi_elec, rho0
        real(dp), intent(in) :: psi_axis, psi_edge, v0
        type(potato_invariant_t), intent(out) :: invariant

        invariant = potato_invariant_t()
        call copy_input_metadata(p, xi, phi_elec, psi_axis, psi_edge, v0, &
            invariant)
        if (.not. local_state_is_valid(p, xi, bmod, h_phi, psi, rho0, &
                psi_axis, psi_edge, v0)) return
        invariant%h0 = p**2 + phi_elec
        invariant%j_perp = p**2*(1.0_dp - xi**2)/bmod
        invariant%psi_star = psi + rho0*p*xi*h_phi
        if (abs(phi_elec) > POTENTIAL_TOLERANCE) then
            invariant%status = POTATO_INVARIANT_NONZERO_POTENTIAL
            return
        end if
        invariant%status = POTATO_INVARIANT_SUCCESS
    end subroutine potato_invariant_from_local_state

    pure subroutine copy_input_metadata(p, xi, phi_elec, psi_axis, psi_edge, &
            v0, invariant)
        real(dp), intent(in) :: p, xi, phi_elec, psi_axis, psi_edge, v0
        type(potato_invariant_t), intent(inout) :: invariant

        invariant%p = p
        invariant%xi = xi
        invariant%phi_elec = phi_elec
        invariant%psi_axis = psi_axis
        invariant%psi_edge = psi_edge
        invariant%v0 = v0
        if (xi > 0.0_dp) invariant%sigma = 1
        if (xi < 0.0_dp) invariant%sigma = -1
    end subroutine copy_input_metadata

    pure logical function local_state_is_valid(p, xi, bmod, h_phi, psi, &
            rho0, psi_axis, psi_edge, v0)
        real(dp), intent(in) :: p, xi, bmod, h_phi, psi, rho0
        real(dp), intent(in) :: psi_axis, psi_edge, v0

        local_state_is_valid = .false.
        if (.not. ieee_is_finite(p)) return
        if (.not. ieee_is_finite(xi)) return
        if (.not. ieee_is_finite(bmod)) return
        if (.not. ieee_is_finite(h_phi)) return
        if (.not. ieee_is_finite(psi)) return
        if (.not. ieee_is_finite(rho0)) return
        if (.not. ieee_is_finite(psi_axis)) return
        if (.not. ieee_is_finite(psi_edge)) return
        if (.not. ieee_is_finite(v0)) return
        if (p <= 0.0_dp) return
        if (abs(xi) > 1.0_dp) return
        if (bmod <= 0.0_dp) return
        if (v0 <= 0.0_dp) return
        if (psi_edge == psi_axis) return
        local_state_is_valid = .true.
    end function local_state_is_valid
end module potato_invariants_mod
