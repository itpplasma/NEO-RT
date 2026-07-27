module box_counting_status_mod
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer, parameter :: box_count_skipped_zero_weight = -1
    integer, parameter :: box_count_ok = 0
    integer, parameter :: box_count_nonpositive_bounce = 1
    integer, parameter :: box_count_vode_mxstep = 2
    integer, parameter :: box_count_vode_failed = 3
    integer, parameter :: box_count_nonfinite = 4
    integer, parameter :: box_count_negative_residence = 5
    integer, parameter :: box_count_residence_mismatch = 6

contains

    pure integer function box_count_vode_status(istate)
        integer, intent(in) :: istate

        if (istate == 2) then
            box_count_vode_status = box_count_ok
        elseif (istate == -1) then
            box_count_vode_status = box_count_vode_mxstep
        else
            box_count_vode_status = box_count_vode_failed
        endif
    end function box_count_vode_status

    pure integer function box_count_residence_status(taub, tau, tau_outside)
        real(8), intent(in) :: taub
        real(8), intent(in) :: tau(:)
        real(8), intent(in) :: tau_outside

        real(8) :: residence_sum, tolerance

        if (.not. ieee_is_finite(taub) .or. &
            .not. all(ieee_is_finite(tau)) .or. &
            .not. ieee_is_finite(tau_outside)) then
            box_count_residence_status = box_count_nonfinite
            return
        endif
        if (taub <= 0d0) then
            box_count_residence_status = box_count_nonpositive_bounce
            return
        endif

        tolerance = max(1d-10*abs(taub), &
            1024d0*epsilon(taub)*max(abs(taub), 1d0))
        if (any(tau < -tolerance) .or. tau_outside < -tolerance) then
            box_count_residence_status = box_count_negative_residence
            return
        endif

        residence_sum = sum(tau) + tau_outside
        if (abs(residence_sum-taub) > tolerance) then
            box_count_residence_status = box_count_residence_mismatch
            return
        endif

        box_count_residence_status = box_count_ok
    end function box_count_residence_status

end module box_counting_status_mod

subroutine linspace(a, b, cnt, out)
    implicit none

    real(8), intent(in) :: a, b
    integer :: i, cnt
    real(8), intent(inout) :: out(cnt)
    real(8) :: delta

    delta = (b-a)/(cnt-1)
    do i=1,cnt
        out(i) = a + delta*(i-1)
    end do
end subroutine linspace

subroutine timestep_vode(n, tau, z, vz)
    implicit none
    ! See velo for details
    integer(4), intent(in) :: n ! number of equations
    real(8), intent(in)    :: tau ! Time
    real(8), intent(in)    :: z(n) ! Phase-position
    real(8), intent(out)   :: vz(n) ! Phase-velocity
    call velo(tau, z, vz)
end subroutine timestep_vode

subroutine time_in_box(z, cnt, sbox, taub, tau, tau_outside, status)
    ! Returns time spent in boxes
    use dvode_f90_m, only: vode_opts, set_opts, dvode_f90, get_stats
    use field_sub, only: psif
    use field_eq_mod, only: psi_axis, psi_sep
    use orbit_dim_mod, only: neqm
    use box_counting_status_mod, only: box_count_ok, &
        box_count_nonpositive_bounce, box_count_nonfinite, &
        box_count_residence_status, box_count_vode_status
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    real(8), intent(in)    :: z(neqm) ! Starting position
    integer, intent(in)    :: cnt ! Box boundary count
    real(8), intent(in)    :: sbox(cnt) ! Box boundaries
    real(8), intent(in)    :: taub ! Bounce time
    real(8), intent(out)   :: tau(cnt) ! Time in each box
    real(8), intent(out)   :: tau_outside ! Time at s_pol > final box edge
    integer, intent(out)   :: status ! box_counting_status_mod code

    real(8) :: smid
    real(8) :: y(neqm)

    integer(4) :: k, nsample
    real(8) :: ti

    real(8) :: atol(neqm), rtol, tout
    integer(4) :: itask, istate, method_flag
    type (vode_opts) :: options
        real(8) :: bmod, phi_elec, s
        real(8) :: sold, told
        integer(4) :: sind ! s index

        external timestep_vode

        tau = 0d0
        tau_outside = 0d0
        if (.not. ieee_is_finite(taub)) then
            status = box_count_nonfinite
            return
        elseif (taub <= 0d0) then
            ! An interpolant overshoot through zero is a failed disposition, not a
            ! zero torque deposition.
            status = box_count_nonpositive_bounce
            return
        else
            status = box_count_ok
        endif

        y = z

        rtol = 1d-12
        atol = 1d-13
        itask = 1
        istate = 1
        method_flag = 10
        options = set_opts(method_flag=method_flag, abserr_vector=atol, &
            relerr=rtol, mxstep=2000)

        nsample = min(256, max(64, size(sbox)))
        ti = 0d0

        call get_bmod_and_Phi(z(1:3), bmod, phi_elec)

        s = abs((psif-psi_axis)/(psi_sep-psi_axis))
        sold = s
        told = 0d0
        do k = 1,nsample
            tout = taub*dble(k)/dble(nsample)
            call dvode_f90( &
                timestep_vode, neqm, y, ti, tout, itask, istate, options)
            status = box_count_vode_status(istate)
            if (status /= box_count_ok) then
                if (istate == -1) then
                    print *, 'time_in_box: VODE exceeded mxstep at k =', k, &
                        ' ti =', ti, ' taub =', taub
                else
                    print *, 'time_in_box: VODE failed with istate =', istate, &
                        ' k =', k, ' ti =', ti, ' taub =', taub
                endif
                tau = 0d0
                tau_outside = 0d0
                return
            end if

            call get_bmod_and_Phi(y(1:3), bmod, phi_elec)

            s = abs((psif-psi_axis)/(psi_sep-psi_axis))
            smid = 0.5d0*(sold+s)
            if (.not. ieee_is_finite(smid)) then
                status = box_count_nonfinite
                tau = 0d0
                tau_outside = 0d0
                return
            elseif (smid > sbox(size(sbox))) then
                tau_outside = tau_outside + ti-told
            else
                sind = radial_box_index(smid)
                tau(sind) = tau(sind) + ti-told
            endif
            sold = s
            told = ti
        end do

        status = box_count_residence_status(taub, tau, tau_outside)
        if (status /= box_count_ok) then
            tau = 0d0
            tau_outside = 0d0
        endif

    contains

        integer function radial_box_index(sval)
            real(8), intent(in) :: sval
            integer :: ibox

            radial_box_index = size(sbox)
            do ibox = 1,size(sbox)
                if (sval <= sbox(ibox)) then
                    radial_box_index = ibox
                    exit
                endif
            enddo
        end function radial_box_index
    end subroutine time_in_box
