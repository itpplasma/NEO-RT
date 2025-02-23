module neort_profiles
    use iso_fortran_env, only: dp => real64
    use polylag_3, only: plag1d, indef, mp
    use util, only: ev, qi, mi, mu, qe

    implicit none

    real(dp) vth, dvthds, M_t, dM_tds

    real(8) :: Om_tE, dOm_tEds, efac
    real(dp) ni1, ni2, Ti1, Ti2, Te, dni1ds, dni2ds, dTi1ds, dTi2ds, dTeds

    contains

    subroutine init_profile_input(s)
        ! Init s profile for finite orbit width boxes in radial s
        use do_magfie_mod, only: bfac, orbit_mode_transp
        use neort_profiles, only: M_t, dM_tds
        use spline

        ! For splining electric precession frequency
        real(8), allocatable :: Mt_spl_coeff(:, :)

        real(8), allocatable :: data(:, :)
        real(8) :: splineval(3)

        call readdata("profile.in", 3, data)

        allocate (Mt_spl_coeff(size(data(:, 1)), 5))

        Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
        splineval = spline_val_0(Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        if (orbit_mode_transp > 0) then
            call init_profile_finite_width(data)
        end if
    end subroutine init_profile_input

    subroutine init_plasma_input(s)
        real(dp), intent(in) :: s

        real(dp), parameter :: pmass = 1.6726d-24

        real(dp) :: amb, am1, am2, Zb, Z1, Z2, dchichi, slowrate, dchichi_norm, slowrate_norm
        real(dp) :: v0, ebeam
        real(dp), dimension(:, :), allocatable :: plasma(:, :)
        integer, dimension(mp) :: indu
        real(dp), dimension(mp) :: xp, fp
        real(dp) :: dxm1
        integer :: nplasma, i

        ! read plasma file
        open (1, file="plasma.in")
        read (1, *)
        read (1, *) nplasma, am1, am2, Z1, Z2
        read (1, *)
        allocate (plasma(nplasma, 6))
        do i = 1, nplasma
            read (1, *) plasma(i, :)
        end do
        dxm1 = 1.d0/(plasma(2, 1) - plasma(1, 1))
        close (1)

        ! interpolate to s value
        call indef(s, plasma(1, 1), dxm1, nplasma, indu)

        xp = plasma(indu, 1)
        fp = plasma(indu, 2)
        call plag1d(s, fp, dxm1, xp, ni1, dni1ds)
        fp = plasma(indu, 3)
        call plag1d(s, fp, dxm1, xp, ni2, dni2ds)
        fp = plasma(indu, 4)
        call plag1d(s, fp, dxm1, xp, Ti1, dTi1ds)
        fp = plasma(indu, 5)
        call plag1d(s, fp, dxm1, xp, Ti2, dTi2ds)
        fp = plasma(indu, 6)
        call plag1d(s, fp, dxm1, xp, Te, dTeds)

        qi = Z1*qe
        mi = am1*mu
        vth = sqrt(2d0*Ti1*ev/mi)
        dvthds = 0.5d0*sqrt(2d0*ev/(mi*Ti1))*dTi1ds
        v0 = vth
        amb = 2d0
        Zb = 1d0
        ebeam = amb*pmass*v0**2/(2d0*ev)

        call loacol_nbi(amb, am1, am2, Zb, Z1, Z2, ni1, ni2, Ti1, Ti2, Te, &
                        ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)

    end subroutine init_plasma_input

    subroutine init_profiles(use_thermodynamic_profiles)
        logical, intent(in) :: use_thermodynamic_profiles

        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0

        if (use_thermodynamic_profiles) then
            A1 = dni1ds/ni1 - qi/(Ti1*ev)*psi_pr/(q*c)*Om_tE - 3d0/2d0*dTi1ds/Ti1
            A2 = dTi1ds/Ti1
        end if
    end subroutine init_profiles

    subroutine init_profile_finite_width(data)
        use driftorbit, only: sbox, taubins, fluxint_box, torque_int_box

        real(8), allocatable, intent(in) :: data(:, :)

        integer :: k

        allocate (sbox(size(data, 1) + 1))
        allocate (taubins(size(sbox) + 1))
        sbox(1) = 0.0
        sbox(size(sbox)) = 1.0
        do k = 1, (size(data, 1) - 1)
            sbox(k + 1) = (data(k, 1) + data(k + 1, 1))/2d0
        end do
        allocate (fluxint_box(2, size(sbox) + 1))
        allocate (torque_int_box(size(sbox) + 1))
    end subroutine init_profile_finite_width

end module neort_profiles
