program test_single_well_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_config, only: config_t
    use neort_lib, only: neort_init, neort_compute_no_splines
    use neort_datatypes, only: transport_data_t
    use neort_lorentz_native, only: validate_lorentz_geometry, sample_lorentz_orbit
    use neort_profiles, only: vth
    use driftorbit, only: etatp, etadt
    use util, only: pi, mi
    implicit none
    real(dp), parameter :: field_t = 2.0_dp, major_m = 5.0_dp, minor_m = 0.5_dp
    real(dp), parameter :: surface = 0.5_dp, rate = 2.0e5_dp
    real(dp), parameter :: radius = minor_m*sqrt(surface), epsilon = radius/major_m
    real(dp), parameter :: field_g = field_t*1.0e4_dp
    real(dp), parameter :: g_covar = field_g*major_m*100.0_dp
    real(dp), parameter :: i_covar = field_g*(radius*100.0_dp)**2/(200*major_m)
    real(dp) :: delta, x(96), w(96), eta, tau, diffusion, omega, action
    real(dp) :: expected_tau, expected_diffusion, energy, mu, step, action_plus
    real(dp) :: action_minus, relative, max_error, action_error
    complex(dp) :: drive, expected_drive
    type(config_t) :: config
    type(transport_data_t) :: unused
    character(64) :: field_path
    integer :: specimen, branch, info

    call gauss_legendre(size(x), x, w)
    call write_perturbation('single_well_pert.bc')
    config%s = surface
    config%M_t = 0.01_dp
    config%qs = 1.0_dp
    config%ms = 2.014_dp
    config%vth = 4.0e7_dp
    config%epsmn = 1.0e-3_dp
    config%m0 = -6
    config%mph = 3
    config%comptorque = .false.
    config%magdrift = .false.
    config%noshear = .true.
    config%inp_swi = 8
    config%inp_swi_pert = 9
    config%pertfile = .true.
    config%vsteps = 1
    config%mth_max_abs = 0
    max_error = 0.0_dp
    action_error = 0.0_dp
    do specimen = 1, 2
        delta = merge(0.1_dp, -0.1_dp, specimen == 1)*epsilon
        write(field_path, '(a,i0,a)') 'single_well_', specimen, '.bc'
        call write_chart(trim(field_path))
        call neort_init(config, trim(field_path), 'single_well_pert.bc')
        call neort_compute_no_splines(unused)
        call validate_lorentz_geometry(info)
        if (info /= 0) error stop 'certified two-harmonic single well rejected'
        do branch = 1, 3
            eta = merge((etatp + etadt)/2, 0.4_dp*etatp, branch == 1)
            call independent_orbit(vth, eta, expected_tau, expected_drive, &
                                   expected_diffusion, action)
            call sample_lorentz_orbit(1.0_dp, branch, eta, rate, 0.0_dp, &
                                     tau, diffusion, omega, drive)
            relative = max(abs(tau/expected_tau - 1.0_dp), &
                           abs(drive - expected_drive)/abs(expected_drive), &
                           abs(diffusion/expected_diffusion - 1.0_dp))
            max_error = max(max_error, relative)
            if (relative > 3.0e-7_dp) then
                print *, 'two-harmonic orbit error:', specimen, branch, relative
                error stop 'native period/complex Hamiltonian/diffusion oracle failed'
            end if
            if (abs(aimag(drive)) < 0.1_dp*abs(drive)) &
                error stop 'complex source oracle lost its independent phase'
            ! Differentiate the independently integrated positive action at
            ! fixed magnetic moment mu=eta*E, not at fixed pitch eta.
            energy = mi*vth**2/2
            mu = eta*energy
            step = 1.0e-4_dp*energy
            call independent_orbit(sqrt(2*(energy + step)/mi), &
                    mu/(energy + step), tau, drive, diffusion, action_plus)
            call independent_orbit(sqrt(2*(energy - step)/mi), &
                    mu/(energy - step), tau, drive, diffusion, action_minus)
            relative = abs(pi*(action_plus - action_minus)/step/expected_tau - 1)
            action_error = max(action_error, relative)
            if (relative > 3.0e-7_dp) error stop '2pi*dJ/dE differs from period'
        end do
    end do
    print *, 'two-harmonic native orbit and action errors:', max_error, action_error
    print *, 'test_single_well_native PASSED'

contains

    pure function magnetic_field(theta) result(bmod)
        real(dp), intent(in) :: theta
        real(dp) :: bmod
        bmod = field_g*(1 - epsilon*cos(theta) + delta*cos(2*theta))
    end function magnetic_field

    subroutine independent_orbit(speed, pitch, period, hamiltonian, diff, action_)
        real(dp), intent(in) :: speed, pitch
        real(dp), intent(out) :: period, diff, action_
        complex(dp), intent(out) :: hamiltonian
        real(dp) :: lo, hi, turn, alpha, theta, jacobian, bmod, htheta
        real(dp) :: weight, parallel_square, inverse_integral
        complex(dp) :: local_h
        integer :: point, iteration
        logical :: trapped

        trapped = pitch*magnetic_field(pi) > 1.0_dp
        turn = pi
        if (trapped) then
            lo = 0.0_dp
            hi = pi
            do iteration = 1, 60
                turn = (lo + hi)/2
                if (pitch*magnetic_field(turn) > 1.0_dp) then
                    hi = turn
                else
                    lo = turn
                end if
            end do
            turn = (lo + hi)/2
        end if
        period = 0.0_dp
        hamiltonian = (0.0_dp, 0.0_dp)
        inverse_integral = 0.0_dp
        action_ = 0.0_dp
        do point = 1, size(x)
            if (trapped) then
                ! theta=turn*sin(alpha) regularizes either simple turn for
                ! this monotone well; no pure-cosine inverse is used.
                alpha = pi*x(point)/2
                theta = turn*sin(alpha)
                jacobian = turn*cos(alpha)
            else
                theta = pi*(1 + x(point))
                jacobian = 1.0_dp
            end if
            bmod = magnetic_field(theta)
            htheta = bmod/(2*g_covar + i_covar)
            parallel_square = 1 - pitch*bmod
            weight = pi*w(point)*jacobian/(speed*htheta*sqrt(parallel_square))
            ! The trapped integral includes both visits to [-turn,+turn].
            ! The same pi factor is the passing [0,2pi] change of variable.
            local_h = mi*speed**2/2*(2 - pitch*bmod)*1.0e-3_dp*field_g/bmod* &
                (cmplx(1.0_dp, 0.3_dp, dp) + cmplx(0.2_dp, -0.1_dp, dp)* &
                 exp(cmplx(0.0_dp, theta, dp)))
            period = period + weight
            hamiltonian = hamiltonian + weight*local_h
            inverse_integral = inverse_integral + weight*(1/bmod - pitch)
            action_ = action_ + weight*mi*speed**2*parallel_square/(2*pi)
        end do
        hamiltonian = hamiltonian/period
        diff = 2*rate*pitch*inverse_integral
    end subroutine independent_orbit

    subroutine write_chart(path)
        character(*), intent(in) :: path
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]
        real(dp) :: local_radius, bph, bth
        integer :: unit, k
        open(newunit=unit, file=path, status='replace')
        write(unit, '(a)') 'CC manufactured analytic two-harmonic single well'
        write(unit, '(a)') 'CC constant q=2, epsilon>4*abs(delta), min theta=0'
        write(unit, '(a)') 'CC m0b n0b nsurf nper flux a R'
        write(unit, '(a)') 'CC'
        write(unit, '(a)') 'CC'
        write(unit, *) 2, 0, 3, 1, -pi*minor_m**2*field_t, minor_m, major_m
        do k = 1, size(surfaces)
            local_radius = minor_m*sqrt(surfaces(k))
            bph = -field_t*major_m*1.0e6_dp/0.2_dp
            bth = -field_t*local_radius**2*0.5_dp/major_m*1.0e6_dp/0.2_dp
            write(unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write(unit, '(a)') 'CC units A A Pa m3'
            write(unit, *) surfaces(k), 0.5_dp, bph, bth, 0.0_dp, 0.0_dp
            write(unit, '(a)') 'CC m n rmnc zmnc vmnc bmnc'
            write(unit, *) 0, 0, major_m, 0.0_dp, 0.0_dp, field_t
            write(unit, *) 1, 0, local_radius, local_radius, 0.0_dp, &
                           -field_t*local_radius/major_m
            write(unit, *) 2, 0, 0.0_dp, 0.0_dp, 0.0_dp, field_t*delta
        end do
        close(unit)
    end subroutine write_chart

    subroutine write_perturbation(path)
        character(*), intent(in) :: path
        real(dp), parameter :: surfaces(3) = [0.25_dp, 0.5_dp, 0.75_dp]
        integer :: unit, k
        open(newunit=unit, file=path, status='replace')
        write(unit, '(a)') 'CC independent complex perturbation'
        write(unit, '(a)') 'CC complex m=-6,-5 source independent of the well shape'
        write(unit, '(a)') 'CC'
        write(unit, '(a)') 'CC'
        write(unit, '(a)') 'CC'
        write(unit, *) 1, 0, 3, 1, 0.0_dp, 0.0_dp, 0.0_dp
        do k = 1, size(surfaces)
            write(unit, '(a)') 'CC s iota Jpol Itor pprime sqrtg00'
            write(unit, '(a)') 'CC'
            write(unit, *) surfaces(k), 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write(unit, '(a)') 'CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write(unit, *) -6, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, field_t*1.0e-3_dp, -field_t*3.0e-4_dp
            write(unit, *) -5, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, field_t*2.0e-4_dp, field_t*1.0e-4_dp
        end do
        close(unit)
    end subroutine write_perturbation
end program test_single_well_native
