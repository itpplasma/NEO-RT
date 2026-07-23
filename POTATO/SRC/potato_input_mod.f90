module potato_input_mod
    use field_sub, only: read_field_input
    use input_files, only: convexfile
    use iso_fortran_env, only: output_unit
    use wall_loss_mod, only: load_wall
    implicit none

    ! Calculation type: 1=orbits, 2=equilibrium profiles, 3=resonant torque
    integer :: itest_type = 3

    ! Species parameters
    double precision :: E_alpha = 5d3       ! Reference energy [eV]
    double precision :: A_alpha = 2d0       ! Mass number
    double precision :: Z_alpha = 1d0       ! Charge number

    ! Flux surface
    double precision :: rho_pol = 0.6d0     ! Poloidal radius
    double precision :: rho_pol_max = 0.9d0 ! Max rho_pol for Poincare cut

    ! Scaling factors
    double precision :: scalfac_energy = 1d0
    double precision :: scalfac_efield = 1d0

    ! Orbit integration
    double precision :: Rmax_orbit = 200d0  ! Max R for orbit integration
    integer :: ntimstep = 30                ! Number of time steps

    ! Poincare cut
    integer :: npoicut = 10000

    ! Mode numbers for resonant torque
    integer :: m_min = -5
    integer :: m_max = 5
    integer :: n_tor = 1

    ! Energy grid
    integer :: nenerg = 60
    double precision :: thermen_max = 6d0   ! Max kinetic energy [T units]
    ! Lowest slice starts at this kinetic energy [T units] when set.
    double precision :: enkin_min_over_temp = 0d0

    ! Box counting
    integer :: nbox = 100

    ! Adaptive J_perp integration
    logical :: adaptive_jperp = .true.
    integer :: npoi_init = 51
    integer :: nlagr_sampling = 3
    double precision :: eps_sampling = 1d-2
    integer :: itermax_sampling = 5
    logical :: clip_resonance_classes = .true.

    ! Orbit plotting
    double precision :: toten_plot = 1d0
    double precision :: perpinv_plot = 4.5d-5

    ! Single-orbit trace (itest_type=4): start point on the poloidal plane and
    ! pitch cosine of one guiding-center orbit, traced over one bounce period.
    double precision :: orbit_Rstart = 210d0   ! start major radius [cm]
    double precision :: orbit_Zstart = 0d0     ! start height [cm]
    double precision :: orbit_lambda = 0.4d0   ! pitch cosine v_par/v at start

    ! Frequency radial scan (itest_type=5): trace one bounce at each of freq_n
    ! midplane start radii from freq_Rmin to freq_Rmax (pitch orbit_lambda),
    ! and write rho_pol, omega_b, omega_phi to freq_scan.dat.
    double precision :: freq_Rmin = 178d0      ! inner start radius [cm]
    double precision :: freq_Rmax = 221d0      ! outer start radius [cm], into SOL
    integer          :: freq_n = 40            ! number of surfaces

    ! Resonance probe diagnostic
    double precision :: probe_rho_pol = 0.9d0
    double precision :: probe_ux = 1.5d0
    double precision :: probe_eta = 4.1d-5
    double precision :: probe_toten = -1.0d0
    double precision :: probe_perpinv = -1.0d0
    integer :: probe_m = 0
    integer :: probe_n = 2

    ! Canonical frequency plot
    double precision :: enkin_over_temp = 1d0

    ! Profile input file
    character(len=256) :: profile_file = 'profile_poly.in'

    ! Edge extension: when .true. the orbit integrator accepts guiding-center
    ! points outside the separatrix (scrape-off layer), so bananas whose tips
    ! cross the last closed flux surface can be traced in the extended
    ! equilibrium field. Sets field_eq_mod::allow_sol. Default .false.
    logical :: edge_extension = .false.

    ! Accepted for old case files.
    logical :: plot_poicut = .false.
    logical :: plot_equilibrium = .false.

    namelist /potato_nml/ &
        itest_type, E_alpha, A_alpha, Z_alpha, &
        rho_pol, rho_pol_max, scalfac_energy, scalfac_efield, &
        Rmax_orbit, ntimstep, npoicut, plot_poicut, plot_equilibrium, &
        m_min, m_max, n_tor, &
        nenerg, thermen_max, enkin_min_over_temp, nbox, &
        adaptive_jperp, npoi_init, nlagr_sampling, eps_sampling, &
        itermax_sampling, clip_resonance_classes, &
        toten_plot, perpinv_plot, enkin_over_temp, &
        profile_file, edge_extension, &
        orbit_Rstart, orbit_Zstart, orbit_lambda, &
        freq_Rmin, freq_Rmax, freq_n, &
        probe_rho_pol, probe_ux, probe_eta, probe_toten, probe_perpinv, &
        probe_m, probe_n

contains

    subroutine read_potato_input(filename)
        character(len=*), intent(in) :: filename
        integer :: iunit, ierr
        logical :: field_input_exists, file_exists

        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, 'WARNING: Input file ' // trim(filename) &
                // ' not found, using defaults'
            return
        endif

        open(newunit=iunit, file=filename, status='old', &
            action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, 'WARNING: Cannot open ' // trim(filename) &
                // ', using defaults'
            return
        endif

        read(iunit, nml=potato_nml, iostat=ierr)
        if (ierr /= 0) then
            print *, 'WARNING: Error reading namelist from ' &
                // trim(filename) // ', some defaults may apply'
        endif

        close(iunit)

        if (.not. potato_input_is_valid()) then
            error stop 'invalid POTATO input: resonant-torque rho_pol_max must ' &
                // 'extend strictly beyond rho_pol and remain below one'
        endif

        inquire(file='field_divB0.inp', exist=field_input_exists)
        if (field_input_exists) then
            call read_field_input('field_divB0.inp')
            call load_wall(trim(convexfile))
        else
            call load_wall('convexwall.dat')
        endif
    end subroutine read_potato_input

    logical function potato_input_is_valid()
        potato_input_is_valid = 0.d0 < rho_pol .and. rho_pol <= rho_pol_max &
            .and. rho_pol_max < 1.d0
        if (itest_type == 3) then
            potato_input_is_valid = potato_input_is_valid &
                .and. rho_pol < rho_pol_max
        endif
    end function potato_input_is_valid

    subroutine print_potato_input(iunit)
        integer, intent(in) :: iunit

        write(iunit, '(A)') '=== POTATO Input Parameters ==='
        write(iunit, '(A,I0)') '  itest_type       = ', itest_type
        write(iunit, '(A,ES12.5)') '  E_alpha          = ', E_alpha
        write(iunit, '(A,ES12.5)') '  A_alpha          = ', A_alpha
        write(iunit, '(A,ES12.5)') '  Z_alpha          = ', Z_alpha
        write(iunit, '(A,ES12.5)') '  rho_pol          = ', rho_pol
        write(iunit, '(A,ES12.5)') '  rho_pol_max      = ', rho_pol_max
        write(iunit, '(A,ES12.5)') '  scalfac_energy   = ', scalfac_energy
        write(iunit, '(A,ES12.5)') '  scalfac_efield   = ', scalfac_efield
        write(iunit, '(A,ES12.5)') '  Rmax_orbit       = ', Rmax_orbit
        write(iunit, '(A,I0)') '  ntimstep         = ', ntimstep
        write(iunit, '(A,I0)') '  npoicut          = ', npoicut
        write(iunit, '(A,I0)') '  m_min            = ', m_min
        write(iunit, '(A,I0)') '  m_max            = ', m_max
        write(iunit, '(A,I0)') '  n_tor            = ', n_tor
        write(iunit, '(A,I0)') '  nenerg           = ', nenerg
        write(iunit, '(A,ES12.5)') '  thermen_max      = ', thermen_max
        write(iunit, '(A,ES12.5)') '  enkin_min_over_temp = ', enkin_min_over_temp
        write(iunit, '(A,I0)') '  nbox             = ', nbox
        write(iunit, '(A,L1)') '  adaptive_jperp   = ', adaptive_jperp
        write(iunit, '(A,I0)') '  npoi_init        = ', npoi_init
        write(iunit, '(A,I0)') '  nlagr_sampling   = ', nlagr_sampling
        write(iunit, '(A,ES12.5)') '  eps_sampling     = ', eps_sampling
        write(iunit, '(A,I0)') '  itermax_sampling = ', itermax_sampling
        write(iunit, '(A,L1)') '  clip_resonance_classes = ', clip_resonance_classes
        write(iunit, '(A,ES12.5)') '  toten_plot       = ', toten_plot
        write(iunit, '(A,ES12.5)') '  perpinv_plot     = ', perpinv_plot
        write(iunit, '(A,ES12.5)') '  enkin_over_temp  = ', enkin_over_temp
        write(iunit, '(A,A)') '  profile_file     = ', trim(profile_file)
        write(iunit, '(A,L1)') '  edge_extension   = ', edge_extension
        write(iunit, '(A,ES12.5)') '  orbit_Rstart     = ', orbit_Rstart
        write(iunit, '(A,ES12.5)') '  orbit_Zstart     = ', orbit_Zstart
        write(iunit, '(A,ES12.5)') '  orbit_lambda     = ', orbit_lambda
        write(iunit, '(A,ES12.5)') '  freq_Rmin        = ', freq_Rmin
        write(iunit, '(A,ES12.5)') '  freq_Rmax        = ', freq_Rmax
        write(iunit, '(A,I0)') '  freq_n           = ', freq_n
        write(iunit, '(A,ES12.5)') '  probe_rho_pol    = ', probe_rho_pol
        write(iunit, '(A,ES12.5)') '  probe_ux         = ', probe_ux
        write(iunit, '(A,ES12.5)') '  probe_eta        = ', probe_eta
        write(iunit, '(A,ES12.5)') '  probe_toten      = ', probe_toten
        write(iunit, '(A,ES12.5)') '  probe_perpinv    = ', probe_perpinv
        write(iunit, '(A,I0)') '  probe_m          = ', probe_m
        write(iunit, '(A,I0)') '  probe_n          = ', probe_n
        write(iunit, '(A)') '================================'
    end subroutine print_potato_input

end module potato_input_mod
