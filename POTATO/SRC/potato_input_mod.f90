module potato_input_mod
    use iso_fortran_env, only: output_unit
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

    ! Box counting
    integer :: nbox = 100

    ! Adaptive J_perp integration
    logical :: adaptive_jperp = .true.
    integer :: npoi_init = 51
    integer :: nlagr_sampling = 3
    double precision :: eps_sampling = 1d-2
    integer :: itermax_sampling = 5

    ! Orbit plotting
    double precision :: toten_plot = 1d0
    double precision :: perpinv_plot = 4.5d-5

    ! Canonical frequency plot
    double precision :: enkin_over_temp = 1d0

    ! Profile input file
    character(len=256) :: profile_file = 'profile_poly.in'

    namelist /potato_nml/ &
        itest_type, E_alpha, A_alpha, Z_alpha, &
        rho_pol, rho_pol_max, scalfac_energy, scalfac_efield, &
        Rmax_orbit, ntimstep, npoicut, &
        m_min, m_max, n_tor, &
        nenerg, thermen_max, nbox, &
        adaptive_jperp, npoi_init, nlagr_sampling, eps_sampling, &
        itermax_sampling, &
        toten_plot, perpinv_plot, enkin_over_temp, &
        profile_file

contains

    subroutine read_potato_input(filename)
        character(len=*), intent(in) :: filename
        integer :: iunit, ierr
        logical :: file_exists

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
    end subroutine read_potato_input

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
        write(iunit, '(A,I0)') '  nbox             = ', nbox
        write(iunit, '(A,L1)') '  adaptive_jperp   = ', adaptive_jperp
        write(iunit, '(A,I0)') '  npoi_init        = ', npoi_init
        write(iunit, '(A,I0)') '  nlagr_sampling   = ', nlagr_sampling
        write(iunit, '(A,ES12.5)') '  eps_sampling     = ', eps_sampling
        write(iunit, '(A,I0)') '  itermax_sampling = ', itermax_sampling
        write(iunit, '(A,ES12.5)') '  toten_plot       = ', toten_plot
        write(iunit, '(A,ES12.5)') '  perpinv_plot     = ', perpinv_plot
        write(iunit, '(A,ES12.5)') '  enkin_over_temp  = ', enkin_over_temp
        write(iunit, '(A,A)') '  profile_file     = ', trim(profile_file)
        write(iunit, '(A)') '================================'
    end subroutine print_potato_input

end module potato_input_mod
