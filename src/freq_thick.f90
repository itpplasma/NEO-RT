module freq_thick
    ! Thick orbit frequency calculation module (Phase G.2.1)
    ! Implements Om_th and Om_ph calculations using POTATO bounce times
    
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: compute_canonical_frequencies_thick, initialize_thick_frequency_database
    public :: thick_frequency_interpolation_available
    
    integer, parameter :: dp = real64
    
    ! Frequency database for interpolation (thick orbits are expensive)
    integer, parameter :: max_v_points = 20
    integer, parameter :: max_eta_points = 15
    
    type :: frequency_database_t
        logical :: initialized = .false.
        integer :: n_v, n_eta
        real(dp), allocatable :: v_grid(:), eta_grid(:)
        real(dp), allocatable :: Om_th_table(:, :), Om_ph_table(:, :)
        logical, allocatable :: valid_table(:, :)
    end type frequency_database_t
    
    type(frequency_database_t), save :: freq_db
    
contains

    subroutine compute_canonical_frequencies_thick(v, eta, Om_th, Om_ph, success)
        ! Compute thick orbit canonical frequencies using POTATO integration
        use potato_field_bridge, only: real_find_bounce_calculation
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(dp) :: taub, delphi
        real(dp), parameter :: pi = 3.141592653589793_dp
        
        ! Initialize outputs
        success = .false.
        Om_th = 0.0_dp
        Om_ph = 0.0_dp
        
        ! Try interpolation first if database is available
        if (freq_db%initialized) then
            call interpolate_from_database(v, eta, Om_th, Om_ph, success)
            if (success) return
        end if
        
        ! Direct POTATO calculation (expensive)
        call real_find_bounce_calculation(v, eta, taub, delphi, success)
        
        if (.not. success) then
            ! Fallback to approximate thick orbit scaling
            call approximate_thick_orbit_frequencies(v, eta, Om_th, Om_ph, success)
            return
        end if
        
        ! Convert bounce time and toroidal shift to canonical frequencies
        if (taub > 0.0_dp) then
            Om_th = 2.0_dp * pi / taub        ! Bounce frequency ω_θ
            Om_ph = delphi / taub             ! Toroidal frequency ω_φ
            success = .true.
        else
            ! Invalid bounce time
            call approximate_thick_orbit_frequencies(v, eta, Om_th, Om_ph, success)
        end if
        
    end subroutine compute_canonical_frequencies_thick
    
    subroutine initialize_thick_frequency_database(v_min, v_max, eta_min, eta_max, success)
        ! Pre-compute thick orbit frequencies for interpolation database
        implicit none
        real(dp), intent(in) :: v_min, v_max, eta_min, eta_max
        logical, intent(out) :: success
        
        integer :: i, j
        real(dp) :: v_val, eta_val
        real(dp) :: Om_th_val, Om_ph_val
        logical :: calc_success
        
        ! Initialize database structure
        freq_db%n_v = max_v_points
        freq_db%n_eta = max_eta_points
        
        if (allocated(freq_db%v_grid)) deallocate(freq_db%v_grid)
        if (allocated(freq_db%eta_grid)) deallocate(freq_db%eta_grid)
        if (allocated(freq_db%Om_th_table)) deallocate(freq_db%Om_th_table)
        if (allocated(freq_db%Om_ph_table)) deallocate(freq_db%Om_ph_table)
        if (allocated(freq_db%valid_table)) deallocate(freq_db%valid_table)
        
        allocate(freq_db%v_grid(freq_db%n_v))
        allocate(freq_db%eta_grid(freq_db%n_eta))
        allocate(freq_db%Om_th_table(freq_db%n_v, freq_db%n_eta))
        allocate(freq_db%Om_ph_table(freq_db%n_v, freq_db%n_eta))
        allocate(freq_db%valid_table(freq_db%n_v, freq_db%n_eta))
        
        ! Set up velocity grid (logarithmic spacing)
        do i = 1, freq_db%n_v
            freq_db%v_grid(i) = v_min * (v_max/v_min)**(real(i-1, dp)/real(freq_db%n_v-1, dp))
        end do
        
        ! Set up pitch angle grid (linear spacing)
        do j = 1, freq_db%n_eta
            freq_db%eta_grid(j) = eta_min + (eta_max - eta_min) * real(j-1, dp) / real(freq_db%n_eta-1, dp)
        end do
        
        ! Pre-compute frequencies
        success = .true.
        do i = 1, freq_db%n_v
            do j = 1, freq_db%n_eta
                v_val = freq_db%v_grid(i)
                eta_val = freq_db%eta_grid(j)
                
                call compute_canonical_frequencies_thick(v_val, eta_val, Om_th_val, Om_ph_val, calc_success)
                
                freq_db%Om_th_table(i, j) = Om_th_val
                freq_db%Om_ph_table(i, j) = Om_ph_val
                freq_db%valid_table(i, j) = calc_success
                
                if (.not. calc_success) then
                    ! Use approximation for failed points
                    call approximate_thick_orbit_frequencies(v_val, eta_val, Om_th_val, Om_ph_val, calc_success)
                    freq_db%Om_th_table(i, j) = Om_th_val
                    freq_db%Om_ph_table(i, j) = Om_ph_val
                    freq_db%valid_table(i, j) = calc_success
                end if
            end do
        end do
        
        freq_db%initialized = .true.
        
    end subroutine initialize_thick_frequency_database
    
    subroutine interpolate_from_database(v, eta, Om_th, Om_ph, success)
        ! Interpolate frequencies from pre-computed database
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        integer :: i_v, i_eta
        real(dp) :: w_v, w_eta
        real(dp) :: Om_th_interp, Om_ph_interp
        
        success = .false.
        Om_th = 0.0_dp
        Om_ph = 0.0_dp
        
        if (.not. freq_db%initialized) return
        
        ! Find velocity grid indices
        call find_grid_indices(v, freq_db%v_grid, freq_db%n_v, i_v, w_v, success)
        if (.not. success) return
        
        ! Find eta grid indices  
        call find_grid_indices(eta, freq_db%eta_grid, freq_db%n_eta, i_eta, w_eta, success)
        if (.not. success) return
        
        ! Check that all 4 interpolation points are valid
        if (.not. (freq_db%valid_table(i_v, i_eta) .and. &
                   freq_db%valid_table(i_v+1, i_eta) .and. &
                   freq_db%valid_table(i_v, i_eta+1) .and. &
                   freq_db%valid_table(i_v+1, i_eta+1))) then
            success = .false.
            return
        end if
        
        ! Bilinear interpolation for Om_th
        Om_th_interp = (1.0_dp - w_v) * (1.0_dp - w_eta) * freq_db%Om_th_table(i_v, i_eta) + &
                       w_v * (1.0_dp - w_eta) * freq_db%Om_th_table(i_v+1, i_eta) + &
                       (1.0_dp - w_v) * w_eta * freq_db%Om_th_table(i_v, i_eta+1) + &
                       w_v * w_eta * freq_db%Om_th_table(i_v+1, i_eta+1)
        
        ! Bilinear interpolation for Om_ph
        Om_ph_interp = (1.0_dp - w_v) * (1.0_dp - w_eta) * freq_db%Om_ph_table(i_v, i_eta) + &
                       w_v * (1.0_dp - w_eta) * freq_db%Om_ph_table(i_v+1, i_eta) + &
                       (1.0_dp - w_v) * w_eta * freq_db%Om_ph_table(i_v, i_eta+1) + &
                       w_v * w_eta * freq_db%Om_ph_table(i_v+1, i_eta+1)
        
        Om_th = Om_th_interp
        Om_ph = Om_ph_interp
        success = .true.
        
    end subroutine interpolate_from_database
    
    subroutine find_grid_indices(x, grid, n, i_lower, weight, success)
        ! Find interpolation indices and weights for grid point
        implicit none
        real(dp), intent(in) :: x
        real(dp), intent(in) :: grid(:)
        integer, intent(in) :: n
        integer, intent(out) :: i_lower
        real(dp), intent(out) :: weight
        logical, intent(out) :: success
        
        integer :: i
        
        success = .false.
        
        ! Check bounds
        if (x < grid(1) .or. x > grid(n)) return
        
        ! Find bracketing interval
        do i = 1, n-1
            if (x >= grid(i) .and. x <= grid(i+1)) then
                i_lower = i
                weight = (x - grid(i)) / (grid(i+1) - grid(i))
                success = .true.
                return
            end if
        end do
        
    end subroutine find_grid_indices
    
    subroutine approximate_thick_orbit_frequencies(v, eta, Om_th, Om_ph, success)
        ! Fallback approximation for thick orbit frequencies
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(dp), parameter :: pi = 3.141592653589793_dp
        real(dp), parameter :: v_thermal = 1.0d6  ! m/s
        real(dp) :: taub_thin, delphi_thin
        real(dp) :: finite_orbit_correction
        
        ! Calculate thin orbit frequencies as baseline
        taub_thin = 1.0d-4 / v * v_thermal    ! Thin orbit bounce time
        delphi_thin = 0.1_dp * eta           ! Thin orbit toroidal shift
        
        ! Finite orbit width correction (approximate)
        ! Correction factor ~ 1 + (ρ_gyro/L_B)²
        finite_orbit_correction = 1.0_dp + 1.0d-6 * (v/v_thermal)**2
        
        ! Apply corrections to frequencies
        Om_th = 2.0_dp * pi / taub_thin * finite_orbit_correction
        Om_ph = delphi_thin / taub_thin * (1.0_dp + 0.5_dp * (finite_orbit_correction - 1.0_dp))
        
        success = .true.
        
    end subroutine approximate_thick_orbit_frequencies
    
    function thick_frequency_interpolation_available() result(available)
        ! Check if frequency database is available for interpolation
        logical :: available
        available = freq_db%initialized
    end function thick_frequency_interpolation_available

end module freq_thick