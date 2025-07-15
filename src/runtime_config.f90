module runtime_config
    ! Runtime configuration for NEO-RT thick orbit integration
    ! Controls whether to use thin orbit (NEO-RT) or thick orbit (POTATO) calculations
    
    implicit none
    
    private
    public :: get_use_thick_orbits, set_use_thick_orbits
    public :: init_runtime_config, get_config_from_namelist
    
    logical :: use_thick_orbits = .false.
    logical :: config_initialized = .false.
    
contains

    function get_use_thick_orbits() result(use_thick)
        ! Get the current thick orbit setting
        logical :: use_thick
        
        if (.not. config_initialized) then
            call init_runtime_config()
        end if
        
        use_thick = use_thick_orbits
    end function get_use_thick_orbits
    
    subroutine set_use_thick_orbits(use_thick)
        ! Set the thick orbit setting
        logical, intent(in) :: use_thick
        
        use_thick_orbits = use_thick
        config_initialized = .true.
    end subroutine set_use_thick_orbits
    
    subroutine init_runtime_config()
        ! Initialize runtime configuration from environment or defaults
        logical :: thick_orbits_enabled
        
        ! Try to read from namelist first
        call get_config_from_namelist(thick_orbits_enabled)
        
        ! Set the configuration
        use_thick_orbits = thick_orbits_enabled
        config_initialized = .true.
        
        ! Print configuration status
        if (use_thick_orbits) then
            print *, 'NEO-RT: Using thick orbit calculations (POTATO)'
        else
            print *, 'NEO-RT: Using thin orbit calculations (standard)'
        end if
        
    end subroutine init_runtime_config
    
    subroutine get_config_from_namelist(thick_orbits_enabled)
        ! Read configuration from namelist file
        logical, intent(out) :: thick_orbits_enabled
        
        ! Namelist variables
        logical :: use_thick_orbit_physics = .false.
        
        ! Namelist definition
        namelist /thick_orbit_config/ use_thick_orbit_physics
        
        ! Try to read from configuration file
        thick_orbits_enabled = .false.  ! Default to thin orbits
        
        ! TODO: Read from actual namelist file when available
        ! For now, check if build system enabled thick orbits
        thick_orbits_enabled = .false.
        
    end subroutine get_config_from_namelist
    
end module runtime_config