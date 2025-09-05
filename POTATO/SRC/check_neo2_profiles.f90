program check_neo2_profiles
    use neo2_reader, only: neo2_data_t, profiles_t, read_neo2_data, polyfit, polyval
    use fortplot

    implicit none

    character(len=100), parameter :: neo2_file = 'neo2_out.h5'

    call check_neo2_profiles_main

    contains

    subroutine check_neo2_profiles_main
        type(neo2_data_t) :: neo2_data
        type(profiles_t) :: profiles
        integer, parameter :: order = 9
        integer, parameter :: species_idx = 2  ! ions
        integer :: ns

        call read_neo2_data(neo2_file, neo2_data)

        ns = size(neo2_data%spol)
        allocate(profiles%n(0:order))
        call polyfit(neo2_data%spol, neo2_data%n_spec(species_idx,:), order, profiles%n)
        allocate(profiles%T(0:order))
        call polyfit(neo2_data%spol, neo2_data%T_spec(species_idx,:), order, profiles%T)
        
        call plot_fits(neo2_data, profiles, species_idx)
    end subroutine check_neo2_profiles_main
    
    subroutine plot_fits(neo2_data, profiles, species_idx)
        use fortplot
        
        type(neo2_data_t), intent(in) :: neo2_data
        type(profiles_t), intent(in) :: profiles
        integer, intent(in) :: species_idx
        
        real(8), allocatable :: spol_fine(:), dens_fit(:), temp_fit(:)
        integer :: n_fine, i
        
        ! Create fine grid for smooth curves
        n_fine = 200
        allocate(spol_fine(n_fine))
        allocate(dens_fit(n_fine))
        allocate(temp_fit(n_fine))
        
        do i = 1, n_fine
            spol_fine(i) = real(i-1, 8) / real(n_fine-1, 8)
        enddo
        
        ! Evaluate fits
        call polyval(spol_fine, profiles%n, dens_fit)
        call polyval(spol_fine, profiles%T, temp_fit)
        
        ! Plot density
        call figure()
        call plot(neo2_data%spol, neo2_data%n_spec(species_idx,:), 'NEO-2 data', 'ro')
        call plot(spol_fine, dens_fit, 'Polynomial fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('Density [1/cm³]')
        call title('Ion density profile')
        call savefig('density_fit.png')
        
        ! Plot temperature
        call figure()
        call plot(neo2_data%spol, neo2_data%T_spec(species_idx,:), 'NEO-2 data', 'ro')
        call plot(spol_fine, temp_fit, 'Polynomial fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('Temperature [eV]')
        call title('Ion temperature profile')
        call savefig('temperature_fit.png')
        
        print *, 'Plots saved to density_fit.png and temperature_fit.png'
        
        deallocate(spol_fine, dens_fit, temp_fit)
    end subroutine plot_fits

end program check_neo2_profiles