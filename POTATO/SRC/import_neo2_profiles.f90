program import_neo2_profiles
    use neo2_reader, only: neo2_data_t, profiles_t, read_neo2_data, polyfit, polyval, polyint
    use fortplot

    implicit none

    character(len=100), parameter :: neo2_file = 'neo2_out.h5'

    call import_neo2_profiles_main

    contains

    subroutine import_neo2_profiles_main
        use phielec_of_psi_mod, only: npolyphi
        
        type(neo2_data_t) :: neo2_data
        type(profiles_t) :: profiles
        real(8), allocatable :: Er_poly(:), T_in_eV(:)
        real(8), parameter :: erg_in_eV = 6.242d11  ! Conversion factor
        integer, parameter :: order = npolyphi  ! Use POTATO's standard order (9)
        integer, parameter :: order_Er = order + 1  ! One order higher for Er
        integer, parameter :: species_idx = 2  ! ions
        integer :: ns, iunit, i

        print *, 'Reading NEO-2 data from: ', trim(neo2_file)
        call read_neo2_data(neo2_file, neo2_data)

        ns = size(neo2_data%spol)
        
        ! Fit density with order 9
        allocate(profiles%n(0:order))
        call polyfit(neo2_data%spol, neo2_data%n_spec(species_idx,:), order, profiles%n)
        
        ! Convert temperature from erg to eV, then fit
        allocate(T_in_eV(ns))
        T_in_eV = neo2_data%T_spec(species_idx,:) * erg_in_eV
        allocate(profiles%T(0:order))
        call polyfit(neo2_data%spol, T_in_eV, order, profiles%T)
        
        ! Fit Er with order 10, then integrate to get phi_e
        allocate(Er_poly(0:order_Er))
        call polyfit(neo2_data%spol, neo2_data%Er, order_Er, Er_poly)
        
        ! Integrate Er to get phi_e (phi_e = -integral(Er * dspol))
        ! Need to negate Er before integrating
        Er_poly = -Er_poly
        allocate(profiles%phi_e(0:order_Er+1))  ! Result will be order 11 after integration
        
        call polyint(Er_poly, profiles%phi_e)
        
        ! Write profile_poly.in file in POTATO format
        open(newunit=iunit, file='profile_poly.in', status='replace', action='write')
        
        ! Line 1: dummy (not used by POTATO)
        write(iunit, *) '# Generated from NEO-2 data by import_neo2_profiles'
        
        ! Line 2: dummy (not used by POTATO)  
        write(iunit, *) '# Polynomial order:', order, ' Species:', species_idx
        
        ! Line 3: density polynomial coefficients (cm^-3 in CGS)
        write(iunit, *) (profiles%n(i), i=0,order)
        
        ! Line 4: dummy (originally for electron temperature)
        write(iunit, *) (profiles%T(i), i=0,order)  ! Could be electron temp if species_idx=1
        
        ! Line 5: ion temperature polynomial coefficients (eV)
        write(iunit, *) (profiles%T(i), i=0,order)
        
        ! Line 6: electrostatic potential polynomial coefficients (statV in CGS)
        ! Only use first order+1 coefficients (0:order) from phi_e
        write(iunit, *) (profiles%phi_e(i), i=0,order)
        
        close(iunit)
        
        print *, 'Successfully wrote profile_poly.in'
        print *, 'Polynomial order:', order
        print *, 'Coefficients for density, temperature, and phi_e saved'
        
        ! Create plots for verification
        call plot_fits(neo2_data, profiles, species_idx)
        call plot_Er_and_phi(neo2_data, Er_poly, profiles%phi_e)
    end subroutine import_neo2_profiles_main

    subroutine plot_fits(neo2_data, profiles, species_idx)
        use fortplot
        
        type(neo2_data_t), intent(in) :: neo2_data
        type(profiles_t), intent(in) :: profiles
        integer, intent(in) :: species_idx
        
        real(8), allocatable :: spol_fine(:), dens_fit(:), temp_fit(:), T_data_eV(:)
        real(8), parameter :: erg_in_eV = 6.242d11
        integer :: n_fine, i
        
        ! Create fine grid for smooth curves
        n_fine = 200
        allocate(spol_fine(n_fine))
        allocate(dens_fit(n_fine))
        allocate(temp_fit(n_fine))
        allocate(T_data_eV(size(neo2_data%spol)))
        
        do i = 1, n_fine
            spol_fine(i) = real(i-1, 8) / real(n_fine-1, 8)
        enddo
        
        ! Evaluate fits
        call polyval(spol_fine, profiles%n, dens_fit)
        call polyval(spol_fine, profiles%T, temp_fit)
        
        ! Convert temperature data to eV for plotting
        T_data_eV = neo2_data%T_spec(species_idx,:) * erg_in_eV
        
        ! Plot density (CGS units)
        call figure()
        call plot(neo2_data%spol, neo2_data%n_spec(species_idx,:), 'NEO-2 data', 'ro')
        call plot(spol_fine, dens_fit, 'Polynomial fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('Density [cm⁻³]')
        call title('Ion density profile')
        call savefig('density_fit.png')
        
        ! Plot temperature in eV
        call figure()
        call plot(neo2_data%spol, T_data_eV, 'NEO-2 data', 'ro')
        call plot(spol_fine, temp_fit, 'Polynomial fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('Temperature [eV]')
        call title('Ion temperature profile')
        call savefig('temperature_fit.png')
        
        print *, 'Plots saved to density_fit.png and temperature_fit.png'
        
        deallocate(spol_fine, dens_fit, temp_fit)
    end subroutine plot_fits
    
    subroutine plot_Er_and_phi(neo2_data, Er_poly, phi_e_poly)
        use fortplot
        
        type(neo2_data_t), intent(in) :: neo2_data
        real(8), intent(in) :: Er_poly(0:), phi_e_poly(0:)
        
        real(8), allocatable :: spol_fine(:), Er_fit(:), phi_fit(:), phi_from_data(:)
        integer :: n_fine, i
        
        ! Create fine grid
        n_fine = 200
        allocate(spol_fine(n_fine))
        allocate(Er_fit(n_fine))
        allocate(phi_fit(n_fine))
        allocate(phi_from_data(size(neo2_data%spol)))
        
        do i = 1, n_fine
            spol_fine(i) = real(i-1, 8) / real(n_fine-1, 8)
        enddo
        
        ! Evaluate fits (remember Er_poly is already negated)
        call polyval(spol_fine, -Er_poly, Er_fit)  ! Negate back for plotting
        call polyval(spol_fine, phi_e_poly, phi_fit)
        
        ! Integrate Er data numerically for comparison
        phi_from_data(1) = 0.0d0
        do i = 2, size(neo2_data%spol)
            phi_from_data(i) = phi_from_data(i-1) - 0.5d0 * (neo2_data%Er(i) + neo2_data%Er(i-1)) &
                             * (neo2_data%spol(i) - neo2_data%spol(i-1))
        enddo
        
        ! Plot Er (CGS: statV/cm)
        call figure()
        call plot(neo2_data%spol, neo2_data%Er, 'NEO-2 data', 'ro')
        call plot(spol_fine, Er_fit, 'Polynomial fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('Er [statV/cm]')
        call title('Radial electric field')
        call savefig('Er_fit.png')
        
        ! Plot phi_e (CGS: statV = erg/esu)
        call figure()
        call plot(neo2_data%spol, phi_from_data, 'Integrated NEO-2 Er', 'ro')
        call plot(spol_fine, phi_fit, 'Polynomial from integrated fit', 'b-')
        call xlabel('Normalized poloidal flux')
        call ylabel('φₑ [statV]')
        call title('Electrostatic potential')
        call savefig('phi_e_fit.png')
        
        print *, 'Additional plots saved to Er_fit.png and phi_e_fit.png'
        
        deallocate(spol_fine, Er_fit, phi_fit, phi_from_data)
    end subroutine plot_Er_and_phi

end program import_neo2_profiles